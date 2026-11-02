//! Humidity-consistency train, plan C4: one validator, one default and one impossible-air rule at
//! every input boundary.
//!
//! A relative humidity is a PERCENT in `[0, 100]` on every Rust and app surface (a FRACTION in
//! `[0, 1]` only on the SI JSON wire). A sensor reading within 5 points of that range (S5b,
//! `HUMIDITY_SENSOR_TOLERANCE_PERCENT`) is the saturated or dry air it reads and flies exactly as
//! the end; a value beyond the band, a non-finite one, or air that cannot exist is REJECTED
//! with a field-named error. 0.49.0 clamped it instead on most
//! surfaces (150 flew as 100 %, -5 as 0 %), and every solve capped impossible air while
//! `atmosphere.density_altitude` refused it. The one validator is `RelativeHumidity` and the one
//! impossible-air rule is `StationAir::new`'s (the density-altitude service's own predicate).
//!
//! The card, hold-curve, corridor, bridge and C ABI surfaces range-check the humidity at their own
//! parse step and refuse impossible air through the solver's one rule; both refusals are pinned
//! here per surface, so a surface that stopped reaching either check is seen.
//!
//! The sensor-tolerance band (step S5b, owner decision 2026-11-01): weather is an input source,
//! never a gate, so a READING within 5 points of the range (a Kestrel at 101.7 % in fog, a
//! phone at -2 % in desert air) is the saturated or dry air it reads, on every surface, and
//! flies that end to the bit. Only a value beyond the band (`BAD_PERCENTS`) is refused.
//!
//! Public API only. The CLI surfaces are in `tests/humidity_cli.rs`, the WASM terminal in
//! `src/wasm_tests.rs`, and the source guard in `tests/humidity_guard.rs`.

use ballistics_engine::atmosphere::{AtmoSegment, StationAtmospherePresence};
use ballistics_engine::{
    calculate_zero_angle_with_conditions, run_monte_carlo, AtmosphericConditions, BallisticInputs,
    DragModel, MonteCarloParams, RelativeHumidity, StationAir, TrajectorySolver, WindConditions,
    DEFAULT_HUMIDITY_PERCENT,
};

/// Percents every boundary must refuse: beyond the 5-point sensor band on either side, by a
/// little and by a lot, and not finite. (Before S5b this held -5, -0.0001 and 100.0001, which are
/// now sensor readings: see [`SENSOR_READINGS`].)
const BAD_PERCENTS: [f64; 7] = [
    -5.000_000_1,
    150.0,
    -20.0,
    105.000_000_1,
    f64::NAN,
    f64::INFINITY,
    f64::NEG_INFINITY,
];

/// Sensor readings within 5 points of the range, each with the end it reads as: saturated air
/// above 100 and dry air below 0 (a negative zero included). Each must fly the end to the bit.
const SENSOR_READINGS: [(f64, f64); 7] = [
    (100.000_000_1, 100.0),
    (101.7, 100.0),
    (105.0, 100.0),
    (-0.0, 0.0),
    (-0.0001, 0.0),
    (-2.0, 0.0),
    (-5.0, 0.0),
];

/// Both ends of the range are humidities, and so is a negative zero.
const GOOD_PERCENTS: [f64; 3] = [0.0, 100.0, -0.0];

/// Air that cannot exist, as `(temperature C, pressure hPa, humidity %)`: saturated air at the
/// boiling point needs more vapour pressure than the air's total pressure, and at 95 C and
/// 500 hPa the saturation pressure alone is above the total.
const IMPOSSIBLE_AIR: [(f64, f64, f64); 2] = [(100.0, 1013.25, 100.0), (95.0, 500.0, 100.0)];

/// The `AtmosphereErrorReason::ImpossibleAir` text: the shared refusal, word for word.
const IMPOSSIBLE_AIR_REASON: &str = "no such air exists";

/// Each [`IMPOSSIBLE_AIR`] case at half the humidity: the same temperature and pressure, and air
/// that exists, so a surface that refuses the first and flies this one is refusing the AIR.
fn possible_counterpart((temperature, pressure, humidity): (f64, f64, f64)) -> (f64, f64, f64) {
    (temperature, pressure, humidity / 2.0)
}

/// Asserts `message` is the one impossible-air refusal, for the surface named by `label`.
fn assert_impossible_air_refusal(label: &str, message: &str) {
    assert!(
        message.contains(IMPOSSIBLE_AIR_REASON),
        "{label}: refused with StationAir's own reason: {message}"
    );
}

/// Whether `message` is `RelativeHumidity`'s own refusal: the range rule, or (for NaN and the
/// infinities) the finiteness rule, each stating the sensor band it enforces (R2), derived here
/// from the constant.
fn names_the_humidity_rule(message: &str) -> bool {
    let t = ballistics_engine::HUMIDITY_SENSOR_TOLERANCE_PERCENT;
    let band = format!(
        "percent from 0 to 100 (a sensor reading from {} to {} is read as 0 or 100)",
        -t,
        100.0 + t
    );
    message.contains(&format!("relative humidity must be a {band}"))
        || message.contains(&format!("relative humidity must be a finite {band}"))
}

/// A .308 168 gr G7 load, short enough to keep the suite fast.
fn load() -> BallisticInputs {
    BallisticInputs {
        bc_value: 0.243,
        bc_type: DragModel::G7,
        bullet_mass: 0.01134,
        muzzle_velocity: 823.0,
        bullet_diameter: 0.00782,
        bullet_length: 0.031,
        muzzle_angle: 0.002,
        ..BallisticInputs::default()
    }
}

fn atmosphere(temperature: f64, pressure: f64, humidity: f64) -> AtmosphericConditions {
    AtmosphericConditions {
        temperature,
        pressure,
        humidity,
        altitude: 0.0,
    }
}

/// A solver whose station temperature and pressure fly exactly as given.
fn solver(inputs: BallisticInputs, atmosphere: AtmosphericConditions) -> TrajectorySolver {
    let mut solver = TrajectorySolver::new_with_station_presence(
        inputs,
        WindConditions::default(),
        atmosphere,
        StationAtmospherePresence::BOTH,
    );
    solver.set_max_range(300.0);
    solver
}

/// `(use_rk4, use_adaptive_rk45)` for Euler, fixed RK4 and adaptive RK45.
const METHODS: [(bool, bool); 3] = [(false, false), (true, false), (true, true)];

fn with_method(method: (bool, bool)) -> BallisticInputs {
    BallisticInputs {
        use_rk4: method.0,
        use_adaptive_rk45: method.1,
        ..load()
    }
}

// --- the library solver -----------------------------------------------------------------------

#[test]
fn trajectory_solver_rejects_an_out_of_range_humidity_and_never_clamps_it() {
    for method in METHODS {
        for humidity in BAD_PERCENTS {
            let error = solver(with_method(method), atmosphere(15.0, 1013.25, humidity))
                .solve()
                .expect_err(&format!("{method:?}: RH {humidity} must be refused"));
            let message = error.to_string();
            assert!(
                message.contains("atmosphere.humidity"),
                "{method:?}: RH {humidity}: the error must name the field: {message}"
            );
        }
        for humidity in GOOD_PERCENTS {
            solver(with_method(method), atmosphere(15.0, 1013.25, humidity))
                .solve()
                .unwrap_or_else(|e| panic!("{method:?}: RH {humidity} is valid: {e}"));
        }
    }
}

#[test]
fn trajectory_solver_refuses_impossible_air_by_the_one_shared_rule() {
    for method in METHODS {
        for (temperature, pressure, humidity) in IMPOSSIBLE_AIR {
            let message = solver(
                with_method(method),
                atmosphere(temperature, pressure, humidity),
            )
            .solve()
            .expect_err(&format!(
                "{method:?}: {temperature} C, {pressure} hPa, {humidity} % cannot exist"
            ))
            .to_string();
            assert!(
                message.contains(IMPOSSIBLE_AIR_REASON) && message.contains("atmosphere."),
                "{method:?}: the solve must refuse it with StationAir's own reason: {message}"
            );
        }
    }
    // The same air one step drier exists, and flies.
    solver(load(), atmosphere(95.0, 1013.25, 50.0))
        .solve()
        .expect("95 C / 1013.25 hPa / 50 % exists");
    // S5b: the sensor band does not relax it. Saturated air at the boiling point read as a
    // sensor's 102 % is the 100 % that cannot exist.
    for reading in [101.7, 105.0] {
        let message = solver(load(), atmosphere(100.0, 1013.25, reading))
            .solve()
            .expect_err("a reading of air that cannot exist")
            .to_string();
        assert!(
            message.contains(IMPOSSIBLE_AIR_REASON),
            "{reading}: {message}"
        );
    }
}

#[test]
fn trajectory_solver_refuses_a_non_positive_resolved_station_pressure() {
    for pressure in [0.0, -1013.25] {
        let message = solver(load(), atmosphere(15.0, pressure, 50.0))
            .solve()
            .expect_err(&format!("{pressure} hPa is not air"))
            .to_string();
        assert!(
            message.contains("atmosphere.pressure"),
            "the error must name the field: {message}"
        );
    }
}

#[test]
fn trajectory_solver_validates_every_downrange_atmosphere_zone() {
    // (zone, the field its error must name). Zone 1 is the bad one, after a good zone 0.
    let bad_zones: [(AtmoSegment, &str); 7] = [
        (
            (15.0, 1013.25, 150.0, 600.0),
            "atmo_sock.zones[1].humidity_percent",
        ),
        (
            (15.0, 1013.25, -5.000_000_1, 600.0),
            "atmo_sock.zones[1].humidity_percent",
        ),
        (
            (15.0, 1013.25, f64::NAN, 600.0),
            "atmo_sock.zones[1].humidity_percent",
        ),
        (
            (100.0, 1013.25, 100.0, 600.0),
            "atmo_sock.zones[1].humidity_percent",
        ),
        ((15.0, 0.0, 50.0, 600.0), "atmo_sock.zones[1].pressure_hpa"),
        (
            (-300.0, 1013.25, 50.0, 600.0),
            "atmo_sock.zones[1].temperature_c",
        ),
        ((15.0, 1013.25, 50.0, f64::NAN), "atmo_sock.zones"),
    ];
    for (zone, field) in bad_zones {
        let mut solver = solver(load(), atmosphere(15.0, 1013.25, 50.0));
        solver.set_atmo_segments(vec![(15.0, 1013.25, 50.0, 100.0), zone]);
        let message = solver
            .solve()
            .expect_err(&format!("zone {zone:?} must be refused"))
            .to_string();
        assert!(message.contains(field), "zone {zone:?}: {message}");
    }
    // Both ends of the humidity range are valid zones.
    let mut solver = solver(load(), atmosphere(15.0, 1013.25, 50.0));
    solver.set_atmo_segments(vec![
        (15.0, 1013.25, 0.0, 100.0),
        (15.0, 1013.25, 100.0, 600.0),
    ]);
    solver.solve().expect("0 % and 100 % zones are valid");
}

#[test]
fn zero_and_monte_carlo_inherit_the_refusal() {
    let zero = calculate_zero_angle_with_conditions(
        load(),
        100.0,
        0.0,
        WindConditions::default(),
        atmosphere(15.0, 1013.25, 150.0),
    );
    assert!(zero.is_err(), "the zero solve must refuse RH 150: {zero:?}");

    let params = MonteCarloParams {
        num_simulations: 3,
        target_distance: Some(100.0),
        ..MonteCarloParams::default()
    };
    for humidity in [150.0, -5.000_000_1] {
        let result = run_monte_carlo(
            BallisticInputs {
                humidity_percent: humidity,
                ..load()
            },
            params.clone(),
        );
        assert!(
            result.is_err(),
            "the Monte Carlo must refuse RH {humidity}, not fly it clamped"
        );
    }
}

#[test]
fn the_bc_estimator_names_the_air_it_refuses() {
    use ballistics_engine::{estimate_bc_fit, estimate_bc_fit_with_station_presence, BcFitMode};
    // Bore-referenced drops a G7 load at 800 m/s fits.
    let points = [(100.0, 0.08), (200.0, 0.34), (300.0, 0.82)];
    let fit = |atmosphere: AtmosphericConditions| {
        estimate_bc_fit_with_station_presence(
            800.0,
            0.0109,
            0.00782,
            &points,
            DragModel::G7,
            BcFitMode::Drop,
            atmosphere,
            None,
            0.05,
            StationAtmospherePresence::BOTH,
        )
    };
    // 0.49.0 and the first S5 round swallowed every trial's refusal and reported "Unable to
    // estimate BC from provided data", which blames the data for the air.
    for air in IMPOSSIBLE_AIR {
        let (temperature, pressure, humidity) = air;
        let message = fit(atmosphere(temperature, pressure, humidity))
            .expect_err(&format!("{air:?}: no such air exists"))
            .to_string();
        assert_impossible_air_refusal(&format!("estimate_bc_fit {air:?}"), &message);
        let (temperature, pressure, humidity) = possible_counterpart(air);
        fit(atmosphere(temperature, pressure, humidity))
            .expect("the half-humidity air exists, and the data fits");
    }
    for humidity in BAD_PERCENTS {
        for (name, result) in [
            ("with presence", fit(atmosphere(15.0, 1013.25, humidity))),
            (
                "presence-blind",
                estimate_bc_fit(
                    800.0,
                    0.0109,
                    0.00782,
                    &points,
                    DragModel::G7,
                    BcFitMode::Drop,
                    atmosphere(15.0, 1013.25, humidity),
                    None,
                    0.05,
                ),
            ),
        ] {
            let message = result
                .expect_err(&format!("{name}: RH {humidity} must be refused"))
                .to_string();
            assert!(
                message.contains("atmosphere.humidity") && names_the_humidity_rule(&message),
                "{name}: RH {humidity}: the one validator's field-named error: {message}"
            );
        }
    }
}

// --- one impossible-air rule ------------------------------------------------------------------

#[test]
fn impossible_air_is_refused_identically_by_the_station_the_solve_solve_v1_and_density_altitude() {
    for (temperature, pressure, humidity) in IMPOSSIBLE_AIR {
        let station = StationAir::new(
            0.0,
            temperature,
            pressure,
            RelativeHumidity::try_percent(humidity).unwrap(),
        )
        .expect_err("StationAir::new refuses it");
        assert!(station.to_string().contains(IMPOSSIBLE_AIR_REASON));

        let wire_atmosphere = serde_json::json!({
            "temperature_k": temperature + 273.15,
            "pressure_pa": pressure * 100.0,
            "relative_humidity": humidity / 100.0,
        });
        let solve = ballistics_engine::solve_json::decode_solve_request_v1(
            &serde_json::json!({
                "schema_version": 1,
                "projectile": {"mass_kg": 0.01134, "diameter_m": 0.00782, "length_m": 0.031,
                               "drag_model": "G7", "ballistic_coefficient": 0.243},
                "rifle": {"muzzle_velocity_mps": 823.0},
                "shot": {"max_range_m": 300.0},
                "atmosphere": wire_atmosphere,
                "wind": {}, "solver": {}, "effects": {}, "sampling": {}
            })
            .to_string(),
        )
        .map_err(|e| serde_json::to_value(&e).unwrap())
        .and_then(|request| {
            ballistics_engine::solve_v1(request).map_err(|e| serde_json::to_value(&e).unwrap())
        });
        let solve_error = solve.expect_err("solve v1 must refuse impossible air");
        assert!(
            solve_error["error"]["message"]
                .as_str()
                .unwrap_or_default()
                .contains(IMPOSSIBLE_AIR_REASON),
            "solve v1 refuses it for the shared reason: {solve_error}"
        );

        let density_altitude =
            ballistics_engine::atmosphere_service::decode_density_altitude_request_v1(
                &serde_json::json!({"atmosphere": wire_atmosphere}),
            )
            .and_then(|request| {
                ballistics_engine::atmosphere_service::density_altitude_v1(&request)
            })
            .map_err(|e| serde_json::to_value(&e).unwrap());
        let density_altitude_error = density_altitude.expect_err("density_altitude refuses it too");
        // R2: refused identically, envelope and all: one code, one path, one message.
        assert_eq!(
            solve_error["error"]["code"], "invalid_value",
            "{solve_error}"
        );
        assert_eq!(
            solve_error["error"]["path"], "$.atmosphere.relative_humidity",
            "{solve_error}"
        );
        assert_eq!(
            solve_error["error"], density_altitude_error["error"],
            "the solve and density_altitude refuse the same air the same way"
        );
    }
}

#[test]
fn solve_v1_keeps_refusing_a_humidity_fraction_outside_zero_to_one() {
    // Beyond the 0.05 sensor band (S5b moved 1.0001, a saturated-air reading, to the band test).
    for fraction in [1.5, -0.1, 50.0, 1.050_000_1, -0.050_000_1] {
        let decoded = ballistics_engine::solve_json::decode_solve_request_v1(
            &serde_json::json!({
                "schema_version": 1,
                "projectile": {"mass_kg": 0.01134, "diameter_m": 0.00782, "length_m": 0.031,
                               "drag_model": "G7", "ballistic_coefficient": 0.243},
                "rifle": {"muzzle_velocity_mps": 823.0},
                "shot": {"max_range_m": 300.0},
                "atmosphere": {"relative_humidity": fraction},
                "wind": {}, "solver": {}, "effects": {}, "sampling": {}
            })
            .to_string(),
        );
        let refused = match decoded {
            Err(error) => format!("{error:?}"),
            Ok(request) => match ballistics_engine::solve_v1(request) {
                Err(error) => serde_json::to_string(&error).unwrap(),
                Ok(_) => panic!("relative_humidity {fraction} must be refused"),
            },
        };
        assert!(
            refused.contains("relative_humidity"),
            "{fraction}: the error names the wire field: {refused}"
        );
    }
}

// --- card_service, HoldCurveLoad, CorridorLoad --------------------------------------------------

/// The card fixture on the wire: `humidity` omitted.
fn card_request_json() -> serde_json::Value {
    serde_json::json!({
        "units": "imperial",
        "muzzle_velocity": 2600.0,
        "ballistic_coefficient": 0.243,
        "mass": 175.0,
        "diameter": 0.308,
        "drag_model": "g7",
        "zero_distance": 100.0,
        "start": 100.0,
        "end": 300.0,
        "step": 100.0
    })
}

fn card_request(humidity: Option<f64>) -> ballistics_engine::card_service::CardRequestV1 {
    let mut request = card_request_json();
    if let Some(humidity) = humidity.filter(|h| h.is_finite()) {
        request["humidity"] = serde_json::json!(humidity);
    }
    let mut request: ballistics_engine::card_service::CardRequestV1 =
        serde_json::from_value(request).expect("the fixture decodes");
    if let Some(humidity) = humidity.filter(|h| !h.is_finite()) {
        request.humidity = humidity;
    }
    request
}

#[test]
fn card_service_rejects_an_out_of_range_humidity() {
    use ballistics_engine::card_service::{come_ups_v1, range_table_v1, CardServiceError};
    for humidity in BAD_PERCENTS {
        for (name, result) in [
            (
                "come_ups",
                come_ups_v1(&card_request(Some(humidity))).map(|_| ()),
            ),
            (
                "range_table",
                range_table_v1(&card_request(Some(humidity))).map(|_| ()),
            ),
        ] {
            match result {
                Err(CardServiceError::InvalidRequest(message)) => assert!(
                    message.contains("humidity"),
                    "{name} RH {humidity}: the error names the field: {message}"
                ),
                other => panic!("{name} RH {humidity} must be an invalid request: {other:?}"),
            }
        }
    }
    for humidity in GOOD_PERCENTS {
        come_ups_v1(&card_request(Some(humidity))).expect("a valid humidity");
    }
}

#[test]
fn card_service_omitted_humidity_is_the_engine_default_bit_for_bit() {
    use ballistics_engine::card_service::come_ups_v1;
    let omitted = serde_json::to_string(&come_ups_v1(&card_request(None)).unwrap()).unwrap();
    let typed =
        serde_json::to_string(&come_ups_v1(&card_request(Some(DEFAULT_HUMIDITY_PERCENT))).unwrap())
            .unwrap();
    assert_eq!(omitted, typed);
    assert_eq!(card_request(None).humidity, DEFAULT_HUMIDITY_PERCENT);
}

fn hold_load(humidity: f64) -> ballistics_engine::hold_curve::HoldCurveLoad {
    ballistics_engine::hold_curve::HoldCurveLoad {
        velocity_mps: 800.0,
        bc: 0.223,
        mass_kg: 0.0109,
        diameter_m: 0.00782,
        drag_model: DragModel::G7,
        sight_height_m: 0.045,
        zero_distance_m: 100.0,
        temperature_c: 15.0,
        pressure_hpa: 1013.25,
        humidity,
        altitude_m: 0.0,
        wind_speed_mps: 0.0,
        wind_direction_deg: 90.0,
        station_presence: StationAtmospherePresence::NEITHER,
    }
}

#[test]
fn hold_curve_load_rejects_an_out_of_range_humidity_percent() {
    use ballistics_engine::hold_curve::HoldCurve;
    for humidity in BAD_PERCENTS {
        let message = HoldCurve::solve(&hold_load(humidity), 300.0)
            .err()
            .unwrap_or_else(|| panic!("RH {humidity} must be refused"))
            .to_string();
        assert!(
            message.contains("humidity") && names_the_humidity_rule(&message),
            "RH {humidity}: one validator's field-named error: {message}"
        );
    }
    for humidity in GOOD_PERCENTS {
        HoldCurve::solve(&hold_load(humidity), 300.0).expect("a valid humidity");
    }
}

fn corridor_request(humidity_pct: f64) -> ballistics_engine::wind_scenarios::RobustHoldRequest {
    use ballistics_engine::wind_scenarios::{
        parse_wind_scenario_set, CorridorLoad, RobustHoldRequest,
    };
    RobustHoldRequest {
        scenarios: parse_wind_scenario_set(
            r#"{"version":1,"scenarios":[{"name":"a","segments":["5:90:1000"]}]}"#,
            ballistics_engine::cli_api::UnitSystem::Imperial,
        )
        .expect("the fixture parses"),
        ranges_m: vec![200.0],
        target: None,
        load: CorridorLoad {
            muzzle_velocity_mps: 823.0,
            bc: 0.243,
            drag_model: DragModel::G7,
            mass_kg: 0.01134,
            diameter_m: 0.00782,
            bullet_length_m: 0.031,
            zero_distance_m: 100.0,
            sight_height_m: 0.038,
            temperature_c: 15.0,
            pressure_hpa: 1013.25,
            humidity_pct,
            altitude_m: 0.0,
            station_presence: StationAtmospherePresence::NEITHER,
        },
    }
}

#[test]
fn corridor_load_enforces_its_documented_zero_to_one_hundred_percent() {
    use ballistics_engine::wind_scenarios::solve_robust_hold;
    for humidity in BAD_PERCENTS {
        let message = solve_robust_hold(&corridor_request(humidity))
            .expect_err(&format!("RH {humidity} must be refused"))
            .to_string();
        assert!(
            message.contains("humidity_pct") && names_the_humidity_rule(&message),
            "RH {humidity}: the field-named error: {message}"
        );
    }
    for humidity in GOOD_PERCENTS {
        solve_robust_hold(&corridor_request(humidity)).expect("a valid humidity");
    }
}

/// The card fixture on the wire at the given station air, in degrees C and hPa.
fn card_request_json_in_air(
    (temperature, pressure, humidity): (f64, f64, f64),
) -> serde_json::Value {
    let mut request = card_request_json();
    request["temperature"] = serde_json::json!(temperature);
    request["temperature_unit"] = serde_json::json!("celsius");
    request["pressure"] = serde_json::json!(pressure);
    request["pressure_unit"] = serde_json::json!("hpa");
    request["humidity"] = serde_json::json!(humidity);
    request
}

fn card_request_in_air(air: (f64, f64, f64)) -> ballistics_engine::card_service::CardRequestV1 {
    serde_json::from_value(card_request_json_in_air(air)).expect("the fixture decodes")
}

#[test]
fn card_service_refuses_impossible_air_by_the_one_shared_rule() {
    use ballistics_engine::card_service::{come_ups_v1, range_table_v1};
    for air in IMPOSSIBLE_AIR {
        let request = card_request_in_air(air);
        for (name, result) in [
            ("come_ups", come_ups_v1(&request).map(|_| ())),
            ("range_table", range_table_v1(&request).map(|_| ())),
        ] {
            let message = result
                .expect_err(&format!("{name} {air:?}: no such air exists"))
                .to_string();
            assert_impossible_air_refusal(&format!("{name} {air:?}"), &message);
        }
        let possible = card_request_in_air(possible_counterpart(air));
        come_ups_v1(&possible).expect("the half-humidity air exists, and flies");
        range_table_v1(&possible).expect("the half-humidity air exists, and flies");
    }
}

fn hold_load_in_air(
    (temperature, pressure, humidity): (f64, f64, f64),
) -> ballistics_engine::hold_curve::HoldCurveLoad {
    ballistics_engine::hold_curve::HoldCurveLoad {
        temperature_c: temperature,
        pressure_hpa: pressure,
        station_presence: StationAtmospherePresence::BOTH,
        ..hold_load(humidity)
    }
}

#[test]
fn hold_curve_load_refuses_impossible_air_by_the_one_shared_rule() {
    use ballistics_engine::hold_curve::HoldCurve;
    for air in IMPOSSIBLE_AIR {
        let message = HoldCurve::solve(&hold_load_in_air(air), 300.0)
            .err()
            .unwrap_or_else(|| panic!("{air:?}: no such air exists"))
            .to_string();
        assert_impossible_air_refusal(&format!("HoldCurveLoad {air:?}"), &message);
        HoldCurve::solve(&hold_load_in_air(possible_counterpart(air)), 300.0)
            .expect("the half-humidity air exists, and flies");
    }
}

fn corridor_request_in_air(
    (temperature, pressure, humidity): (f64, f64, f64),
) -> ballistics_engine::wind_scenarios::RobustHoldRequest {
    let mut request = corridor_request(humidity);
    request.load.temperature_c = temperature;
    request.load.pressure_hpa = pressure;
    request.load.station_presence = StationAtmospherePresence::BOTH;
    request
}

#[test]
fn corridor_load_refuses_impossible_air_by_the_one_shared_rule() {
    use ballistics_engine::wind_scenarios::solve_robust_hold;
    for air in IMPOSSIBLE_AIR {
        let message = solve_robust_hold(&corridor_request_in_air(air))
            .expect_err(&format!("CorridorLoad {air:?}: no such air exists"))
            .to_string();
        assert_impossible_air_refusal(&format!("CorridorLoad {air:?}"), &message);
        solve_robust_hold(&corridor_request_in_air(possible_counterpart(air)))
            .expect("the half-humidity air exists, and flies");
    }
}

// --- the bridge -------------------------------------------------------------------------------

#[cfg(feature = "bridge")]
fn bridge(command: &str, request: serde_json::Value) -> serde_json::Value {
    let envelope = serde_json::json!({"api_version": 1, "command": command, "request": request});
    serde_json::from_str(&ballistics_engine::bridge::bridge_call(
        &envelope.to_string(),
    ))
    .unwrap()
}

#[cfg(feature = "bridge")]
#[test]
fn bridge_card_commands_reject_an_out_of_range_humidity() {
    for command in ["card.come_ups", "card.range_table", "card.wind"] {
        for humidity in [150.0, -5.000_000_1] {
            let mut request = card_request_json();
            request["humidity"] = serde_json::json!(humidity);
            if command == "card.wind" {
                request["wind_speeds"] = serde_json::json!([5.0, 10.0]);
            }
            let out = bridge(command, request);
            assert_eq!(out["ok"], false, "{command} RH {humidity} must fail: {out}");
            let message = out["error"]["message"].as_str().unwrap_or_default();
            assert!(message.contains("humidity"), "{command}: {message}");
        }
    }
}

#[cfg(feature = "bridge")]
#[test]
fn bridge_card_commands_refuse_impossible_air_by_the_one_shared_rule() {
    for command in ["card.come_ups", "card.range_table", "card.wind"] {
        for air in IMPOSSIBLE_AIR {
            let request = |air| {
                let mut request = card_request_json_in_air(air);
                if command == "card.wind" {
                    request["wind_speeds"] = serde_json::json!([5.0, 10.0]);
                }
                request
            };
            let out = bridge(command, request(air));
            assert_eq!(out["ok"], false, "{command} {air:?} must fail: {out}");
            assert_impossible_air_refusal(
                &format!("{command} {air:?}"),
                out["error"]["message"].as_str().unwrap_or_default(),
            );
            let out = bridge(command, request(possible_counterpart(air)));
            assert_eq!(
                out["ok"], true,
                "{command}: the half-humidity air flies: {out}"
            );
        }
    }
}

#[cfg(feature = "bridge")]
#[test]
fn bridge_true_commands_reject_through_the_shared_validator() {
    let model = |humidity: f64| {
        serde_json::json!({
            "muzzle_velocity_fps": 2700.0, "ballistic_coefficient": 0.243,
            "drag_model": "g7", "mass_gr": 168.0, "diameter_in": 0.308,
            "zero_distance_yd": 100.0, "sight_height_in": 2.0,
            "temperature_f": 59.0, "pressure_inhg": 29.92,
            "humidity_pct": humidity, "altitude_ft": 0.0
        })
    };
    for humidity in [150.0, -5.000_000_1] {
        let out = bridge(
            "true.dsf",
            serde_json::json!({
                "model": model(humidity), "range_yd": 600.0,
                "observed_drop": 4.8, "drop_unit": "mil"
            }),
        );
        assert_eq!(out["ok"], false, "true.dsf RH {humidity} must fail: {out}");
        let message = out["error"]["message"].as_str().unwrap_or_default();
        // R2: the refusal states the band it enforces.
        assert!(
            message.contains("humidity_pct = ") && names_the_humidity_rule(message),
            "true.dsf: the one validator's error: {message}"
        );
    }
}

#[cfg(feature = "bridge")]
#[test]
fn bridge_profile_validate_flags_an_out_of_range_stored_humidity() {
    let profile = |humidity: f64| {
        serde_json::json!({
            "name": "humid", "velocity": 2700.0, "bc": 0.243, "mass": 175.0,
            "diameter": 0.308, "drag_model": "G7", "humidity": humidity
        })
    };
    for humidity in [150.0, -5.000_000_1] {
        let out = bridge("profile.validate", profile(humidity));
        assert_eq!(out["ok"], true, "{out}");
        assert_eq!(out["result"]["valid"], false, "RH {humidity}: {out}");
        let warnings = out["result"]["warnings"].to_string();
        assert!(warnings.contains("humidity"), "{warnings}");
    }
    let out = bridge("profile.validate", profile(100.0));
    assert_eq!(out["result"]["valid"], true, "100 % is a humidity: {out}");
}

// --- the C ABI --------------------------------------------------------------------------------

#[cfg(feature = "ffi")]
mod ffi {
    use super::*;
    use ballistics_engine::ffi::{
        ballistics_calculate_trajectory, ballistics_calculate_zero_angle,
        ballistics_free_monte_carlo_results, ballistics_free_trajectory_result,
        ballistics_monte_carlo, FFIAtmosphericConditions, FFIBallisticInputs, FFIMonteCarloParams,
    };

    fn inputs() -> FFIBallisticInputs {
        FFIBallisticInputs {
            muzzle_velocity: 823.0,
            muzzle_angle: 0.002,
            bc_value: 0.475,
            bullet_mass: 0.01088,
            bullet_diameter: 0.00782,
            bc_type: 0,
            sight_height: 0.05,
            target_distance: 300.0,
            temperature: 15.0,
            twist_rate: BallisticInputs::default().twist_rate,
            is_twist_right: 1,
            shooting_angle: 0.0,
            altitude: 0.0,
            latitude: f64::NAN,
            azimuth_angle: 0.0,
            use_rk4: 1,
            use_adaptive_rk45: 0,
            enable_wind_shear: 0,
            enable_trajectory_sampling: 0,
            sample_interval: 0.0,
            enable_pitch_damping: 0,
            enable_precession_nutation: 0,
            enable_spin_drift: 0,
            enable_magnus: 0,
            enable_coriolis: 0,
            shot_azimuth: 0.0,
            cant_angle: 0.0,
            zero_poi_vertical: 0.0,
            zero_poi_horizontal: 0.0,
            sight_offset_lateral: 0.0,
        }
    }

    fn air(humidity: f64) -> FFIAtmosphericConditions {
        FFIAtmosphericConditions {
            temperature: 15.0,
            pressure: 1013.25,
            humidity,
            altitude: 0.0,
        }
    }

    /// The impact velocity, or `None` for a null result.
    fn trajectory(atmosphere: Option<&FFIAtmosphericConditions>) -> Option<u64> {
        let inputs = inputs();
        let atmosphere = atmosphere.map_or(std::ptr::null(), |a| a as *const _);
        unsafe {
            let result =
                ballistics_calculate_trajectory(&inputs, std::ptr::null(), atmosphere, 300.0, 1.0);
            if result.is_null() {
                return None;
            }
            let bits = (*result).impact_velocity.to_bits();
            ballistics_free_trajectory_result(result);
            Some(bits)
        }
    }

    fn monte_carlo_succeeds(atmosphere: &FFIAtmosphericConditions) -> bool {
        let inputs = inputs();
        let params = FFIMonteCarloParams {
            num_simulations: 3,
            velocity_std_dev: 0.0,
            angle_std_dev: 0.0,
            bc_std_dev: 0.0,
            wind_speed_std_dev: 0.0,
            target_distance: 100.0,
            base_wind_speed: 0.0,
            base_wind_direction: 0.0,
            azimuth_std_dev: 0.0,
        };
        unsafe {
            let results = ballistics_monte_carlo(&inputs, atmosphere, &params);
            let ok = !results.is_null();
            ballistics_free_monte_carlo_results(results);
            ok
        }
    }

    #[test]
    fn the_c_abi_single_solve_returns_null_for_a_humidity_it_cannot_fly() {
        for humidity in BAD_PERCENTS {
            assert_eq!(
                trajectory(Some(&air(humidity))),
                None,
                "RH {humidity}: null, not a clamped flight"
            );
            let inputs = inputs();
            let zero = unsafe {
                ballistics_calculate_zero_angle(&inputs, std::ptr::null(), &air(humidity), 100.0)
            };
            assert!(
                zero.is_nan(),
                "RH {humidity}: the zero angle is NaN, got {zero}"
            );
        }
        for humidity in GOOD_PERCENTS {
            assert!(trajectory(Some(&air(humidity))).is_some(), "RH {humidity}");
        }
    }

    #[test]
    fn the_c_abi_monte_carlo_returns_null_for_a_humidity_it_cannot_fly() {
        for humidity in BAD_PERCENTS {
            assert!(
                !monte_carlo_succeeds(&air(humidity)),
                "RH {humidity}: null, not a clamped Monte Carlo"
            );
        }
        for humidity in GOOD_PERCENTS {
            assert!(monte_carlo_succeeds(&air(humidity)), "RH {humidity}");
        }
    }

    fn air_at((temperature, pressure, humidity): (f64, f64, f64)) -> FFIAtmosphericConditions {
        FFIAtmosphericConditions {
            temperature,
            pressure,
            humidity,
            altitude: 0.0,
        }
    }

    #[test]
    fn the_c_abi_refuses_impossible_air_by_the_one_shared_rule() {
        for air in IMPOSSIBLE_AIR {
            assert_eq!(
                trajectory(Some(&air_at(air))),
                None,
                "{air:?}: the single solve is null, not a capped flight"
            );
            let inputs = inputs();
            let zero = unsafe {
                ballistics_calculate_zero_angle(&inputs, std::ptr::null(), &air_at(air), 100.0)
            };
            assert!(zero.is_nan(), "{air:?}: the zero angle is NaN, got {zero}");
            assert!(
                !monte_carlo_succeeds(&air_at(air)),
                "{air:?}: the Monte Carlo is null, not a capped run"
            );
            let possible = air_at(possible_counterpart(air));
            assert!(trajectory(Some(&possible)).is_some(), "{air:?} at half RH");
            assert!(monte_carlo_succeeds(&possible), "{air:?} at half RH");
        }
    }

    /// S5b: the C ABI reads a sensor reading within 5 points of the range as the end it reads
    /// as: the single solve and the zero angle are the end's to the bit, and the Monte Carlo
    /// flies it.
    #[test]
    fn the_c_abi_flies_a_sensor_reading_as_the_end_it_reads_as() {
        let zero = |humidity: f64| {
            let inputs = inputs();
            unsafe {
                ballistics_calculate_zero_angle(&inputs, std::ptr::null(), &air(humidity), 100.0)
            }
            .to_bits()
        };
        for (reading, end) in SENSOR_READINGS {
            let flown = trajectory(Some(&air(reading)));
            assert!(flown.is_some(), "RH {reading} flies");
            assert_eq!(flown, trajectory(Some(&air(end))), "RH {reading}");
            assert_eq!(zero(reading), zero(end), "zero angle at RH {reading}");
            assert!(
                monte_carlo_succeeds(&air(reading)),
                "Monte Carlo at RH {reading}"
            );
        }
    }

    #[test]
    fn a_null_c_atmosphere_is_the_engine_default_humidity_bit_for_bit() {
        assert_eq!(
            trajectory(None),
            trajectory(Some(&air(DEFAULT_HUMIDITY_PERCENT))),
            "a null atmosphere flies 15 C, 1013.25 hPa and the one default humidity"
        );
    }
}

// --- .a7p import --------------------------------------------------------------------------------

#[cfg(feature = "profile-import")]
mod a7p {
    use ballistics_engine::profile_import::{map_a7p_to_profile, parse_a7p, wrap_payload};

    fn enc_varint(mut v: u64, out: &mut Vec<u8>) {
        loop {
            let byte = (v & 0x7f) as u8;
            v >>= 7;
            if v == 0 {
                out.push(byte);
                return;
            }
            out.push(byte | 0x80);
        }
    }
    fn enc_i32(number: u32, value: i64, out: &mut Vec<u8>) {
        enc_varint(u64::from(number) << 3, out);
        enc_varint(value as u64, out);
    }
    fn enc_str(number: u32, s: &str, out: &mut Vec<u8>) {
        enc_varint((u64::from(number) << 3) | 2, out);
        enc_varint(s.len() as u64, out);
        out.extend_from_slice(s.as_bytes());
    }
    fn enc_bytes(number: u32, payload: &[u8], out: &mut Vec<u8>) {
        enc_varint((u64::from(number) << 3) | 2, out);
        enc_varint(payload.len() as u64, out);
        out.extend_from_slice(payload);
    }

    /// A G7 .a7p with field 17 (`c_zero_air_humidity`, an int32 percent) set to `humidity`, or
    /// absent when `None`.
    fn a7p(humidity: Option<i64>) -> Vec<u8> {
        let mut p = Vec::new();
        enc_str(1, "HUMID 6.5CM", &mut p);
        enc_i32(9, 50, &mut p);
        enc_i32(10, 800, &mut p);
        enc_i32(11, 8230, &mut p);
        enc_i32(15, 15, &mut p);
        enc_i32(16, 10132, &mut p);
        if let Some(humidity) = humidity {
            enc_i32(17, humidity, &mut p);
        }
        enc_i32(20, 264, &mut p);
        enc_i32(21, 1400, &mut p);
        enc_i32(22, 1360, &mut p);
        enc_i32(24, 1, &mut p);
        let mut packed = Vec::new();
        enc_varint(10_000, &mut packed);
        enc_bytes(26, &packed, &mut p);
        let mut row = Vec::new();
        enc_i32(1, 3260, &mut row);
        enc_i32(2, 8230, &mut row);
        enc_bytes(27, &row, &mut p);
        let mut payload = Vec::new();
        enc_bytes(1, &p, &mut payload);
        wrap_payload(&payload)
    }

    #[test]
    fn an_out_of_range_a7p_humidity_is_refused_at_import() {
        // Beyond the 5-point sensor band (S5b: -5 and 101 are readings, the next test's).
        for humidity in [150, -6, 106] {
            let doc = parse_a7p(&a7p(Some(humidity))).expect("the envelope parses");
            let error = map_a7p_to_profile(&doc, None, None)
                .err()
                .unwrap_or_else(|| panic!("RH {humidity} must be refused at import"));
            // R2: the refusal names the file's field and states the band it enforces.
            assert!(
                error.starts_with("c_zero_air_humidity = ")
                    && super::names_the_humidity_rule(&error),
                "RH {humidity}: the field-named error: {error}"
            );
        }
        for humidity in [0, 100, 37] {
            let doc = parse_a7p(&a7p(Some(humidity))).unwrap();
            let profile = map_a7p_to_profile(&doc, None, None).unwrap().profile;
            assert_eq!(profile.humidity, humidity as f64);
        }
    }

    /// S5b: a file's humidity within 5 points of the range is a sensor reading, imported as the
    /// end it reads as.
    #[test]
    fn an_a7p_humidity_within_five_points_of_the_range_imports_as_the_end() {
        for (humidity, end) in [
            (101, 100.0),
            (103, 100.0),
            (105, 100.0),
            (-1, 0.0),
            (-5, 0.0),
        ] {
            let doc = parse_a7p(&a7p(Some(humidity))).unwrap();
            let profile = map_a7p_to_profile(&doc, None, None)
                .unwrap_or_else(|e| panic!("RH {humidity} imports: {e}"))
                .profile;
            assert_eq!(
                profile.humidity.to_bits(),
                f64::to_bits(end),
                "RH {humidity}"
            );
        }
    }

    /// The deliberate exception to the one 50 % default: a proto3 int32 cannot tell an absent
    /// field from a 0, and the format writes 0 % by omitting the field (the a7p library, the
    /// accepted oracle, does the same), so an absent `c_zero_air_humidity` is 0 %.
    #[test]
    fn an_absent_a7p_humidity_is_zero_percent_by_the_formats_own_meaning() {
        let doc = parse_a7p(&a7p(None)).unwrap();
        let profile = map_a7p_to_profile(&doc, None, None).unwrap().profile;
        assert_eq!(profile.humidity.to_bits(), 0.0f64.to_bits());
    }
}

// --- defaults ---------------------------------------------------------------------------------

#[test]
fn every_library_default_is_the_one_default_humidity() {
    assert_eq!(DEFAULT_HUMIDITY_PERCENT, 50.0);
    assert_eq!(
        RelativeHumidity::DEFAULT.percent(),
        DEFAULT_HUMIDITY_PERCENT
    );
    assert_eq!(
        AtmosphericConditions::default().humidity.to_bits(),
        DEFAULT_HUMIDITY_PERCENT.to_bits()
    );
    assert_eq!(
        BallisticInputs::default().humidity_percent.to_bits(),
        DEFAULT_HUMIDITY_PERCENT.to_bits()
    );
    let profile: ballistics_engine::profile::ProfileData = serde_json::from_value(
        serde_json::json!({"name": "p", "velocity": 2700.0, "bc": 0.243, "mass": 175.0,
                           "diameter": 0.308, "drag_model": "G7"}),
    )
    .unwrap();
    assert_eq!(
        profile.humidity.to_bits(),
        DEFAULT_HUMIDITY_PERCENT.to_bits()
    );
}

// --- the online client --------------------------------------------------------------------------

/// The query string `ApiClient` sends for one request, captured by a one-shot loopback server
/// that answers 500 (the client's error is expected; only the request line matters).
#[cfg(feature = "online")]
fn captured_query(request: &ballistics_engine::api_client::TrajectoryRequest) -> String {
    use std::io::{BufRead, BufReader, Write};
    let listener = std::net::TcpListener::bind("127.0.0.1:0").expect("bind loopback");
    let port = listener.local_addr().unwrap().port();
    let server = std::thread::spawn(move || {
        let (stream, _) = listener.accept().expect("the client connects");
        let mut reader = BufReader::new(stream.try_clone().unwrap());
        let mut request_line = String::new();
        reader.read_line(&mut request_line).unwrap();
        let mut line = String::new();
        while reader.read_line(&mut line).unwrap() > 2 {
            line.clear();
        }
        let mut stream = stream;
        stream
            .write_all(b"HTTP/1.1 500 Internal Server Error\r\nContent-Length: 0\r\nConnection: close\r\n\r\n")
            .unwrap();
        request_line
    });
    let client =
        ballistics_engine::api_client::ApiClient::new(&format!("http://127.0.0.1:{port}"), 5);
    let _ = client.calculate_trajectory(request);
    server.join().unwrap()
}

#[cfg(feature = "online")]
#[test]
fn the_online_client_sends_the_humidity_at_full_precision() {
    for humidity in [37.375, 12.3456, 100.0, 0.0, 49.95] {
        let request = ballistics_engine::api_client::TrajectoryRequestBuilder::new()
            .bc_value(0.475)
            .bc_type("G1")
            .bullet_mass(10.9)
            .muzzle_velocity(823.0)
            .target_distance(500.0)
            .humidity(humidity)
            .build()
            .unwrap();
        let request_line = captured_query(&request);
        let sent = request_line
            .split(['?', '&', ' '])
            .find_map(|pair| pair.strip_prefix("humidity="))
            .unwrap_or_else(|| panic!("no humidity in {request_line:?}"));
        assert_eq!(
            sent.parse::<f64>().unwrap().to_bits(),
            humidity.to_bits(),
            "RH {humidity} must round-trip exactly, not to one decimal: sent {sent:?}"
        );
    }
}

// --- the sensor band, surface by surface (S5b) --------------------------------------------------

/// Every trajectory point of `solver`'s solve, as bits.
fn flight_bits(mut solver: TrajectorySolver) -> Vec<u64> {
    solver.set_max_range(300.0);
    let result = solver.solve().expect("a sensor reading flies");
    result
        .points
        .iter()
        .flat_map(|p| {
            [
                p.time,
                p.position.x,
                p.position.y,
                p.position.z,
                p.velocity_magnitude,
            ]
        })
        .map(f64::to_bits)
        .collect()
}

#[test]
fn trajectory_solver_flies_a_sensor_reading_as_the_end_it_reads_as() {
    for method in METHODS {
        for (reading, end) in SENSOR_READINGS {
            assert_eq!(
                flight_bits(solver(
                    with_method(method),
                    atmosphere(15.0, 1013.25, reading)
                )),
                flight_bits(solver(with_method(method), atmosphere(15.0, 1013.25, end))),
                "{method:?}: station RH {reading} must fly {end}"
            );
            // The same reading in a downrange zone.
            let zoned = |humidity: f64| {
                let mut solver = solver(with_method(method), atmosphere(15.0, 1013.25, 50.0));
                solver.set_atmo_segments(vec![
                    (15.0, 1013.25, 50.0, 100.0),
                    (25.0, 1000.0, humidity, 600.0),
                ]);
                flight_bits(solver)
            };
            assert_eq!(
                zoned(reading),
                zoned(end),
                "{method:?}: zone RH {reading} must fly {end}"
            );
        }
    }
}

#[test]
fn zero_monte_carlo_and_the_bc_estimator_read_a_sensor_reading_as_the_end() {
    use ballistics_engine::estimate_bc_fit;
    use ballistics_engine::monte_carlo::solve_trajectory_for_monte_carlo_with_station_presence;
    let points = [(100.0, 0.08), (200.0, 0.34), (300.0, 0.82)];
    let estimate = |humidity: f64| {
        let fit = estimate_bc_fit(
            800.0,
            0.0109,
            0.00782,
            &points,
            DragModel::G7,
            ballistics_engine::BcFitMode::Drop,
            atmosphere(15.0, 1013.25, humidity),
            None,
            0.05,
        )
        .expect("a sensor reading fits");
        format!("{fit:?}")
    };
    let zero = |humidity: f64| {
        calculate_zero_angle_with_conditions(
            load(),
            100.0,
            0.0,
            WindConditions::default(),
            atmosphere(15.0, 1013.25, humidity),
        )
        .expect("a sensor reading zeroes")
        .to_bits()
    };
    let lean = |humidity: f64| {
        let output = solve_trajectory_for_monte_carlo_with_station_presence(
            &BallisticInputs {
                humidity_percent: humidity,
                target_distance: 300.0,
                ..load()
            },
            StationAtmospherePresence::BOTH,
        )
        .expect("a sensor reading flies the lean Monte Carlo");
        format!("{output:?}")
    };
    let params = MonteCarloParams {
        num_simulations: 3,
        target_distance: Some(100.0),
        ..MonteCarloParams::default()
    };
    for (reading, end) in SENSOR_READINGS {
        assert_eq!(zero(reading), zero(end), "zero at RH {reading}");
        assert_eq!(lean(reading), lean(end), "lean Monte Carlo at RH {reading}");
        assert_eq!(
            estimate(reading),
            estimate(end),
            "estimate_bc_fit at RH {reading}"
        );
        run_monte_carlo(
            BallisticInputs {
                humidity_percent: reading,
                ..load()
            },
            params.clone(),
        )
        .unwrap_or_else(|e| panic!("the Monte Carlo flies RH {reading}: {e}"));
    }
}

/// A solve-json v1 request at the given atmosphere, solved, as the response's JSON or the
/// refusal's.
fn solve_v1_json(atmosphere: serde_json::Value) -> Result<String, String> {
    let request = ballistics_engine::solve_json::decode_solve_request_v1(
        &serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.01134, "diameter_m": 0.00782, "length_m": 0.031,
                           "drag_model": "G7", "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0},
            "shot": {"max_range_m": 300.0, "zero_distance_m": 100.0},
            "atmosphere": atmosphere,
            "wind": {}, "solver": {}, "effects": {}, "sampling": {}
        })
        .to_string(),
    )
    .map_err(|e| format!("{e:?}"))?;
    ballistics_engine::solve_v1(request)
        .map(|response| serde_json::to_string(&response).unwrap())
        .map_err(|e| serde_json::to_string(&e).unwrap())
}

/// Fractions within 0.05 of the SI wire's range, each with the end it reads as.
const SENSOR_FRACTIONS: [(f64, f64); 6] = [
    (1.000_000_1, 1.0),
    (1.03, 1.0),
    (1.05, 1.0),
    (-0.0, 0.0),
    (-0.02, 0.0),
    (-0.05, 0.0),
];

#[test]
fn solve_v1_and_density_altitude_read_a_sensor_fraction_as_the_end_it_reads_as() {
    use ballistics_engine::atmosphere_service::{
        decode_density_altitude_request_v1, density_altitude_v1,
    };
    let density_altitude = |fraction: f64| {
        let request = decode_density_altitude_request_v1(&serde_json::json!({
            "atmosphere": {"temperature_k": 300.0, "pressure_pa": 95_000.0,
                           "relative_humidity": fraction}
        }))
        .expect("decodes");
        serde_json::to_string(&density_altitude_v1(&request).expect("a sensor reading")).unwrap()
    };
    for (reading, end) in SENSOR_FRACTIONS {
        let at = |fraction: f64| {
            solve_v1_json(serde_json::json!({
                "temperature_k": 300.0, "pressure_pa": 95_000.0, "relative_humidity": fraction
            }))
            .unwrap_or_else(|e| panic!("relative_humidity {fraction} flies: {e}"))
        };
        // The whole response, the resolved request's echo included: the reading resolves to
        // the end, and flies it.
        assert_eq!(
            at(reading),
            at(end),
            "solve v1 at relative_humidity {reading}"
        );
        assert_eq!(
            density_altitude(reading),
            density_altitude(end),
            "density_altitude at relative_humidity {reading}"
        );
    }
    // Beyond the band: the same path and code as before S5b.
    for fraction in [1.06, 1.050_000_1, -0.06, -0.050_000_1] {
        let refused = solve_v1_json(serde_json::json!({"relative_humidity": fraction}))
            .expect_err("beyond the band");
        let refused: serde_json::Value = serde_json::from_str(&refused).unwrap();
        assert_eq!(refused["error"]["code"], "invalid_value", "{refused}");
        assert_eq!(
            refused["error"]["path"], "$.atmosphere.relative_humidity",
            "{refused}"
        );
    }
}

#[test]
fn card_hold_curve_and_corridor_read_a_sensor_reading_as_the_end_it_reads_as() {
    use ballistics_engine::card_service::{come_ups_v1, range_table_v1};
    use ballistics_engine::hold_curve::HoldCurve;
    use ballistics_engine::wind_scenarios::solve_robust_hold;
    let card = |humidity: f64| {
        let request = card_request(Some(humidity));
        (
            serde_json::to_string(&come_ups_v1(&request).expect("come-ups flies")).unwrap(),
            serde_json::to_string(&range_table_v1(&request).expect("range table flies")).unwrap(),
        )
    };
    let hold = |humidity: f64| {
        let curve = HoldCurve::solve(&hold_load(humidity), 300.0).expect("the hold curve flies");
        let mut bits = vec![curve.reach_m().to_bits()];
        for range in curve.sample_ranges_m() {
            if let Some(point) = curve.at_range(range) {
                bits.extend(
                    [
                        point.drop_mil,
                        point.wind_mil,
                        point.velocity_mps,
                        point.time_s,
                    ]
                    .map(f64::to_bits),
                );
            }
        }
        bits
    };
    let corridor = |humidity: f64| {
        serde_json::to_string(&solve_robust_hold(&corridor_request(humidity)).expect("flies"))
            .unwrap()
    };
    for (reading, end) in SENSOR_READINGS {
        assert_eq!(card(reading), card(end), "card at RH {reading}");
        assert_eq!(hold(reading), hold(end), "hold curve at RH {reading}");
        assert_eq!(corridor(reading), corridor(end), "corridor at RH {reading}");
    }
}

#[cfg(feature = "bridge")]
#[test]
fn the_bridge_reads_a_sensor_reading_as_the_end_it_reads_as() {
    let true_dsf = |humidity: f64| {
        bridge(
            "true.dsf",
            serde_json::json!({
                "model": {
                    "muzzle_velocity_fps": 1600.0, "ballistic_coefficient": 0.243,
                    "drag_model": "g7", "mass_gr": 168.0, "diameter_in": 0.308,
                    "zero_distance_yd": 100.0, "sight_height_in": 2.0,
                    "temperature_f": 59.0, "pressure_inhg": 29.92,
                    "humidity_pct": humidity, "altitude_ft": 0.0
                },
                "range_yd": 200.0, "observed_drop": 12.0, "drop_unit": "mil"
            }),
        )
    };
    // The apps' live shot path: solve-json v1 through the bridge's `solve`, and the
    // density-altitude command, at a wire fraction in the band.
    let wire = |command: &str, fraction: f64| {
        let atmosphere = serde_json::json!({
            "temperature_k": 300.0, "pressure_pa": 95_000.0, "relative_humidity": fraction
        });
        let request = if command == "solve" {
            serde_json::json!({
                "schema_version": 1,
                "projectile": {"mass_kg": 0.01134, "diameter_m": 0.00782, "length_m": 0.031,
                               "drag_model": "G7", "ballistic_coefficient": 0.243},
                "rifle": {"muzzle_velocity_mps": 823.0},
                "shot": {"max_range_m": 300.0, "zero_distance_m": 100.0},
                "atmosphere": atmosphere,
                "wind": {}, "solver": {}, "effects": {}, "sampling": {}
            })
        } else {
            serde_json::json!({"atmosphere": atmosphere})
        };
        bridge(command, request)
    };
    for (reading, end) in SENSOR_FRACTIONS {
        for command in ["solve", "atmosphere.density_altitude"] {
            let out = wire(command, reading);
            assert_eq!(out["ok"], true, "{command} at {reading} flies: {out}");
            assert_eq!(
                out,
                wire(command, end),
                "{command} at relative_humidity {reading}"
            );
        }
    }
    for fraction in [1.06, -0.06] {
        for command in ["solve", "atmosphere.density_altitude"] {
            let out = wire(command, fraction);
            assert_eq!(
                out["ok"], false,
                "{command} at {fraction} is refused: {out}"
            );
            assert!(
                out.to_string().contains("$.atmosphere.relative_humidity"),
                "{command} at {fraction}: located at the field: {out}"
            );
        }
    }
    for (reading, end) in SENSOR_READINGS {
        for command in ["card.come_ups", "card.range_table", "card.wind"] {
            let at = |humidity: f64| {
                let mut request = card_request_json();
                request["humidity"] = serde_json::json!(humidity);
                if command == "card.wind" {
                    request["wind_speeds"] = serde_json::json!([5.0, 10.0]);
                }
                bridge(command, request)
            };
            let out = at(reading);
            assert_eq!(out["ok"], true, "{command} RH {reading} flies: {out}");
            assert_eq!(out, at(end), "{command} at RH {reading}");
        }
        let out = true_dsf(reading);
        assert_eq!(out["ok"], true, "true.dsf RH {reading} flies: {out}");
        assert_eq!(out, true_dsf(end), "true.dsf at RH {reading}");
        let out = bridge(
            "profile.validate",
            serde_json::json!({
                "name": "humid", "velocity": 2700.0, "bc": 0.243, "mass": 175.0,
                "diameter": 0.308, "drag_model": "G7", "humidity": reading
            }),
        );
        assert_eq!(
            out["result"]["valid"], true,
            "a stored reading of {reading} is valid: {out}"
        );
    }
}

/// The online client sends a sensor reading as the end it reads as, and refuses a value beyond
/// the band before it sends anything.
#[cfg(feature = "online")]
#[test]
fn the_online_client_sends_a_sensor_reading_as_the_end_it_reads_as() {
    let request = |humidity: f64| {
        ballistics_engine::api_client::TrajectoryRequestBuilder::new()
            .bc_value(0.475)
            .bc_type("G1")
            .bullet_mass(10.9)
            .muzzle_velocity(823.0)
            .target_distance(500.0)
            .humidity(humidity)
            .build()
            .unwrap()
    };
    for (reading, end) in SENSOR_READINGS {
        let request_line = captured_query(&request(reading));
        let sent = request_line
            .split(['?', '&', ' '])
            .find_map(|pair| pair.strip_prefix("humidity="))
            .unwrap_or_else(|| panic!("no humidity in {request_line:?}"));
        assert_eq!(
            sent.parse::<f64>().unwrap().to_bits(),
            f64::to_bits(end),
            "RH {reading} is sent as {end}: sent {sent:?}"
        );
    }
    // The velocity-truing POST body too.
    for (reading, end) in SENSOR_READINGS {
        let body = captured_true_velocity_body(reading);
        let sent: serde_json::Value = serde_json::from_str(&body).expect("a JSON body");
        assert_eq!(
            sent["humidity"].as_f64().map(f64::to_bits),
            Some(f64::to_bits(end)),
            "RH {reading} is sent as {end}: {body}"
        );
    }
    // Nothing listens here: a request that got as far as the network would fail differently.
    let client = ballistics_engine::api_client::ApiClient::new("http://127.0.0.1:9", 1);
    for humidity in BAD_PERCENTS {
        match client.calculate_trajectory(&request(humidity)) {
            Err(ballistics_engine::api_client::ApiError::RequestError(message)) => assert!(
                message.contains("humidity"),
                "RH {humidity}: the refusal names the field: {message}"
            ),
            other => panic!("RH {humidity} must be refused before sending: {other:?}"),
        }
        match client.true_velocity(&true_velocity_request(humidity)) {
            Err(ballistics_engine::api_client::ApiError::RequestError(message)) => assert!(
                message.contains("humidity"),
                "RH {humidity}: the refusal names the field: {message}"
            ),
            other => panic!("RH {humidity} must be refused before sending: {other:?}"),
        }
    }
}

#[cfg(feature = "online")]
fn true_velocity_request(humidity: f64) -> ballistics_engine::api_client::TrueVelocityRequest {
    ballistics_engine::api_client::TrueVelocityRequest {
        measured_drop_mil: 3.6,
        range_yd: 500.0,
        bc: 0.475,
        drag_model: "G1".to_string(),
        weight_gr: 168.0,
        caliber: 0.308,
        zero_range_yd: Some(100.0),
        chrono_velocity_fps: None,
        altitude_ft: None,
        temperature_f: None,
        pressure_inhg: None,
        temperature_supplied: None,
        pressure_supplied: None,
        humidity: Some(humidity),
        sight_height_in: None,
        use_bc_enhancement: None,
    }
}

/// The JSON body `ApiClient::true_velocity` POSTs for a request at `humidity`, captured by a
/// one-shot loopback server that answers 500.
#[cfg(feature = "online")]
fn captured_true_velocity_body(humidity: f64) -> String {
    use std::io::{BufRead, BufReader, Read, Write};
    let listener = std::net::TcpListener::bind("127.0.0.1:0").expect("bind loopback");
    let port = listener.local_addr().unwrap().port();
    let server = std::thread::spawn(move || {
        let (stream, _) = listener.accept().expect("the client connects");
        let mut reader = BufReader::new(stream.try_clone().unwrap());
        let mut line = String::new();
        let mut length = 0usize;
        loop {
            line.clear();
            reader.read_line(&mut line).unwrap();
            if line.trim().is_empty() {
                break;
            }
            if let Some((name, value)) = line.split_once(':') {
                if name.eq_ignore_ascii_case("content-length") {
                    length = value.trim().parse().unwrap();
                }
            }
        }
        let mut body = vec![0u8; length];
        reader.read_exact(&mut body).unwrap();
        let mut stream = stream;
        stream
            .write_all(b"HTTP/1.1 500 Internal Server Error\r\nContent-Length: 0\r\nConnection: close\r\n\r\n")
            .unwrap();
        String::from_utf8(body).unwrap()
    });
    let client =
        ballistics_engine::api_client::ApiClient::new(&format!("http://127.0.0.1:{port}"), 5);
    let _ = client.true_velocity(&true_velocity_request(humidity));
    server.join().unwrap()
}
