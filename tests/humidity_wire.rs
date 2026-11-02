//! Humidity-consistency train, plan C5 (step S6): the SI JSON wire.
//!
//! The SI JSON wire family is solve v1 (the CLI's `solve-json`, the bridge's `solve`, the MCP
//! `solve` tool), the perturbation kernel, `explain`, the error budget, the resolved request's
//! reverse conversion and `atmosphere.density_altitude`. It keeps its ONE humidity field,
//! `relative_humidity`, a FRACTION beside `temperature_k`, `pressure_pa` and `latitude_rad`.
//! There is no percent alias, and every conversion to the percent the engine flies goes through
//! `RelativeHumidity::try_fraction`.
//!
//! Each wire command is pinned here three ways:
//!
//! - a round trip: the fraction sent is the fraction echoed (to the bit), the echo sent back
//!   flies the same bits, and the air flown is `StationAir`'s at `try_fraction(f)`;
//! - the sensor band (S5b, the owner's decision of 2026-11-01): a reading up to 0.05 outside
//!   0-1, such as a phone or a Kestrel above 100 % in fog, is the saturated or dry air it reads,
//!   and flies exactly as 1 or 0, so a high reading never fails a calculation;
//! - the refusal: a value beyond the band (a percent sent as a fraction, a typo) is refused at
//!   `$.atmosphere.relative_humidity`, never clamped, and impossible air is refused for the one
//!   shared reason.
//!
//! The bridge's per-command unit table (`include/ballistics_bridge.h`, `src/bridge/mod.rs`) is
//! checked against what each bridge command actually parses, command by command. The MCP
//! transport is pinned in `src/mcp_command.rs` (the module lives in the binary).

use ballistics_engine::atmosphere_service::{
    decode_density_altitude_request_v1, density_altitude_v1, DensityAltitudeResponseV1,
};
use ballistics_engine::error_budget::error_budget;
use ballistics_engine::explain::explain_difference;
use ballistics_engine::perturbation::{
    axis_meta, central_difference, evaluate, read_axis, with_axis, AxisKind, AxisValue, InputAxis,
    InputGroup, KernelError, Observation,
};
use ballistics_engine::solve_json::SolveErrorV1;
use ballistics_engine::{
    decode_solve_request_v1, solve_v1, AtmosphereErrorReason, RelativeHumidity,
    ResolvedSolveRequestV1, SolveErrorCodeV1, SolveRequestV1, SolveSuccessV1, StationAir,
};
use serde_json::{json, Value};

#[cfg(all(feature = "bridge", feature = "pdf", not(target_arch = "wasm32")))]
#[path = "support/pdf_text.rs"]
mod pdf_text;

/// The station every wire test flies: warm and a little high, so humidity matters.
const TEMPERATURE_K: f64 = 300.0;
const PRESSURE_PA: f64 = 95_000.0;

/// In-range fractions, ends included, and one with a long mantissa: each must come back with
/// its bits unchanged.
const FRACTIONS: [f64; 6] = [0.0, 0.3, 0.375, 0.5, 0.123_456_789_012_345_6, 1.0];

/// Sensor readings within 0.05 of 0-1, each with the end it reads as.
const SENSOR_FRACTIONS: [(f64, f64); 5] = [
    (1.000_000_1, 1.0),
    (1.03, 1.0),
    (1.05, 1.0),
    (-0.02, 0.0),
    (-0.05, 0.0),
];

/// Fractions beyond the band: a percent sent by mistake (50, 30), a typo, the next steps out.
const BAD_FRACTIONS: [f64; 6] = [50.0, 30.0, 1.5, -0.1, 1.050_000_1, -0.050_000_1];

const WIRE_PATH: &str = "$.atmosphere.relative_humidity";

fn solve_json(relative_humidity: Option<f64>) -> Value {
    let mut atmosphere = json!({"temperature_k": TEMPERATURE_K, "pressure_pa": PRESSURE_PA});
    if let Some(fraction) = relative_humidity {
        atmosphere["relative_humidity"] = json!(fraction);
    }
    json!({
        "schema_version": 1,
        "projectile": {"mass_kg": 0.01134, "diameter_m": 0.00782, "length_m": 0.031,
                       "drag_model": "G7", "ballistic_coefficient": 0.243},
        "rifle": {"muzzle_velocity_mps": 823.0},
        "shot": {"max_range_m": 600.0, "zero_distance_m": 100.0},
        "atmosphere": atmosphere,
        "wind": {"speed_mps": 4.0, "direction_from_rad": std::f64::consts::FRAC_PI_2},
        "solver": {}, "effects": {}, "sampling": {"interval_m": 100.0}
    })
}

fn request(relative_humidity: Option<f64>) -> SolveRequestV1 {
    decode_solve_request_v1(&solve_json(relative_humidity).to_string()).expect("decodes")
}

fn solved(relative_humidity: f64) -> SolveSuccessV1 {
    solve_v1(request(Some(relative_humidity)))
        .unwrap_or_else(|e| panic!("relative_humidity {relative_humidity} flies: {e:?}"))
}

fn refused(request: SolveRequestV1) -> SolveErrorV1 {
    solve_v1(request).expect_err("must be refused").error
}

fn resolved(relative_humidity: f64) -> ResolvedSolveRequestV1 {
    solved(relative_humidity).resolved_request
}

/// The flown part of a response, serialized: everything a re-submission must reproduce (the
/// notices are about what the REQUEST omitted, which a resolved echo never does).
fn flight(success: &SolveSuccessV1) -> String {
    serde_json::to_string(&(
        &success.resolved_request,
        &success.summary,
        &success.samples,
    ))
    .unwrap()
}

/// The station `StationAir` builds at the wire's fraction, by the one conversion.
fn station(relative_humidity: f64) -> StationAir {
    StationAir::new(
        0.0,
        TEMPERATURE_K - 273.15,
        PRESSURE_PA / 100.0,
        RelativeHumidity::try_fraction(relative_humidity).expect("a wire fraction"),
    )
    .expect("air that exists")
}

fn relative_error(actual: f64, expected: f64) -> f64 {
    ((actual - expected) / expected).abs()
}

// --- solve v1 -----------------------------------------------------------------------------------

#[test]
fn solve_v1_echoes_the_fraction_it_flies_and_its_echo_flies_the_same_bits() {
    for fraction in FRACTIONS {
        let success = solved(fraction);
        let echoed = success.resolved_request.atmosphere.relative_humidity;
        assert_eq!(echoed.to_bits(), fraction.to_bits(), "echo of {fraction}");

        // The air flown is StationAir's at try_fraction(f): every sample's Mach is its speed
        // over that station's sound speed. A fraction read in another unit (0.3 as 0.3 %)
        // moves the sound speed by about 1e-4, far outside this tolerance.
        let sound = station(fraction).speed_of_sound_mps();
        for sample in &success.samples {
            if sample.mach > 0.0 {
                assert!(
                    relative_error(sample.speed_mps / sample.mach, sound) < 1e-12,
                    "{fraction}: Mach at {} m is against {} m/s, the station's is {sound}",
                    sample.distance_m,
                    sample.speed_mps / sample.mach
                );
            }
        }

        // The resolved echo, sent back as a request through its JSON text (what an app that
        // stores and re-sends it does), flies the same bits.
        let rebuilt = SolveRequestV1::from(&success.resolved_request);
        assert_eq!(
            rebuilt.atmosphere.relative_humidity.map(f64::to_bits),
            Some(fraction.to_bits()),
            "the reverse conversion carries {fraction}"
        );
        let text = serde_json::to_string(&rebuilt).unwrap();
        let resent = decode_solve_request_v1(&text).expect("the rebuilt request decodes");
        assert_eq!(
            resent.atmosphere.relative_humidity.map(f64::to_bits),
            Some(fraction.to_bits()),
            "the JSON text carries {fraction}"
        );
        let again = solve_v1(resent).expect("the echo flies");
        assert_eq!(flight(&again), flight(&success), "re-solve at {fraction}");
    }
}

#[test]
fn solve_v1_reads_a_sensor_fraction_as_the_end_and_echoes_the_end() {
    for (reading, end) in SENSOR_FRACTIONS {
        let success = solved(reading);
        assert_eq!(
            success
                .resolved_request
                .atmosphere
                .relative_humidity
                .to_bits(),
            end.to_bits(),
            "{reading} resolves to {end}"
        );
        assert_eq!(
            serde_json::to_string(&success).unwrap(),
            serde_json::to_string(&solved(end)).unwrap(),
            "{reading} flies as {end}"
        );
    }
}

#[test]
fn solve_v1_refuses_a_fraction_beyond_the_band_at_its_path() {
    for fraction in BAD_FRACTIONS {
        let error = refused(request(Some(fraction)));
        assert_eq!(error.code, SolveErrorCodeV1::InvalidValue, "{fraction}");
        assert_eq!(error.path(), Some(WIRE_PATH), "{fraction}");
    }
}

#[test]
fn solve_v1_omitted_humidity_is_the_one_default_and_names_its_unit() {
    let omitted = solve_v1(request(None)).expect("flies");
    assert_eq!(
        flight(&omitted),
        flight(&solved(RelativeHumidity::DEFAULT.fraction())),
        "an omitted relative_humidity flies the default"
    );
    let notice = omitted
        .assumptions
        .iter()
        .find(|notice| notice.path.as_deref() == Some(WIRE_PATH))
        .expect("the default is announced at its path");
    assert_eq!(notice.message, "Relative humidity defaulted to 0.5 (50%).");
}

// --- the perturbation kernel --------------------------------------------------------------------

fn observations(request: &SolveRequestV1) -> Vec<u64> {
    evaluate(request, &[300.0, 600.0])
        .expect("evaluates")
        .iter()
        .flat_map(|o: &Observation| {
            [o.drop_m, o.windage_m, o.time_s, o.velocity_mps].map(f64::to_bits)
        })
        .collect()
}

#[test]
fn the_perturbation_axis_is_the_wire_fraction_in_and_out() {
    let meta = axis_meta(InputAxis::RelativeHumidity);
    assert_eq!(meta.group, InputGroup::Atmosphere);
    match meta.kind {
        AxisKind::Continuous { unit, .. } => {
            assert_eq!(unit, "", "the axis is the dimensionless wire fraction")
        }
        AxisKind::Categorical => panic!("relative humidity is continuous"),
    }

    for fraction in FRACTIONS {
        let base = resolved(fraction);
        assert_eq!(
            read_axis(&base, InputAxis::RelativeHumidity),
            Some(AxisValue::Scalar(fraction)),
            "read_axis gives the echoed fraction"
        );
    }

    // Moving the axis is sending that fraction: the rebuilt request flies the same bits as a
    // request sent at it, and a sensor reading moves it to the end.
    let base = resolved(0.3);
    for (moved, flies_as) in [
        (0.6, 0.6),
        (0.0, 0.0),
        (1.0, 1.0),
        (1.03, 1.0),
        (-0.02, 0.0),
    ] {
        let rebuilt = with_axis(&base, InputAxis::RelativeHumidity, AxisValue::Scalar(moved))
            .expect("rebuilds");
        assert_eq!(rebuilt.atmosphere.relative_humidity, Some(moved));
        assert_eq!(
            observations(&rebuilt),
            observations(&request(Some(flies_as))),
            "the axis moved to {moved}"
        );
    }

    // A percent where the fraction belongs is refused as the solve refuses it.
    for moved in [50.0, 1.5, -0.1] {
        let rebuilt = with_axis(&base, InputAxis::RelativeHumidity, AxisValue::Scalar(moved))
            .expect("rebuilds");
        match evaluate(&rebuilt, &[300.0]) {
            Err(KernelError::Solve { code, .. }) => {
                assert_eq!(code, SolveErrorCodeV1::InvalidValue, "{moved}")
            }
            other => panic!("the axis moved to {moved} must be refused: {other:?}"),
        }
    }
}

#[test]
fn the_humidity_derivative_is_per_unit_fraction() {
    let base = resolved(0.5);
    let derivative = central_difference(&base, InputAxis::RelativeHumidity, &[600.0], None)
        .expect("differentiates")
        .remove(0);
    let h = derivative.step_used;
    assert!(h > 0.0 && h < 0.01, "a step in fraction units: {h}");
    let drop_at =
        |fraction: f64| evaluate(&request(Some(fraction)), &[600.0]).expect("evaluates")[0].drop_m;
    let by_hand = (drop_at(0.5 + h) - drop_at(0.5 - h)) / (2.0 * h);
    assert!(
        relative_error(derivative.d_drop_d_x, by_hand) < 1e-9,
        "d drop / d fraction {} vs {by_hand} by hand",
        derivative.d_drop_d_x
    );

    // The error budget takes its sigma in the same unit: a 0.05 sigma (5 points) is that
    // derivative times 0.05.
    let report =
        error_budget(&base, &[(InputAxis::RelativeHumidity, 0.05)], &[600.0]).expect("budgets");
    assert!(report.unavailable_sources.is_empty(), "{report:?}");
    let row = &report.rows[0];
    let source = &row.sources[0];
    assert_eq!(source.axis, InputAxis::RelativeHumidity);
    assert_eq!(source.sigma, 0.05);
    assert_eq!(source.d_drop_d_x.to_bits(), derivative.d_drop_d_x.to_bits());
    assert!(
        relative_error(row.sigma_drop_m, (derivative.d_drop_d_x * 0.05).abs()) < 1e-12,
        "sigma_drop {} for a 0.05 sigma",
        row.sigma_drop_m
    );
}

#[test]
fn the_error_budget_refuses_a_base_humidity_beyond_the_band() {
    let mut base = resolved(0.5);
    base.atmosphere.relative_humidity = 50.0;
    let refused = error_budget(&base, &[(InputAxis::Temperature, 1.0)], &[300.0]);
    match refused {
        Err(KernelError::Solve { code, .. }) => assert_eq!(code, SolveErrorCodeV1::InvalidValue),
        Ok(report) => assert!(
            !report.unavailable_sources.is_empty(),
            "a base the solve refuses yields no budget: {report:?}"
        ),
        Err(other) => panic!("refused for the humidity, not {other:?}"),
    }
}

// --- explain ------------------------------------------------------------------------------------

#[test]
fn explain_puts_a_humidity_difference_in_the_atmosphere_group_in_wire_units() {
    let ranges = [300.0, 600.0];
    let a = resolved(0.3);
    let same = explain_difference(&a, &resolved(0.3), &ranges).expect("explains");
    for row in &same.rows {
        assert_eq!(row.total.drop_m, 0.0, "{row:?}");
        assert_eq!(row.total.windage_m, 0.0, "{row:?}");
    }

    let b = resolved(0.6);
    let report = explain_difference(&a, &b, &ranges).expect("explains");
    let at_a = evaluate(&request(Some(0.3)), &ranges).unwrap();
    let at_b = evaluate(&request(Some(0.6)), &ranges).unwrap();
    for (index, row) in report.rows.iter().enumerate() {
        let total = at_b[index].drop_m - at_a[index].drop_m;
        assert_eq!(row.total.drop_m.to_bits(), total.to_bits(), "{row:?}");
        assert!(total != 0.0, "humidity moves the drop at {} m", row.range_m);
        for contribution in &row.contributions {
            let expected = if contribution.group == InputGroup::Atmosphere {
                total
            } else {
                0.0
            };
            assert_eq!(
                contribution.delta.drop_m.to_bits(),
                expected.to_bits(),
                "{:?} at {} m",
                contribution.group,
                row.range_m
            );
        }
        assert_eq!(row.interaction_remainder.drop_m, 0.0, "{row:?}");
    }

    // A sensor reading resolves to the end before explain sees it: nothing to explain.
    let report = explain_difference(&resolved(1.0), &resolved(1.03), &ranges).expect("explains");
    for row in &report.rows {
        assert_eq!(row.total.drop_m, 0.0, "{row:?}");
    }

    // A resolved request carrying a percent where the fraction belongs is refused.
    let mut percent = resolved(0.3);
    percent.atmosphere.relative_humidity = 30.0;
    match explain_difference(&a, &percent, &ranges) {
        Err(KernelError::Solve { code, .. }) => assert_eq!(code, SolveErrorCodeV1::InvalidValue),
        other => panic!("a percent on the wire must be refused: {other:?}"),
    }
}

// --- atmosphere.density_altitude ----------------------------------------------------------------

fn density_altitude_at(atmosphere: Value) -> Result<DensityAltitudeResponseV1, SolveErrorV1> {
    let request = decode_density_altitude_request_v1(&json!({ "atmosphere": atmosphere }))
        .map_err(|envelope| envelope.error)?;
    density_altitude_v1(&request).map_err(|envelope| envelope.error)
}

fn density_altitude(relative_humidity: f64) -> Result<DensityAltitudeResponseV1, SolveErrorV1> {
    density_altitude_at(json!({
        "temperature_k": TEMPERATURE_K, "pressure_pa": PRESSURE_PA,
        "relative_humidity": relative_humidity
    }))
}

#[test]
fn density_altitude_echoes_the_fraction_and_describes_the_one_station() {
    for fraction in FRACTIONS {
        let response = density_altitude(fraction).expect("describes");
        assert_eq!(
            response.atmosphere.relative_humidity.to_bits(),
            fraction.to_bits(),
            "echo of {fraction}"
        );
        // The density reported is the station's the solve flies, at try_fraction(f).
        assert_eq!(
            response.air_density_kg_m3.to_bits(),
            station(fraction).density_kg_m3().to_bits(),
            "the density at {fraction} is StationAir's"
        );
        // The resolved atmosphere, sent back, describes the same air.
        let echo = serde_json::to_value(&response.atmosphere).unwrap();
        let again = density_altitude_at(echo).expect("the echo describes");
        assert_eq!(
            serde_json::to_string(&(
                &again.atmosphere,
                &again.pressure_altitude,
                &again.density_altitude,
                again.air_density_kg_m3
            ))
            .unwrap(),
            serde_json::to_string(&(
                &response.atmosphere,
                &response.pressure_altitude,
                &response.density_altitude,
                response.air_density_kg_m3
            ))
            .unwrap(),
            "the echo of {fraction}"
        );
    }
    for (reading, end) in SENSOR_FRACTIONS {
        assert_eq!(
            serde_json::to_string(&density_altitude(reading).expect("a reading")).unwrap(),
            serde_json::to_string(&density_altitude(end).unwrap()).unwrap(),
            "{reading} describes the air at {end}"
        );
    }
    for fraction in BAD_FRACTIONS {
        let error = density_altitude(fraction).expect_err("beyond the band");
        assert_eq!(error.code, SolveErrorCodeV1::InvalidValue, "{fraction}");
        assert_eq!(error.path(), Some(WIRE_PATH), "{fraction}");
    }
}

/// One impossible-air rule on the wire: the solve and `atmosphere.density_altitude` refuse the
/// same air for the same reason, `StationAir::new`'s, stated word for word.
#[test]
fn the_wire_refuses_impossible_air_for_the_one_shared_reason() {
    let reason = AtmosphereErrorReason::ImpossibleAir.to_string();
    // Saturated air at the boiling point, and at 100 C and 1000 hPa (the saturation pressure
    // alone is above the total).
    for (temperature_k, pressure_pa) in [(373.15, 101_325.0), (373.15, 100_000.0)] {
        let atmosphere = json!({
            "temperature_k": temperature_k, "pressure_pa": pressure_pa, "relative_humidity": 1.0
        });
        let station = StationAir::new(
            0.0,
            temperature_k - 273.15,
            pressure_pa / 100.0,
            RelativeHumidity::try_fraction(1.0).unwrap(),
        )
        .expect_err("StationAir refuses it");
        assert_eq!(station.reason(), AtmosphereErrorReason::ImpossibleAir);

        let error = density_altitude_at(atmosphere.clone()).expect_err("refused");
        assert_eq!(error.code, SolveErrorCodeV1::InvalidValue);
        assert_eq!(error.path(), Some(WIRE_PATH));
        assert!(
            error.message.contains(&reason),
            "density_altitude states the shared reason: {}",
            error.message
        );

        let density_altitude_message = error.message.clone();
        let mut solve = solve_json(Some(1.0));
        solve["atmosphere"] = atmosphere;
        let request = decode_solve_request_v1(&solve.to_string()).unwrap();
        let error = refused(request);
        // R2: one envelope, not only one reason. The solve refuses this air before any flight
        // as `invalid_value` at the humidity's path, with density_altitude's message word for
        // word (it was `solve_failed` with no path, in percent, from inside the zero search).
        assert_eq!(error.code, SolveErrorCodeV1::InvalidValue);
        assert_eq!(error.path(), Some(WIRE_PATH));
        assert_eq!(
            error.message, density_altitude_message,
            "the solve and density_altitude refuse it in the same words"
        );
        assert!(
            error.message.contains(&reason),
            "the solve states the shared reason: {}",
            error.message
        );

        // The same air at half the humidity exists, and both describe it.
        let half = json!({
            "temperature_k": temperature_k, "pressure_pa": pressure_pa, "relative_humidity": 0.5
        });
        density_altitude_at(half.clone()).expect("half the humidity exists");
        let mut solve = solve_json(Some(0.5));
        solve["atmosphere"] = half;
        solve_v1(decode_solve_request_v1(&solve.to_string()).unwrap()).expect("and flies");
    }
}

// --- the bridge's per-command unit table --------------------------------------------------------

// Host builds only: the bridge lists a command compiled out on wasm32 alone (it reads a file)
// in every host build, so on the host the table and `meta.capabilities` differ only by the
// feature-gated commands in `GATED`.
#[cfg(all(feature = "bridge", not(target_arch = "wasm32")))]
mod bridge_units {
    use super::*;
    use ballistics_engine::bridge::bridge_call;

    /// One row of the per-command humidity table: the command, the humidity field it reads
    /// (`-` for none), the unit (`fraction`, `percent` or `none`) and what an omitted field
    /// means (a number, `required`, or `-`).
    #[derive(Debug, Clone, PartialEq)]
    struct Row {
        command: String,
        field: String,
        unit: String,
        if_omitted: String,
    }

    /// Commands present only in some host builds: behind the `pdf`, `profile-import` and
    /// `profile-export` features. A command absent only on wasm32 is always listed on the
    /// host, where this module runs, so it needs no entry.
    const GATED: [&str; 3] = ["card.pdf", "profile.import_a7p", "profile.export_a7p"];

    fn source(path: &str) -> String {
        std::fs::read_to_string(std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(path))
            .unwrap_or_else(|e| panic!("{path}: {e}"))
    }

    /// The table in `include/ballistics_bridge.h`: the comment lines after the one that heads
    /// it (`command  humidity field  unit  if omitted`), up to the first blank comment line.
    fn header_table() -> Vec<Row> {
        let text = source("include/ballistics_bridge.h");
        let mut rows = Vec::new();
        let mut in_table = false;
        for line in text.lines() {
            let body = line.trim_start().trim_start_matches('*').trim();
            if !in_table {
                in_table = body.split_whitespace().collect::<Vec<_>>()
                    == ["command", "humidity", "field", "unit", "if", "omitted"];
                continue;
            }
            if body.is_empty() {
                break;
            }
            let cells: Vec<&str> = body.split_whitespace().collect();
            assert_eq!(
                cells.len(),
                4,
                "a header table row has four cells: {line:?}"
            );
            rows.push(Row {
                command: cells[0].to_string(),
                field: cells[1].to_string(),
                unit: cells[2].to_string(),
                if_omitted: cells[3].to_string(),
            });
        }
        assert!(
            in_table,
            "include/ballistics_bridge.h has no humidity unit table"
        );
        rows
    }

    /// The same table in `src/bridge/mod.rs`'s module docs, as a markdown table whose cells are
    /// code spans or bare words.
    fn module_doc_table() -> Vec<Row> {
        let text = source("src/bridge/mod.rs");
        let mut rows = Vec::new();
        let mut in_table = false;
        for line in text.lines() {
            let Some(doc) = line.strip_prefix("//!") else {
                if in_table {
                    break;
                }
                continue;
            };
            let doc = doc.trim();
            let cells: Vec<String> = doc
                .trim_matches('|')
                .split('|')
                .map(|cell| cell.trim().trim_matches('`').to_string())
                .collect();
            if !in_table {
                in_table = doc.starts_with('|')
                    && cells == ["command", "humidity field", "unit", "if omitted"];
                continue;
            }
            if !doc.starts_with('|') {
                break;
            }
            if cells.iter().all(|cell| cell.chars().all(|c| c == '-')) {
                continue;
            }
            assert_eq!(
                cells.len(),
                4,
                "a module-doc table row has four cells: {line:?}"
            );
            rows.push(Row {
                command: cells[0].clone(),
                field: cells[1].clone(),
                unit: cells[2].clone(),
                if_omitted: cells[3].clone(),
            });
        }
        assert!(in_table, "src/bridge/mod.rs has no humidity unit table");
        rows
    }

    fn call(command: &str, request: Value) -> Value {
        let envelope = json!({"api_version": 1, "command": command, "request": request});
        serde_json::from_str(&bridge_call(&envelope.to_string())).unwrap()
    }

    fn capabilities() -> Vec<String> {
        let out = call("meta.capabilities", Value::Null);
        serde_json::from_value(out["result"]["commands"].clone()).unwrap()
    }

    /// What a command made of a humidity: its result, or the refusal's text. Two results are
    /// compared through [`read_back`].
    #[derive(Debug, Clone, PartialEq)]
    enum Outcome {
        Flew(Value),
        Refused(String),
    }

    impl Outcome {
        fn flew(&self) -> &Value {
            match self {
                Self::Flew(value) => value,
                Self::Refused(message) => panic!("refused: {message}"),
            }
        }

        fn refusal(&self) -> &str {
            match self {
                Self::Refused(message) => message,
                Self::Flew(value) => panic!("flew: {value}"),
            }
        }
    }

    fn outcome(out: Value) -> Outcome {
        if out["ok"] == true {
            Outcome::Flew(out["result"].clone())
        } else {
            Outcome::Refused(out["error"].to_string())
        }
    }

    /// `value` at the dotted `path` of `request`, or the field removed when `None`.
    fn set(request: &mut Value, path: &str, value: Option<f64>) {
        let (parents, leaf) = match path.rsplit_once('.') {
            Some((parents, leaf)) => (Some(parents), leaf),
            None => (None, path),
        };
        let mut target = &mut *request;
        if let Some(parents) = parents {
            for key in parents.split('.') {
                target = &mut target[key];
            }
        }
        let object = target
            .as_object_mut()
            .expect("the field's parent is an object");
        match value {
            Some(value) => {
                object.insert(leaf.to_string(), json!(value));
            }
            None => {
                object.remove(leaf);
            }
        }
    }

    fn card_request() -> Value {
        json!({
            "units": "imperial", "muzzle_velocity": 2600.0, "ballistic_coefficient": 0.243,
            "mass": 175.0, "diameter": 0.308, "drag_model": "g7", "zero_distance": 100.0,
            "start": 100.0, "end": 300.0, "step": 100.0,
            "temperature": 95.0, "pressure": 28.5
        })
    }

    fn truing_model() -> Value {
        json!({
            "muzzle_velocity_fps": 2700.0, "ballistic_coefficient": 0.243,
            "drag_model": "g7", "mass_gr": 168.0, "diameter_in": 0.308,
            "zero_distance_yd": 100.0, "sight_height_in": 2.0,
            "temperature_f": 95.0, "pressure_inhg": 28.5,
            "humidity_pct": 50.0, "altitude_ft": 0.0
        })
    }

    fn profile() -> Value {
        json!({
            "name": "humid", "velocity": 792.0, "bc": 0.381, "mass": 19.439673,
            "diameter": 8.5852, "drag_model": "G7", "twist_rate": 254.0, "sight_height": 90.0,
            "zero_distance": 100.0, "units": "metric", "temperature": 35.0,
            "pressure": 950.0, "humidity": 50.0, "altitude": 0.0, "bullet_length": 45.72
        })
    }

    /// The request a humidity-carrying command is probed with, the humidity field at its
    /// table path.
    fn base_request(command: &str) -> Value {
        match command {
            "solve" => solve_json(Some(0.5)),
            "atmosphere.density_altitude" => json!({"atmosphere": {
                "temperature_k": TEMPERATURE_K, "pressure_pa": PRESSURE_PA,
                "relative_humidity": 0.5
            }}),
            "card.come_ups" | "card.range_table" => card_request(),
            "card.wind" => {
                let mut request = card_request();
                request["wind_speeds"] = json!([5.0, 10.0]);
                request
            }
            "card.pdf" => {
                let mut request = card_request();
                request["pdf"] = json!({"title": "Humid"});
                request
            }
            "profile.validate" | "profile.normalize" | "profile.export_a7p" => profile(),
            "true.fit" => json!({
                "model": truing_model(), "drop_unit": "mil",
                "observations": [
                    {"range_yd": 400.0, "drop": 2.6882, "sigma": 0.1},
                    {"range_yd": 600.0, "drop": 4.8626, "sigma": 0.1},
                    {"range_yd": 800.0, "drop": 7.6598, "sigma": 0.1}
                ],
                "priors": {"muzzle_velocity_fps": {"mean": 2700.0, "sigma": 20.0},
                           "ballistic_coefficient": {"mean": 0.243, "sigma": 0.01}},
                "predictions": []
            }),
            "true.dsf" => {
                let mut model = truing_model();
                model["muzzle_velocity_fps"] = json!(1600.0);
                json!({"model": model, "range_yd": 200.0, "observed_drop": 12.0,
                       "drop_unit": "mil"})
            }
            "true.plan" => json!({
                "model": truing_model(),
                "candidate_ranges_yd": [200.0, 400.0, 600.0, 800.0],
                "observation_count": 3, "minimum_separation_yd": 100.0,
                "measurement_sigma_1sd": 0.1, "drop_unit": "mil"
            }),
            "true.wind" => json!({
                "observations": [{"range_m": 457.2, "miss_right_m": 0.315, "sigma_m": null}],
                "muzzle_velocity_fps": 2700.0, "bc": 0.243, "drag_model": "g7",
                "mass_gr": 168.0, "diameter_in": 0.308, "zero_distance_yd": 100.0,
                "sight_height_in": 2.0, "temperature_f": 95.0, "pressure_inhg": 28.5,
                "humidity_pct": 50.0, "altitude_ft": 0.0,
                "twist": {"rate_in": 11.0, "right_hand": true},
                "earth": null, "called_crosswind_mph": null
            }),
            other => panic!("no probe request for {other}"),
        }
    }

    /// The part of a command's result that a humidity decides, so that two humidities the
    /// command reads as the same air give equal values. Results that echo their input (a
    /// resolved request, a profile) are compared whole: the echo is part of what was read.
    fn read_back(command: &str, result: &Value) -> Value {
        match command {
            // A verdict, not a flight: whether the stored humidity is one, and why not.
            "profile.validate" => json!([result["valid"], result["warnings"]]),
            // The notices say what the REQUEST omitted (the default's own is pinned by the
            // solve v1 tests above), not what flew.
            "solve" | "atmosphere.density_altitude" => {
                let mut flown = result.clone();
                if let Some(object) = flown.as_object_mut() {
                    object.remove("assumptions");
                }
                flown
            }
            // The profile imported. The report beside it describes the FILE, and says what the
            // file held (`103 %`) as well as what was imported from it (`100 %`).
            "profile.import_a7p" => result["profile"].clone(),
            // The printed card, read back as text, and everything reported beside it. The
            // document's bytes are not compared: its footer prints the time it was made.
            "card.pdf" => {
                let mut reported = result.clone();
                if let Some(object) = reported.as_object_mut() {
                    object.remove("pdf_base64");
                    object.remove("byte_length");
                }
                json!({"printed": printed_card(result), "reported": reported})
            }
            _ => result.clone(),
        }
    }

    /// `card.pdf`'s drawn text, as tokens, without the footer's clock time.
    #[cfg(feature = "pdf")]
    fn printed_card(result: &Value) -> Vec<String> {
        let pdf = base64_decode(result["pdf_base64"].as_str().expect("pdf_base64"));
        let is_clock = |token: &str| {
            token.len() == 8
                && token.char_indices().all(|(i, c)| {
                    if i == 2 || i == 5 {
                        c == ':'
                    } else {
                        c.is_ascii_digit()
                    }
                })
        };
        let tokens = pdf_text::tokens(&pdf, "card.pdf");
        assert!(
            tokens.iter().any(|token| is_clock(token)),
            "the footer's clock time is where it was: {tokens:?}"
        );
        tokens
            .into_iter()
            .filter(|token| !is_clock(token))
            .collect()
    }

    #[cfg(not(feature = "pdf"))]
    fn printed_card(_result: &Value) -> Vec<String> {
        unreachable!("card.pdf is not in this build")
    }

    /// The humidity a command's result carries, for the commands whose result is a profile:
    /// an echo for `profile.validate` and `profile.normalize`, the value read from the file for
    /// `profile.import_a7p`.
    fn carried_humidity(command: &str, result: &Value) -> Option<f64> {
        match command {
            "profile.validate" => result["normalized"]["humidity"].as_f64(),
            "profile.normalize" | "profile.import_a7p" => result["profile"]["humidity"].as_f64(),
            _ => None,
        }
    }

    /// `command` with its humidity field at `value` (`None`: omitted), through the bridge.
    /// `profile.validate`'s verdict against a stored humidity is its refusal.
    fn probe(row: &Row, value: Option<f64>) -> Outcome {
        let out = if row.command == "profile.import_a7p" {
            call(
                "profile.import_a7p",
                json!({"a7p_base64": base64(&a7p(value.map(|v| v as i64)))}),
            )
        } else {
            let mut request = base_request(&row.command);
            set(&mut request, &row.field, value);
            call(&row.command, request)
        };
        match outcome(out) {
            Outcome::Flew(result)
                if row.command == "profile.validate" && result["valid"] != true =>
            {
                Outcome::Refused(result["warnings"].to_string())
            }
            other => other,
        }
    }

    /// Check one row against what its command does: the unit (a percent sent as a fraction is
    /// refused, and a percent's command reads 30 as 30 % and 0.3 as 0.3 %, never guessing the
    /// unit from the size), the sensor band in that unit, a value beyond the band refused, and
    /// what an omitted field means.
    fn check(row: &Row) {
        let label = row.command.as_str();
        let (typical, saturated, reading, beyond) = match row.unit.as_str() {
            "fraction" => (0.3, 1.0, 1.03, 30.0),
            "percent" => (30.0, 100.0, 103.0, 150.0),
            other => panic!("{label}: unit {other:?}"),
        };
        let at = |value: f64| probe(row, Some(value));
        let same =
            |a: &Outcome, b: &Outcome| read_back(label, a.flew()) == read_back(label, b.flew());
        let from_file = row.field.starts_with("a7p:");

        // profile.normalize re-serializes a stored document: it reads nothing, and echoes the
        // percent it was given unchanged (profile.validate is the command that judges it).
        if label == "profile.normalize" {
            for value in [typical, 0.3, reading, beyond] {
                let echoed = carried_humidity(label, at(value).flew());
                assert_eq!(echoed, Some(value), "{label} echoes {value}");
            }
            let omitted = carried_humidity(label, probe(row, None).flew());
            assert_eq!(
                omitted.map(|h| h.to_string()),
                Some(row.if_omitted.clone()),
                "{label}: an omitted humidity is {}",
                row.if_omitted
            );
            return;
        }

        // The unit.
        let typical_outcome = at(typical);
        if let Some(carried) = carried_humidity(label, typical_outcome.flew()) {
            assert_eq!(carried, typical, "{label} carries {typical} as {typical}");
        }
        if row.unit == "fraction" {
            // A percent sent where the fraction belongs is refused, not read as 30 %.
            let refusal = at(30.0);
            assert!(
                refusal.refusal().contains("relative_humidity"),
                "{label}: 30 is refused at its field: {refusal:?}"
            );
        } else {
            // 0.3 is 0.3 %, a humidity, and another one than 30 %. (A verdict cannot show which
            // humidity it read, only that 0.3 is one.)
            let low = at(0.3);
            if label != "profile.validate" {
                assert!(
                    !same(&low, &typical_outcome),
                    "{label}: 0.3 is not read as 30"
                );
            }
            low.flew();
        }

        // The sensor band: a reading just past saturation is saturated air. A file field holds
        // a whole percent, so its reading is whole.
        let reading = if from_file { reading.round() } else { reading };
        let read = at(reading);
        assert!(
            same(&read, &at(saturated)),
            "{label}: {reading} reads as {saturated}"
        );
        if from_file {
            assert_eq!(
                carried_humidity(label, read.flew()),
                Some(saturated),
                "{label} imports {reading} as {saturated}"
            );
        }

        // Beyond the band: refused, and the refusal names a humidity.
        let refusal = at(beyond);
        assert!(
            refusal.refusal().to_lowercase().contains("humidity"),
            "{label}: {beyond} is refused for its humidity: {refusal:?}"
        );

        // Omitted.
        let omitted = probe(row, None);
        if row.if_omitted == "required" {
            let field = row.field.rsplit('.').next().unwrap();
            assert!(
                omitted.refusal().contains(field),
                "{label}: an omitted humidity is refused for the missing {field}: {omitted:?}"
            );
            return;
        }
        let default: f64 = row.if_omitted.parse().unwrap_or_else(|_| {
            panic!(
                "{label}: if omitted {:?} is a number or `required`",
                row.if_omitted
            )
        });
        assert!(
            same(&omitted, &at(default)),
            "{label}: omitted is {default}"
        );
        if let Some(carried) = carried_humidity(label, omitted.flew()) {
            assert_eq!(
                carried, default,
                "{label}: an omitted humidity is {default}"
            );
        }
        if row.unit == "fraction" {
            assert_eq!(default, RelativeHumidity::DEFAULT.fraction(), "{label}");
        } else if from_file {
            // The format's own meaning: proto3 cannot tell an absent int32 from 0.
            assert_eq!(default, 0.0, "{label}");
        } else {
            assert_eq!(default, RelativeHumidity::DEFAULT.percent(), "{label}");
        }
    }

    #[test]
    fn the_header_and_the_module_docs_carry_the_same_unit_table() {
        assert_eq!(header_table(), module_doc_table());
    }

    #[test]
    fn the_unit_table_names_every_command_once_and_follows_the_wire_rule() {
        let rows = header_table();
        let mut named: Vec<&str> = rows.iter().map(|row| row.command.as_str()).collect();
        named.sort_unstable();
        let before = named.len();
        named.dedup();
        assert_eq!(named.len(), before, "a command is listed twice: {rows:?}");
        for command in capabilities() {
            assert!(
                named.contains(&command.as_str()),
                "{command} is a bridge command, so the unit table must say whether it takes a \
                 humidity and in what unit"
            );
        }
        let available = capabilities();
        for row in &rows {
            assert!(
                available.contains(&row.command) || GATED.contains(&row.command.as_str()),
                "{} is in the table but is no bridge command",
                row.command
            );
            // The wire rule: a field named relative_humidity is a fraction, every other
            // humidity is a percent.
            match row.unit.as_str() {
                "fraction" => assert!(row.field.ends_with("relative_humidity"), "{row:?}"),
                "percent" => assert!(
                    !row.field.ends_with("relative_humidity") && row.field.contains("humidity"),
                    "{row:?}"
                ),
                "none" => assert_eq!((row.field.as_str(), row.if_omitted.as_str()), ("-", "-")),
                other => panic!("unit {other:?}: {row:?}"),
            }
        }
    }

    #[test]
    fn every_command_parses_humidity_as_its_table_row_says() {
        let available = capabilities();
        let mut checked = 0;
        for row in header_table() {
            if row.unit == "none" || !available.contains(&row.command) {
                continue;
            }
            check(&row);
            checked += 1;
        }
        assert!(checked >= 11, "only {checked} commands were probed");
    }

    // --- an .a7p file, for profile.import_a7p's probe --------------------------------------------

    #[cfg(feature = "profile-import")]
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

    #[cfg(feature = "profile-import")]
    fn enc_i32(number: u32, value: i64, out: &mut Vec<u8>) {
        enc_varint(u64::from(number) << 3, out);
        enc_varint(value as u64, out);
    }

    #[cfg(feature = "profile-import")]
    fn enc_bytes(number: u32, payload: &[u8], out: &mut Vec<u8>) {
        enc_varint((u64::from(number) << 3) | 2, out);
        enc_varint(payload.len() as u64, out);
        out.extend_from_slice(payload);
    }

    /// A G7 .a7p with field 17 (`c_zero_air_humidity`, an int32 percent) set to `humidity`,
    /// or absent when `None`.
    #[cfg(feature = "profile-import")]
    fn a7p(humidity: Option<i64>) -> Vec<u8> {
        let mut p = Vec::new();
        enc_bytes(1, b"HUMID 6.5CM", &mut p);
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
        ballistics_engine::profile_import::wrap_payload(&payload)
    }

    #[cfg(not(feature = "profile-import"))]
    fn a7p(_humidity: Option<i64>) -> Vec<u8> {
        unreachable!("profile.import_a7p is not in this build")
    }

    #[cfg(feature = "pdf")]
    fn base64_decode(text: &str) -> Vec<u8> {
        const ALPHABET: &[u8; 64] =
            b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        let (mut acc, mut bits, mut out) = (0u32, 0u32, Vec::new());
        for c in text.bytes().take_while(|&c| c != b'=') {
            let sextet = ALPHABET.iter().position(|&a| a == c).expect("base64") as u32;
            acc = (acc << 6) | sextet;
            bits += 6;
            if bits >= 8 {
                bits -= 8;
                out.push((acc >> bits) as u8);
            }
        }
        out
    }

    fn base64(bytes: &[u8]) -> String {
        const ALPHABET: &[u8; 64] =
            b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        let mut out = String::new();
        for chunk in bytes.chunks(3) {
            let triple = (u32::from(chunk[0]) << 16)
                | (u32::from(chunk.get(1).copied().unwrap_or(0)) << 8)
                | u32::from(chunk.get(2).copied().unwrap_or(0));
            out.push(char::from(ALPHABET[(triple >> 18) as usize & 63]));
            out.push(char::from(ALPHABET[(triple >> 12) as usize & 63]));
            out.push(if chunk.len() > 1 {
                char::from(ALPHABET[(triple >> 6) as usize & 63])
            } else {
                '='
            });
            out.push(if chunk.len() > 2 {
                char::from(ALPHABET[triple as usize & 63])
            } else {
                '='
            });
        }
        out
    }
}

// --- .a7p export of a sensor reading --------------------------------------------------------------

/// A stored humidity within 5 points of 0-100 % (a profile saved from a sensor in fog) is the
/// saturated or dry air it reads, so it exports as that end, exactly as the engine flies it and
/// as `.a7p` import reads such a file. Before S6 the export refused 102 % outright.
#[cfg(feature = "profile-export")]
#[test]
fn a7p_export_writes_a_sensor_reading_as_the_end_it_reads_as() {
    use ballistics_engine::profile::ProfileData;
    use ballistics_engine::profile_export::export_a7p;
    let profile = |humidity: f64| -> ProfileData {
        serde_json::from_value(json!({
            "name": "humid", "velocity": 792.0, "bc": 0.381, "mass": 19.439673,
            "diameter": 8.5852, "drag_model": "G7", "twist_rate": 254.0, "sight_height": 90.0,
            "zero_distance": 100.0, "units": "metric", "temperature": 15.0,
            "pressure": 1000.0, "humidity": humidity, "altitude": 0.0, "bullet_length": 45.72
        }))
        .unwrap()
    };
    let bytes = |humidity: f64| {
        export_a7p(&profile(humidity))
            .unwrap_or_else(|e| panic!("{humidity} % exports: {e}"))
            .bytes
    };
    for (reading, end) in [(100.000_000_1, 100.0), (102.0, 100.0), (105.0, 100.0)] {
        assert_eq!(bytes(reading), bytes(end), "{reading} % exports as {end} %");
    }
    for (reading, end) in [(-0.0001, 0.0), (-2.0, 0.0), (-5.0, 0.0)] {
        assert_eq!(bytes(reading), bytes(end), "{reading} % exports as {end} %");
    }
    for beyond in [150.0, 105.000_000_1, -5.000_000_1] {
        assert!(
            export_a7p(&profile(beyond)).is_err(),
            "{beyond} % is no reading of real air and is refused"
        );
    }
}
