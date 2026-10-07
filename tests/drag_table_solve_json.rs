//! `projectile.drag_table` on the solve-json v1 wire (MBA-1597).
//!
//! The engine has flown user drag decks since MBA-940, but only the CLI (`--drag-table`), the
//! WASM terminal and the array FFI could reach them: `ProjectileV1` is `deny_unknown_fields`,
//! so an app on the JSON bridge could not send a curve at all. Alfredo Mendiola Loyola asked for
//! it so an airgun app can fly GA2 — a pellet drag law distributed with MERO/GPC, whose BCs are
//! measured against that curve. This pins the wire contract:
//!
//!   (a) omitting the field changes nothing, and the key appears nowhere in the response;
//!   (b) `kind` decides what the BC means, with no default, because guessing wrong gives
//!       confidently wrong drops: `"reference"` is a standard curve the BC was measured
//!       against (flown exactly as G1 is flown with a G1 BC), `"projectile"` is the bullet's
//!       own Cd, which ignores the BC exactly as `--drag-table` does;
//!   (c) the bundled G1 and G7 tables sent as `"reference"` curves fly exactly like the
//!       built-in models at the same BC — the BC handling checked against an answer already
//!       trusted;
//!   (d) the resolved request echoes the table and re-solves identically from it, so the
//!       explain / error-budget / tolerance kernels, which perturb from that echo, keep it —
//!       and under the reference kind they see the BC, under the projectile kind they do not;
//!   (e) validation mirrors `--drag-table` and every error names the exact entry.
//!
//! The refusal to combine a table with the offline BC correction (`corrections`) is pinned
//! beside that correction's own bridge tests.

use ballistics_engine::drag::DragTable;
use ballistics_engine::solve_json::{
    decode_solve_request_v1, SolveErrorCodeV1, SolveErrorEnvelopeV1, SolveRequestV1, SolveSuccessV1,
};
use ballistics_engine::solve_v1;
use serde_json::{json, Value};

/// A .308 175 gr at 792 m/s, zeroed at 100 m, solved to 800 m in still ICAO air.
fn rifle_request(drag_model: &str, bc: f64) -> Value {
    json!({
        "schema_version": 1,
        "projectile": {"mass_kg": 0.01134, "diameter_m": 0.00782,
                       "drag_model": drag_model, "ballistic_coefficient": bc},
        "rifle": {"muzzle_velocity_mps": 792.48, "sight_height_m": 0.0508},
        "shot": {"max_range_m": 800.0, "zero_distance_m": 100.0},
        "atmosphere": {"altitude_m": 0.0, "temperature_k": 288.15, "pressure_pa": 101325.0,
                       "relative_humidity": 0.0},
        "wind": {"speed_mps": 0.0, "direction_from_rad": 0.0},
        "solver": {"method": "rk4", "time_step_s": 0.0005},
        "effects": {},
        "sampling": {"interval_m": 50.0}
    })
}

/// A .22 (5.5 mm) 18 gr pellet at 270 m/s, zeroed at 30 m, solved to 100 m — Alfredo's case.
fn pellet_request(drag_model: &str, bc: f64) -> Value {
    json!({
        "schema_version": 1,
        "projectile": {"mass_kg": 0.001166, "diameter_m": 0.0055,
                       "drag_model": drag_model, "ballistic_coefficient": bc},
        "rifle": {"muzzle_velocity_mps": 270.0, "sight_height_m": 0.05},
        "shot": {"max_range_m": 100.0, "zero_distance_m": 30.0},
        "atmosphere": {"altitude_m": 0.0, "temperature_k": 288.15, "pressure_pa": 101325.0,
                       "relative_humidity": 0.0},
        "wind": {"speed_mps": 0.0, "direction_from_rad": 0.0},
        "solver": {"method": "rk4", "time_step_s": 0.0005},
        "effects": {},
        "sampling": {"interval_m": 10.0}
    })
}

fn with_table(mut request: Value, kind: &str, points: Value) -> Value {
    request["projectile"]["drag_table"] = json!({"kind": kind, "points": points});
    request
}

/// The bundled reference table at `data/<name>.csv`, as wire points.
fn bundled_points(csv: &str) -> Value {
    let table = DragTable::from_csv_str(csv).expect("bundled table parses");
    Value::Array(
        table
            .mach_values
            .iter()
            .zip(&table.cd_values)
            .map(|(mach, cd)| json!({"mach": mach, "cd": cd}))
            .collect(),
    )
}

fn g1_points() -> Value {
    bundled_points(include_str!("../data/g1.csv"))
}

fn g7_points() -> Value {
    bundled_points(include_str!("../data/g7.csv"))
}

/// A small invented curve with a transonic rise no built-in model reproduces.
fn invented_points() -> Value {
    json!([
        {"mach": 0.0, "cd": 0.230}, {"mach": 0.5, "cd": 0.220}, {"mach": 0.8, "cd": 0.230},
        {"mach": 1.0, "cd": 0.520}, {"mach": 1.2, "cd": 0.480}, {"mach": 1.5, "cd": 0.400},
        {"mach": 2.0, "cd": 0.330}, {"mach": 2.5, "cd": 0.300}
    ])
}

fn decode(value: &Value) -> Result<SolveRequestV1, SolveErrorEnvelopeV1> {
    decode_solve_request_v1(&serde_json::to_string(value).expect("serialize request"))
}

fn solve(value: &Value) -> SolveSuccessV1 {
    solve_v1(decode(value).expect("request must decode")).expect("request must solve")
}

fn solve_error(value: &Value) -> SolveErrorEnvelopeV1 {
    match decode(value) {
        Err(error) => error,
        Ok(request) => solve_v1(request).expect_err("request must be refused"),
    }
}

fn terminal(response: &SolveSuccessV1) -> (f64, f64) {
    let last = response
        .samples
        .last()
        .expect("a solved trajectory has samples");
    (last.speed_mps, last.drop_m)
}

/// Sectional density in lb/in², the unit BCs are stated in.
fn sectional_density(request: &Value) -> f64 {
    let grains = request["projectile"]["mass_kg"].as_f64().unwrap() / 0.000_064_798_91;
    let inches = request["projectile"]["diameter_m"].as_f64().unwrap() / 0.0254;
    grains / 7000.0 / (inches * inches)
}

fn relative(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

// (a) ---------------------------------------------------------------------------------------

#[test]
fn omitting_the_table_leaves_no_trace_in_the_response() {
    let response = solve(&rifle_request("G7", 0.243));
    let text = serde_json::to_string(&response).unwrap();
    assert!(
        !text.contains("drag_table"),
        "an omitted drag_table must not appear in the response"
    );
}

// (b) ---------------------------------------------------------------------------------------

#[test]
fn kind_is_required() {
    let mut request = rifle_request("G7", 0.243);
    request["projectile"]["drag_table"] = json!({"points": invented_points()});
    let error = solve_error(&request);
    assert_eq!(error.error.code, SolveErrorCodeV1::MissingField);
    assert_eq!(error.error.path(), Some("$.projectile.drag_table.kind"));
}

#[test]
fn projectile_kind_ignores_the_bc() {
    let low = solve(&with_table(
        rifle_request("G7", 0.2),
        "projectile",
        invented_points(),
    ));
    let high = solve(&with_table(
        rifle_request("G7", 0.6),
        "projectile",
        invented_points(),
    ));
    assert_eq!(low.samples, high.samples);
}

#[test]
fn reference_kind_follows_the_bc() {
    let low = terminal(&solve(&with_table(
        rifle_request("G7", 0.2),
        "reference",
        g7_points(),
    )));
    let high = terminal(&solve(&with_table(
        rifle_request("G7", 0.3),
        "reference",
        g7_points(),
    )));
    assert!(
        high.0 > low.0,
        "a higher BC must keep more speed: {high:?} vs {low:?}"
    );
    // `drop_m` is positive below the line of sight.
    assert!(
        high.1 < low.1,
        "a higher BC must drop less: {high:?} vs {low:?}"
    );
}

#[test]
fn reference_kind_is_the_projectile_curve_scaled_by_sd_over_bc() {
    let bc = 0.243;
    let base = rifle_request("G7", bc);
    let scale = sectional_density(&base) / bc;
    let scaled = Value::Array(
        invented_points()
            .as_array()
            .unwrap()
            .iter()
            .map(|p| json!({"mach": p["mach"], "cd": p["cd"].as_f64().unwrap() * scale}))
            .collect(),
    );

    let reference = solve(&with_table(base.clone(), "reference", invented_points()));
    let projectile = solve(&with_table(base, "projectile", scaled));

    assert_eq!(reference.samples.len(), projectile.samples.len());
    for (r, p) in reference.samples.iter().zip(&projectile.samples) {
        assert!(relative(r.speed_mps, p.speed_mps) < 1e-9, "{r:?} vs {p:?}");
        assert!((r.drop_m - p.drop_m).abs() < 1e-9, "{r:?} vs {p:?}");
    }
}

// (c) ---------------------------------------------------------------------------------------

/// Every sample agrees to rounding. The solve-json path flies a built-in model as the bundled
/// reference table's Cd over the BC — the very tables these tests send — so a reference-kind
/// table is not an approximation of a built-in model but the same arithmetic reordered.
fn assert_flies_identically(table: &SolveSuccessV1, builtin: &SolveSuccessV1) {
    assert_eq!(table.samples.len(), builtin.samples.len());
    for (t, b) in table.samples.iter().zip(&builtin.samples) {
        assert!(relative(t.speed_mps, b.speed_mps) < 1e-12, "{t:?} vs {b:?}");
        assert!(
            (t.drop_m - b.drop_m).abs() < 1e-12 * b.drop_m.abs().max(1.0),
            "{t:?} vs {b:?}"
        );
    }
}

#[test]
fn g7_table_as_a_reference_flies_like_the_g7_model() {
    let builtin = solve(&rifle_request("G7", 0.243));
    let table = solve(&with_table(
        rifle_request("G7", 0.243),
        "reference",
        g7_points(),
    ));
    assert_flies_identically(&table, &builtin);
}

#[test]
fn g1_table_as_a_reference_flies_like_the_g1_model_for_a_pellet() {
    let builtin = solve(&pellet_request("G1", 0.028));
    let table = solve(&with_table(
        pellet_request("G1", 0.028),
        "reference",
        g1_points(),
    ));
    assert_flies_identically(&table, &builtin);
}

// (d) ---------------------------------------------------------------------------------------

#[test]
fn resolved_request_echoes_the_table_and_re_solves_identically() {
    let request = with_table(pellet_request("G1", 0.028), "reference", g1_points());
    let response = solve(&request);

    let echoed = serde_json::to_value(&response.resolved_request).unwrap();
    assert_eq!(
        echoed["projectile"]["drag_table"],
        request["projectile"]["drag_table"]
    );

    let rebuilt = SolveRequestV1::from(&response.resolved_request);
    let again = solve_v1(rebuilt).expect("resolved request must re-solve");
    assert_eq!(again.samples, response.samples);
}

/// The kernels behind `explain`, `error-budget` and `tolerance` perturb the resolved request and
/// re-solve it. A reference-kind BC must move the result there as a built-in BC does; a
/// projectile-kind BC must move nothing — which also proves the re-solves kept the table, since
/// on `drag_model`'s curve the BC would matter.
#[test]
fn the_perturbation_kernel_sees_the_table_and_its_kind() {
    use ballistics_engine::perturbation::{central_difference, InputAxis};

    let bc_slope = |kind: &str| {
        let response = solve(&with_table(rifle_request("G7", 0.243), kind, g7_points()));
        central_difference(
            &response.resolved_request,
            InputAxis::BallisticCoefficient,
            &[600.0],
            None,
        )
        .expect("the BC axis is continuous")[0]
            .d_drop_d_x
    };
    let builtin = central_difference(
        &solve(&rifle_request("G7", 0.243)).resolved_request,
        InputAxis::BallisticCoefficient,
        &[600.0],
        None,
    )
    .expect("the BC axis is continuous")[0]
        .d_drop_d_x;

    assert!(
        relative(bc_slope("reference"), builtin) < 1e-6,
        "reference {} vs built-in {builtin}",
        bc_slope("reference")
    );
    assert_eq!(bc_slope("projectile"), 0.0);
}

// (e) ---------------------------------------------------------------------------------------

fn expect_error(request: &Value, code: SolveErrorCodeV1, path: &str) {
    let error = solve_error(request);
    assert_eq!(error.error.code, code, "{}", error.error.message);
    assert_eq!(error.error.path(), Some(path), "{}", error.error.message);
}

#[test]
fn unknown_kind_is_rejected_at_the_kind() {
    let request = with_table(rifle_request("G7", 0.243), "measured", invented_points());
    expect_error(
        &request,
        SolveErrorCodeV1::InvalidValue,
        "$.projectile.drag_table.kind",
    );
}

#[test]
fn null_table_is_not_omission() {
    let mut request = rifle_request("G7", 0.243);
    request["projectile"]["drag_table"] = Value::Null;
    expect_error(
        &request,
        SolveErrorCodeV1::InvalidValue,
        "$.projectile.drag_table",
    );
}

#[test]
fn unknown_fields_are_rejected_with_their_path() {
    let mut request = with_table(rifle_request("G7", 0.243), "reference", g7_points());
    request["projectile"]["drag_table"]["name"] = json!("GA2");
    expect_error(
        &request,
        SolveErrorCodeV1::UnknownField,
        "$.projectile.drag_table.name",
    );

    let mut request = with_table(rifle_request("G7", 0.243), "reference", invented_points());
    request["projectile"]["drag_table"]["points"][1]["cdd"] = json!(0.2);
    expect_error(
        &request,
        SolveErrorCodeV1::UnknownField,
        "$.projectile.drag_table.points[1].cdd",
    );
}

#[test]
fn points_must_be_an_array_of_numbered_objects() {
    let request = with_table(
        rifle_request("G7", 0.243),
        "reference",
        json!({"mach": 0.5}),
    );
    expect_error(
        &request,
        SolveErrorCodeV1::InvalidValue,
        "$.projectile.drag_table.points",
    );

    let request = with_table(
        rifle_request("G7", 0.243),
        "reference",
        json!([{"mach": 0.0, "cd": 0.2}, {"mach": "0.5", "cd": 0.2}]),
    );
    expect_error(
        &request,
        SolveErrorCodeV1::InvalidValue,
        "$.projectile.drag_table.points[1].mach",
    );

    let request = with_table(
        rifle_request("G7", 0.243),
        "reference",
        json!([{"mach": 0.0, "cd": 0.2}, {"mach": 0.5}]),
    );
    expect_error(
        &request,
        SolveErrorCodeV1::MissingField,
        "$.projectile.drag_table.points[1].cd",
    );
}

fn ascending(count: usize) -> Value {
    Value::Array(
        (0..count)
            .map(|i| json!({"mach": i as f64 * 0.001, "cd": 0.2}))
            .collect(),
    )
}

#[test]
fn point_count_matches_the_cli_limits() {
    for count in [0, 1, 4097] {
        let request = with_table(rifle_request("G7", 0.243), "projectile", ascending(count));
        expect_error(
            &request,
            SolveErrorCodeV1::InvalidValue,
            "$.projectile.drag_table.points",
        );
    }
    for count in [2, 4096] {
        let request = with_table(rifle_request("G7", 0.243), "projectile", ascending(count));
        solve(&request);
    }
}

#[test]
fn mach_must_be_non_negative_and_strictly_ascending() {
    let request = with_table(
        rifle_request("G7", 0.243),
        "reference",
        json!([{"mach": -0.1, "cd": 0.2}, {"mach": 0.5, "cd": 0.2}]),
    );
    expect_error(
        &request,
        SolveErrorCodeV1::InvalidValue,
        "$.projectile.drag_table.points[0].mach",
    );

    let request = with_table(
        rifle_request("G7", 0.243),
        "reference",
        json!([{"mach": 0.0, "cd": 0.2}, {"mach": 0.9, "cd": 0.2}, {"mach": 0.9, "cd": 0.3}]),
    );
    expect_error(
        &request,
        SolveErrorCodeV1::InvalidValue,
        "$.projectile.drag_table.points[2].mach",
    );
}

#[test]
fn cd_must_be_positive() {
    for cd in [0.0, -0.2] {
        let request = with_table(
            rifle_request("G7", 0.243),
            "reference",
            json!([{"mach": 0.0, "cd": 0.2}, {"mach": 0.5, "cd": cd}]),
        );
        expect_error(
            &request,
            SolveErrorCodeV1::InvalidValue,
            "$.projectile.drag_table.points[1].cd",
        );
    }
}

/// Rust callers can build a request without the decoder; the service still refuses a bad table.
#[test]
fn the_service_refuses_a_bad_table_that_skipped_the_decoder() {
    let mut request = decode(&with_table(
        rifle_request("G7", 0.243),
        "reference",
        invented_points(),
    ))
    .unwrap();
    request.projectile.drag_table.as_mut().unwrap().points[3].cd = f64::NAN;
    let error = solve_v1(request).expect_err("a NaN Cd must be refused");
    assert_eq!(error.error.code, SolveErrorCodeV1::InvalidValue);
    assert_eq!(
        error.error.path(),
        Some("$.projectile.drag_table.points[3].cd")
    );
}

// Bridge ------------------------------------------------------------------------------------

#[cfg(feature = "bridge")]
#[test]
fn the_bridge_solve_command_accepts_a_table() {
    let envelope = json!({
        "api_version": 1,
        "command": "solve",
        "request": with_table(pellet_request("G1", 0.028), "reference", g1_points()),
    });
    let out: Value = serde_json::from_str(&ballistics_engine::bridge::bridge_call(
        &envelope.to_string(),
    ))
    .unwrap();
    assert_eq!(out["ok"], json!(true), "{out}");
    assert_eq!(
        out["result"]["resolved_request"]["projectile"]["drag_table"]["kind"],
        json!("reference")
    );
}
