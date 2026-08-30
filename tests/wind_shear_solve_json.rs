//! `effects.wind_shear_model` on the solve-json v1 wire (0.36.0).
//!
//! The engine has modelled altitude-dependent wind since MBA-728, and the CLI
//! (`--enable-wind-shear --wind-shear-model ...`) and the C ABI both reach it. The JSON solve
//! path — the mobile embedding surface — hard-coded `enable_wind_shear: false` instead, so an
//! app built on it had no shear at all. This pins the wire contract that closes that gap:
//!
//!   (a) omitting the field is byte-identical to every pre-0.36.0 response, echo included:
//!       the key does not appear anywhere in the envelope;
//!   (b) an explicit `"none"` solves identically to omitting it — the one state that must
//!       NOT quietly become power-law shear, because `cli_api`'s shear branch maps every
//!       unrecognized model name (including `"none"`) to the power law, so deriving
//!       `enable_wind_shear` from the model rather than carrying a separate flag is what
//!       keeps that state unreachable;
//!   (c) every accepted model is applied and echoed back in its canonical spelling, with
//!       the `"ekman"` alias canonicalized to `"ekman_spiral"`;
//!   (d) the shear actually MOVES the numbers: on a lofted shot the terminal windage shifts
//!       by tens of percent, and the logarithmic and power-law profiles differ from each
//!       other, so this cannot pass on a field that merely parses;
//!   (e) an unknown model is a structured `invalid_value` at the exact path, never a silent
//!       fall back to no shear — the failure mode this whole feature exists to avoid;
//!   (f) the model survives the resolved-request round trip;
//!   (g) downrange segments plus shear is a `conflicting_fields` error, matching the CLI and
//!       WASM front ends: the solver's wind sock takes precedence over its shear branch, so
//!       accepting the pair would drop the shear while still echoing a model as applied.

use ballistics_engine::solve_json::{
    decode_solve_request_v1, SolveRequestV1, SolveSuccessV1, WindShearModelV1,
};
use ballistics_engine::solve_v1;

/// A flat 900 m shot: the boundary-layer ratio is floored at 1.0 below the 10 m
/// meteorological reference height, so this trajectory is deliberately one where shear is
/// inert. Used for the contract-shape assertions, never for the "numbers move" one.
fn flat_request_json(effects: &str) -> String {
    format!(
        r#"{{
            "schema_version": 1,
            "projectile": {{"mass_kg": 0.01134, "diameter_m": 0.00782, "length_m": 0.0338,
                            "drag_model": "G7", "ballistic_coefficient": 0.243}},
            "rifle": {{"muzzle_velocity_mps": 823.0, "sight_height_m": 0.0381}},
            "shot": {{"max_range_m": 900.0, "zero_distance_m": 100.0}},
            "atmosphere": {{"altitude_m": 0.0, "temperature_k": 288.15, "pressure_pa": 101325.0,
                            "relative_humidity": 0.5}},
            "wind": {{"speed_mps": 4.4704, "direction_from_rad": 1.5707963267948966}},
            "solver": {{}},
            "effects": {{{effects}}},
            "sampling": {{"interval_m": 100.0}}
        }}"#
    )
}

/// A lofted 2500 m shot whose apex is ~186 m above the muzzle, well clear of the 10 m
/// reference height, so the boundary-layer profile is genuinely engaged. `muzzle_angle_rad`
/// is supplied directly (no zero search) to keep the launch angle identical across models —
/// otherwise a re-converged zero would confound the comparison.
fn lofted_request_json(effects: &str) -> String {
    format!(
        r#"{{
            "schema_version": 1,
            "projectile": {{"mass_kg": 0.01134, "diameter_m": 0.00782, "length_m": 0.0338,
                            "drag_model": "G7", "ballistic_coefficient": 0.243}},
            "rifle": {{"muzzle_velocity_mps": 823.0, "sight_height_m": 0.0381}},
            "shot": {{"max_range_m": 2500.0, "muzzle_angle_rad": 0.15,
                      "ground_threshold_m": -2000.0}},
            "atmosphere": {{"altitude_m": 0.0, "temperature_k": 288.15, "pressure_pa": 101325.0,
                            "relative_humidity": 0.5}},
            "wind": {{"speed_mps": 4.4704, "direction_from_rad": 1.5707963267948966}},
            "solver": {{}},
            "effects": {{{effects}}},
            "sampling": {{"interval_m": 250.0}}
        }}"#
    )
}

fn solve(json: &str) -> SolveSuccessV1 {
    let request = decode_solve_request_v1(json).expect("request must decode");
    solve_v1(request).expect("request must solve")
}

fn model_clause(model: &str) -> String {
    format!(r#""wind_shear_model": "{model}""#)
}

// (a) + (b)
#[test]
fn omitted_field_is_byte_identical_and_explicit_none_solves_the_same() {
    let omitted = solve(&flat_request_json(""));

    // The echo is absent, and so is any trace of the key anywhere in the envelope: a
    // response to a request written before this field existed serializes exactly as it did.
    assert!(omitted.resolved_request.effects.wind_shear_model.is_none());
    let omitted_wire = serde_json::to_string(&omitted).unwrap();
    assert!(
        !omitted_wire.contains("wind_shear_model"),
        "omitting the field must leave the whole envelope byte-identical to pre-0.36.0, \
         including the resolved echo"
    );

    // An explicit "none" differs from omission ONLY by that echo. This is the assertion that
    // would fail if `enable_wind_shear` were set true alongside a "none" model name, because
    // cli_api's shear branch would then quietly run the power law.
    let mut explicit_none = solve(&flat_request_json(&model_clause("none")));
    assert_eq!(
        explicit_none.resolved_request.effects.wind_shear_model,
        Some(WindShearModelV1::None)
    );
    explicit_none.resolved_request.effects.wind_shear_model = None;
    assert_eq!(
        serde_json::to_string(&explicit_none).unwrap(),
        omitted_wire,
        "explicit \"none\" must solve identically to omitting the field, differing only in \
         the resolved echo of the field itself"
    );
}

// (c)
#[test]
fn every_accepted_model_is_applied_and_echoed_canonically() {
    for (requested, expected) in [
        ("none", WindShearModelV1::None),
        ("logarithmic", WindShearModelV1::Logarithmic),
        ("power_law", WindShearModelV1::PowerLaw),
        ("ekman_spiral", WindShearModelV1::EkmanSpiral),
        // Alias in, canonical spelling out.
        ("ekman", WindShearModelV1::EkmanSpiral),
    ] {
        let success = solve(&flat_request_json(&model_clause(requested)));
        assert_eq!(
            success.resolved_request.effects.wind_shear_model,
            Some(expected),
            "`{requested}` must be accepted and echoed as `{}`",
            expected.as_engine_str()
        );
        // The echo is what a front end reads to tell the user what took effect, so it has to
        // survive serialization as the canonical name too.
        let wire = serde_json::to_value(&success.resolved_request).unwrap();
        assert_eq!(
            wire["effects"]["wind_shear_model"],
            serde_json::Value::from(expected.as_engine_str()),
            "`{requested}` must serialize as its canonical spelling"
        );
    }
}

/// The Ekman spiral is accepted for parity with the CLI and the engine's own model names, but
/// this solve path's boundary-layer evaluator has no near-ground closed form for it and leaves
/// the wind at the operative value. That is exactly the "applied but inert" state this feature
/// exists to make visible, so it must carry a warning rather than pass quietly — and no other
/// accepted model may carry it.
#[test]
fn ekman_spiral_is_warned_about_rather_than_silently_inert() {
    let ekman = solve(&lofted_request_json(&model_clause("ekman_spiral")));
    assert!(
        ekman
            .warnings
            .iter()
            .any(|w| w.code == "wind_shear_model_not_modeled"
                && w.path.as_deref() == Some("$.effects.wind_shear_model")),
        "ekman_spiral must warn that it leaves the wind unchanged: {:?}",
        ekman.warnings
    );

    let baseline = solve(&lofted_request_json(&model_clause("none")));
    assert_eq!(
        ekman.samples, baseline.samples,
        "the warning must be telling the truth: ekman_spiral leaves the trajectory unchanged"
    );

    for model in ["none", "logarithmic", "power_law"] {
        let success = solve(&lofted_request_json(&model_clause(model)));
        assert!(
            !success
                .warnings
                .iter()
                .any(|w| w.code == "wind_shear_model_not_modeled"),
            "`{model}` is modelled on this path and must not carry the not-modeled warning"
        );
    }
}

// (d) The point of the whole feature: the numbers move.
#[test]
fn shear_moves_the_trajectory_on_a_lofted_shot() {
    let windage_at_terminal = |effects: &str| -> f64 {
        let success = solve(&lofted_request_json(effects));
        success.samples.last().expect("samples").windage_m
    };

    let unsheared = windage_at_terminal("");
    let power_law = windage_at_terminal(&model_clause("power_law"));
    let logarithmic = windage_at_terminal(&model_clause("logarithmic"));

    // Shear only ever increases the operative wind (the profile ratio is floored at 1.0), and
    // this crosswind pushes windage negative, so both sheared solves must be further from
    // zero than the unsheared one.
    for (name, sheared) in [("power_law", power_law), ("logarithmic", logarithmic)] {
        let shift = (sheared - unsheared).abs();
        assert!(
            shift > 1.0,
            "{name} shear must move terminal windage by more than a metre on a lofted shot, \
             got {sheared} vs unsheared {unsheared} (shift {shift} m)"
        );
        assert!(
            sheared.abs() > unsheared.abs(),
            "{name} shear can only increase the operative wind: {sheared} vs {unsheared}"
        );
    }

    // The two profiles are different physics and must not collapse onto one another — a
    // single shared code path that ignored the model name would pass everything above.
    assert!(
        (power_law - logarithmic).abs() > 1.0e-6,
        "the logarithmic and power-law profiles must produce different answers, got \
         {power_law} and {logarithmic}"
    );

    // A flat shot stays inside the boundary layer's floored region, so the same models are
    // correctly inert there. Pinned so a future change to the profile cannot start perturbing
    // ordinary flat-fire solves unnoticed.
    let flat_unsheared = solve(&flat_request_json("")).samples;
    let flat_sheared = solve(&flat_request_json(&model_clause("power_law"))).samples;
    assert_eq!(
        flat_unsheared, flat_sheared,
        "below the 10 m reference height the profile is floored at 1.0, so a flat-fire solve \
         must be unaffected"
    );
}

// (e)
#[test]
fn unknown_models_are_structured_errors_not_silent_fallbacks() {
    for bogus in [
        "gibberish",
        "",
        "POWER_LAW",
        "power-law",
        // A real engine model deliberately left off the wire: it needs a caller-supplied
        // layer table this contract has no field for, and would degrade to the surface wind.
        "custom_layers",
    ] {
        let error = decode_solve_request_v1(&flat_request_json(&model_clause(bogus)))
            .expect_err("an unknown wind_shear_model must be rejected");
        let rendered = serde_json::to_string(&error).unwrap();
        assert!(
            rendered.contains("invalid_value"),
            "`{bogus}` must be an invalid_value: {rendered}"
        );
        assert!(
            rendered.contains("$.effects.wind_shear_model"),
            "`{bogus}` must be reported at the exact request path: {rendered}"
        );
        assert!(
            rendered.contains("power_law") && rendered.contains("logarithmic"),
            "the error must list the accepted spellings so a typo is self-correcting: \
             {rendered}"
        );
    }

    // Wrong JSON type, same treatment.
    let error = decode_solve_request_v1(&flat_request_json(r#""wind_shear_model": true"#))
        .expect_err("a non-string wind_shear_model must be rejected");
    assert!(serde_json::to_string(&error)
        .unwrap()
        .contains("$.effects.wind_shear_model"));
}

// (f)
#[test]
fn the_model_survives_a_resolved_request_round_trip() {
    for model in ["none", "logarithmic", "power_law", "ekman_spiral"] {
        let first = solve(&lofted_request_json(&model_clause(model)));
        let rebuilt: SolveRequestV1 = (&first.resolved_request).into();
        assert_eq!(
            rebuilt.effects.wind_shear_model,
            first.resolved_request.effects.wind_shear_model,
            "`{model}` must be carried onto the rebuilt request"
        );

        let second = solve_v1(rebuilt).expect("the rebuilt request must solve");
        assert_eq!(
            serde_json::to_value(&first.resolved_request).unwrap(),
            serde_json::to_value(&second.resolved_request).unwrap(),
            "`{model}`: the resolved request changed after a round trip"
        );
        assert_eq!(
            first.samples, second.samples,
            "`{model}`: re-solving a round-tripped request must reproduce the trajectory, \
             not silently drop the shear"
        );
    }

    // An omitted model round-trips as still omitted, rather than materializing an echo.
    let first = solve(&flat_request_json(""));
    let rebuilt: SolveRequestV1 = (&first.resolved_request).into();
    assert!(rebuilt.effects.wind_shear_model.is_none());
    let second = solve_v1(rebuilt).expect("solve");
    assert!(second
        .resolved_request
        .effects
        .wind_shear_model
        .is_none());
}

// (g)
#[test]
fn segmented_wind_and_shear_are_a_conflict_not_a_silently_dropped_model() {
    let segmented = |effects: &str| -> String {
        format!(
            r#"{{
                "schema_version": 1,
                "projectile": {{"mass_kg": 0.01134, "diameter_m": 0.00782, "length_m": 0.0338,
                                "drag_model": "G7", "ballistic_coefficient": 0.243}},
                "rifle": {{"muzzle_velocity_mps": 823.0, "sight_height_m": 0.0381}},
                "shot": {{"max_range_m": 900.0, "zero_distance_m": 100.0}},
                "atmosphere": {{}},
                "wind": {{"segments": [
                    {{"until_distance_m": 450.0, "speed_mps": 3.0, "direction_from_rad": 1.57}},
                    {{"until_distance_m": 900.0, "speed_mps": 6.0, "direction_from_rad": 1.57}}
                ]}},
                "solver": {{}},
                "effects": {{{effects}}},
                "sampling": {{"interval_m": 100.0}}
            }}"#
        )
    };

    for model in ["logarithmic", "power_law", "ekman_spiral"] {
        let request = decode_solve_request_v1(&segmented(&model_clause(model)))
            .expect("the pair is structurally valid; it is refused semantically");
        let error = solve_v1(request).expect_err("segments plus shear must be refused");
        let rendered = serde_json::to_string(&error).unwrap();
        assert!(
            rendered.contains("conflicting_fields")
                && rendered.contains("$.effects.wind_shear_model"),
            "`{model}` with segments must be a conflicting_fields error at the model's path: \
             {rendered}"
        );
    }

    // Segments with no shear asked for stay valid, both omitted and explicitly "none".
    for effects in [String::new(), model_clause("none")] {
        let request = decode_solve_request_v1(&segmented(&effects)).expect("decode");
        solve_v1(request).expect("segmented wind without shear must still solve");
    }
}
