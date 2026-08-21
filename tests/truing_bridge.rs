use ballistics_engine::bridge::bridge_call;
use serde_json::{json, Value};

fn call(command: &str, request: Value) -> Value {
    let envelope = json!({"api_version": 1, "command": command, "request": request});
    serde_json::from_str(&bridge_call(&envelope.to_string())).unwrap()
}

fn model() -> Value {
    json!({
        "muzzle_velocity_fps": 2700.0, "ballistic_coefficient": 0.243,
        "drag_model": "g7", "mass_gr": 168.0, "diameter_in": 0.308,
        "zero_distance_yd": 100.0, "sight_height_in": 2.0,
        "temperature_f": 59.0, "pressure_inhg": 29.92,
        "humidity_pct": 50.0, "altitude_ft": 0.0
    })
}

#[test]
fn true_fit_is_advertised_in_capabilities() {
    let v = call("meta.capabilities", json!(null));
    let names: Vec<&str> = v["result"]["commands"]
        .as_array().unwrap().iter().map(|c| c.as_str().unwrap()).collect();
    assert!(names.contains(&"true.fit"), "capabilities were {names:?}");
}

#[test]
fn true_fit_returns_estimates_with_their_uncertainty() {
    // Observations derived from the engine's own solve on the model below:
    // .308 caliber, 168 gr, 0.308 in dia, G7 BC 0.243, 2700 fps, 100 yd zero, ICAO conditions.
    // Drops (mil) at range: 400 yd -> 2.6882, 600 yd -> 4.8626, 800 yd -> 7.6598
    let v = call("true.fit", json!({
        "model": model(),
        "drop_unit": "mil",
        "observations": [
            {"range_yd": 400.0, "drop": 2.6882, "sigma": 0.1},
            {"range_yd": 600.0, "drop": 4.8626, "sigma": 0.1},
            {"range_yd": 800.0, "drop": 7.6598, "sigma": 0.1}
        ],
        "priors": {"muzzle_velocity_fps": {"mean": 2700.0, "sigma": 20.0},
                   "ballistic_coefficient": {"mean": 0.243, "sigma": 0.01}},
        "predictions": []
    }));
    assert_eq!(v["ok"], true, "response was {v}");
    let r = &v["result"];
    assert!(r["map_muzzle_velocity_fps"].as_f64().unwrap() > 0.0);
    assert!(r["map_ballistic_coefficient"].as_f64().unwrap() > 0.0);

    // The honesty invariant: a point estimate never travels without its uncertainty.
    // Assert that the approximation is available and carries valid intervals for both parameters.
    assert_eq!(r["approximation"]["status"], "available",
        "approximation must be available, got: {}", r["approximation"]["status"]);
    let details = &r["approximation"]["details"];

    // Verify that both parameters have finite intervals.
    let mv_lower = details["muzzle_velocity_interval_95"]["lower"].as_f64().unwrap();
    let mv_upper = details["muzzle_velocity_interval_95"]["upper"].as_f64().unwrap();
    assert!(mv_lower.is_finite() && mv_upper.is_finite(),
        "muzzle velocity interval bounds must be finite");
    assert!(mv_lower < mv_upper, "muzzle velocity lower bound must be less than upper");

    let bc_lower = details["ballistic_coefficient_interval_95"]["lower"].as_f64().unwrap();
    let bc_upper = details["ballistic_coefficient_interval_95"]["upper"].as_f64().unwrap();
    assert!(bc_lower.is_finite() && bc_upper.is_finite(),
        "ballistic coefficient interval bounds must be finite");
    assert!(bc_lower < bc_upper, "ballistic coefficient lower bound must be less than upper");
}

#[test]
fn true_fit_rejects_a_null_payload_and_unknown_fields() {
    let v = call("true.fit", json!(null));
    assert_eq!(v["error"]["code"], "invalid_request");

    let v = call("true.fit", json!({
        "model": model(), "drop_unit": "mil",
        "observations": [{"range_yd": 600.0, "drop": 4.3, "sigma": 0.1, "typo": 1}],
        "priors": {"muzzle_velocity_fps": null, "ballistic_coefficient": null},
        "predictions": []
    }));
    assert_eq!(v["error"]["code"], "invalid_request");
}

#[test]
fn true_wind_fits_a_crosswind_from_an_observed_miss() {
    let v = call("true.wind", json!({
        "observations": [{"range_m": 457.2, "miss_right_m": 0.315, "sigma_m": null}],
        "muzzle_velocity_fps": 2700.0, "bc": 0.243, "drag_model": "g7",
        "mass_gr": 168.0, "diameter_in": 0.308, "zero_distance_yd": 100.0,
        "sight_height_in": 2.0, "temperature_f": 59.0, "pressure_inhg": 29.92,
        "humidity_pct": 50.0, "altitude_ft": 0.0,
        "twist": {"rate_in": 11.0, "right_hand": true},
        "earth": null, "called_crosswind_mph": null
    }));
    assert_eq!(v["ok"], true, "response was {v}");
    assert!(v["result"]["mean_crosswind_mph"].is_number());
    assert_eq!(v["result"]["solutions"].as_array().unwrap().len(), 1);
}

#[test]
fn true_wind_is_advertised_in_capabilities() {
    let v = call("meta.capabilities", json!(null));
    let names: Vec<&str> = v["result"]["commands"]
        .as_array().unwrap().iter().map(|c| c.as_str().unwrap()).collect();
    assert!(names.contains(&"true.wind"), "capabilities were {names:?}");
}