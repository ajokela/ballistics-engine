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
    let v = call("true.fit", json!({
        "model": model(),
        "drop_unit": "mil",
        "observations": [
            {"range_yd": 400.0, "drop": 1.9, "sigma": 0.1},
            {"range_yd": 600.0, "drop": 4.3, "sigma": 0.1},
            {"range_yd": 800.0, "drop": 7.8, "sigma": 0.1}
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
    assert!(!r["approximation"].is_null(), "approximation missing from {r}");
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
