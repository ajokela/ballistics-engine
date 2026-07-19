//! Native CLI contracts for MBA-1346/MBA-1353.

use std::process::Command;

use serde_json::Value;

fn binary() -> &'static str {
    env!("CARGO_BIN_EXE_ballistics")
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(binary())
        .args(args)
        .output()
        .expect("run ballistics CLI")
}

#[test]
fn plan_truing_json_is_constrained_deterministic_and_self_describing() {
    let args = [
        "plan-truing",
        "--velocity",
        "2700",
        "--bc",
        "0.475",
        "--drag-model",
        "g1",
        "--mass",
        "168",
        "--diameter",
        "0.308",
        "--candidate-ranges",
        "200,300,400,500,600,700,800,900",
        "--observation-count",
        "3",
        "--minimum-separation",
        "100",
        "--measurement-resolution",
        "0.03",
        "--output",
        "json",
    ];
    let first = run(&args);
    assert!(
        first.status.success(),
        "{}",
        String::from_utf8_lossy(&first.stderr)
    );
    let second = run(&args);
    assert!(second.status.success());
    assert_eq!(
        first.stdout, second.stdout,
        "planner output must be deterministic"
    );

    let value: Value = serde_json::from_slice(&first.stdout).expect("planner JSON");
    assert_eq!(value["schema_version"], "truing-experiment-plan-v1");
    assert_eq!(value["mode"], "joint_mv_bc");
    assert_eq!(value["measurement_model"]["sigma_1sd"], 0.03);
    assert_eq!(value["measurement_model"]["drop_unit"], "mil");
    assert!(
        value["information"]["minimum_singular_value"]
            .as_f64()
            .unwrap()
            > 0.0
    );
    assert!(value["information"]["weak_axis_fractional_sigma_1sd"].is_number());
    let selected = value["selected_stations"].as_array().unwrap();
    assert_eq!(selected.len(), 3);
    let supplied = [200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0];
    let mut ranges: Vec<f64> = selected
        .iter()
        .map(|station| station["range_yd"].as_f64().unwrap())
        .collect();
    assert!(ranges.iter().all(|range| supplied.contains(range)));
    ranges.sort_by(f64::total_cmp);
    assert!(ranges.windows(2).all(|pair| pair[1] - pair[0] >= 100.0));
    assert!(selected
        .iter()
        .all(|station| station["leave_one_out_information_gain_nats"].is_number()));
}

#[test]
fn short_facility_is_explicitly_mv_only() {
    let output = run(&[
        "plan-truing",
        "-v",
        "2700",
        "-b",
        "0.475",
        "-m",
        "168",
        "-d",
        "0.308",
        "--candidate-ranges",
        "200,250,300",
        "--observation-count",
        "3",
        "--measurement-resolution",
        "0.03",
        "--output",
        "json",
    ]);
    assert!(output.status.success());
    let value: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["mode"], "mv_only");
    assert!(value["warnings"]
        .as_array()
        .unwrap()
        .iter()
        .any(|warning| warning
            .as_str()
            .unwrap_or("")
            .contains("do not reliably separate")));
}

#[test]
fn uncertainty_json_reports_intervals_correlation_and_two_prediction_bands() {
    // Exact forward-model drops for 2700 fps / G1 0.475, produced by the
    // planner's shared nominal solver at 500/600/900 yd.
    let output = run(&[
        "true-velocity",
        "--range",
        "500",
        "--measured-drop",
        "3.1790310753",
        "--observed",
        "600:4.3487384136:0.03",
        "--observed",
        "900:8.8913106904:0.02",
        "--observation-sigma",
        "0.03",
        "--bc",
        "0.45",
        "--drag-model",
        "g1",
        "--mass",
        "168",
        "--diameter",
        "0.308",
        "--predict-range",
        "1000",
        "--prediction-sigma",
        "0.03",
        "--output",
        "json",
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let value: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["schema_version"], "truing-uncertainty-v1");
    assert_eq!(value["approximation"]["available"], true);
    assert!(value["diagnostics"]["map_convergence_criterion"].is_string());
    assert!(value["diagnostics"]["map_scaled_gradient_inf_norm"].is_number());
    assert!(value["diagnostics"]["map_objective_poll_evaluations"].is_number());
    let approximation = &value["approximation"];
    let mv = &approximation["muzzle_velocity"];
    let bc = &approximation["ballistic_coefficient"];
    assert!(mv["lower"].as_f64().unwrap() < 2700.0);
    assert!(mv["upper"].as_f64().unwrap() > 2700.0);
    assert!(bc["lower"].as_f64().unwrap() < 0.475);
    assert!(bc["upper"].as_f64().unwrap() > 0.475);
    assert!(approximation["mv_bc_correlation"].as_f64().unwrap().abs() <= 1.0);
    let band = &value["predictive_bands"][0];
    let latent = &band["latent_model_drop_interval"];
    let future = &band["future_observation_interval"];
    assert!(
        future["standard_deviation"].as_f64().unwrap()
            > latent["standard_deviation"].as_f64().unwrap()
    );
}

#[test]
fn no_uncertainty_flags_keep_the_legacy_schema() {
    let output = run(&[
        "true-velocity",
        "--range",
        "300",
        "--measured-drop",
        "1.3",
        "--observed",
        "600:4.4",
        "--observed",
        "900:9.0",
        "--bc",
        "0.45",
        "--mass",
        "168",
        "--diameter",
        "0.308",
        "--output",
        "json",
    ]);
    assert!(output.status.success());
    let value: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(value.get("schema_version").is_none());
    assert!(value.get("approximation").is_none());
    assert!(value.get("bc_fitted").is_some());
}

#[test]
fn prediction_sigma_requires_a_prediction_range() {
    let output = run(&[
        "true-velocity",
        "--range",
        "500",
        "--measured-drop",
        "3.1790310753",
        "--observed",
        "900:8.8913106904",
        "--observation-sigma",
        "0.03",
        "--prediction-sigma",
        "0.03",
        "--bc",
        "0.475",
        "--mass",
        "168",
        "--diameter",
        "0.308",
    ]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr)
        .contains("--prediction-sigma requires at least one --predict-range"));
}

#[test]
fn planner_csv_includes_measurement_and_information_metadata() {
    let output = run(&[
        "plan-truing",
        "-v",
        "2700",
        "-b",
        "0.475",
        "-m",
        "168",
        "-d",
        "0.308",
        "--candidate-ranges",
        "300,600,900",
        "--observation-count",
        "3",
        "--measurement-resolution",
        "0.03",
        "--output",
        "csv",
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let csv = String::from_utf8_lossy(&output.stdout);
    assert!(csv.contains("measurement_sigma_1sd,measurement_drop_unit"));
    assert!(csv.contains("0.03000000,mil"));
    assert!(csv.contains("minimum_singular_value"));
    assert!(csv.contains("warning"));
}

#[test]
fn uncertainty_csv_includes_intervals_covariance_and_correlation() {
    let output = run(&[
        "true-velocity",
        "--range",
        "500",
        "--measured-drop",
        "3.1790310753",
        "--observed",
        "600:4.3487384136:0.03",
        "--observed",
        "900:8.8913106904:0.02",
        "--observation-sigma",
        "0.03",
        "--bc",
        "0.45",
        "--mass",
        "168",
        "--diameter",
        "0.308",
        "--output",
        "csv",
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let csv = String::from_utf8_lossy(&output.stdout);
    assert!(csv.contains("parameter,estimate,standard_deviation,interval_95_lower"));
    assert!(csv.contains("muzzle_velocity,"));
    assert!(csv.contains("ballistic_coefficient,"));
    assert!(csv.contains("mv_bc_covariance_fps"));
    assert!(csv.contains("mv_bc_correlation"));
    assert!(csv.contains("prior_parameter,mean,sigma_1sd,unit"));
}
