use serde_json::Value;
use std::collections::BTreeSet;
use std::process::Command;

const ABS_TOLERANCE: f64 = 1.0e-10;
const REL_TOLERANCE: f64 = 1.0e-9;

fn run_legacy_trajectory(full: bool) -> Value {
    let mut command = Command::new(env!("CARGO_BIN_EXE_ballistics"));
    command.args([
        "trajectory",
        "--velocity",
        "2700",
        "--bc",
        "0.475",
        "--mass",
        "168",
        "--diameter",
        "0.308",
        "--max-range",
        "25",
        "--ignore-ground-impact",
        "--output",
        "json",
    ]);
    if full {
        command.arg("--full");
    }

    let output = command
        .output()
        .expect("run legacy trajectory JSON command");
    assert!(
        output.status.success(),
        "legacy trajectory command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap_or_else(|error| {
        panic!(
            "legacy trajectory stdout is not JSON: {error}; stdout={}",
            String::from_utf8_lossy(&output.stdout)
        )
    })
}

fn object_keys(value: &Value) -> BTreeSet<String> {
    value
        .as_object()
        .expect("fixture value must be an object")
        .keys()
        .cloned()
        .collect()
}

fn expected_top_level_keys() -> BTreeSet<String> {
    [
        "units",
        "max_range",
        "max_height",
        "time_of_flight",
        "impact_velocity",
        "impact_energy",
        "stability_coefficient",
        "spin_drift",
        "trajectory",
    ]
    .into_iter()
    .map(str::to_owned)
    .collect()
}

fn expected_point_keys() -> BTreeSet<String> {
    ["time", "x", "y", "z", "velocity", "energy"]
        .into_iter()
        .map(str::to_owned)
        .collect()
}

fn finite_number(value: &Value, path: &str) -> f64 {
    let number = value
        .as_f64()
        .unwrap_or_else(|| panic!("{path} must remain a JSON number"));
    assert!(number.is_finite(), "{path} must remain finite");
    number
}

fn assert_close(actual: f64, expected: f64, path: &str) {
    let difference = (actual - expected).abs();
    let limit = ABS_TOLERANCE + REL_TOLERANCE * actual.abs().max(expected.abs());
    assert!(
        difference <= limit,
        "{path}: actual {actual:.17e}, expected {expected:.17e}, difference {difference:.17e} exceeds {limit:.17e}"
    );
}

fn assert_summary_shape(document: &Value) {
    assert_eq!(object_keys(document), expected_top_level_keys());
    assert_eq!(document["units"], "imperial");
    for field in [
        "max_range",
        "max_height",
        "time_of_flight",
        "impact_velocity",
        "impact_energy",
        "stability_coefficient",
    ] {
        finite_number(&document[field], &format!("$.{field}"));
    }
    assert!(
        document["spin_drift"].is_null(),
        "disabled legacy spin drift must remain JSON null"
    );
    assert!(document["trajectory"].is_array());
}

#[test]
fn legacy_trajectory_json_summary_shape_is_unchanged_without_full() {
    let document = run_legacy_trajectory(false);
    assert_summary_shape(&document);
    assert!(
        document["trajectory"]
            .as_array()
            .expect("trajectory array")
            .is_empty(),
        "legacy JSON without --full must retain an empty trajectory array"
    );
}

#[test]
fn legacy_trajectory_json_full_adds_points_without_changing_summary_shape() {
    let summary = run_legacy_trajectory(false);
    let full = run_legacy_trajectory(true);
    assert_summary_shape(&summary);
    assert_summary_shape(&full);

    for field in [
        "max_range",
        "max_height",
        "time_of_flight",
        "impact_velocity",
        "impact_energy",
        "stability_coefficient",
    ] {
        assert_close(
            finite_number(&full[field], &format!("$.{field}")),
            finite_number(&summary[field], &format!("$.{field}")),
            &format!("$.{field}"),
        );
    }

    let points = full["trajectory"]
        .as_array()
        .expect("--full trajectory array");
    assert!(!points.is_empty(), "--full must retain trajectory points");
    for (index, point) in points.iter().enumerate() {
        assert_eq!(
            object_keys(point),
            expected_point_keys(),
            "$.trajectory[{index}] shape changed"
        );
        for field in ["time", "x", "y", "z", "velocity", "energy"] {
            finite_number(&point[field], &format!("$.trajectory[{index}].{field}"));
        }
    }
}
