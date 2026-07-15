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
        // MBA-1315: self-describing units/axes metadata, additive.
        "legend",
    ]
    .into_iter()
    .map(str::to_owned)
    .collect()
}

fn expected_legend_units_keys() -> BTreeSet<String> {
    ["distance", "velocity", "energy"]
        .into_iter()
        .map(str::to_owned)
        .collect()
}

fn expected_legend_axes_keys() -> BTreeSet<String> {
    ["x", "y", "z"].into_iter().map(str::to_owned).collect()
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

fn assert_legend_shape(document: &Value) {
    let legend = &document["legend"];
    let expected_legend_keys: BTreeSet<String> = ["units", "axes"]
        .into_iter()
        .map(str::to_owned)
        .collect();
    assert_eq!(object_keys(legend), expected_legend_keys);

    let units = &legend["units"];
    assert_eq!(object_keys(units), expected_legend_units_keys());
    // MBA-1315: this fixture always runs imperial (default --units), so the concrete labels
    // are pinned. A metric run is covered separately below.
    assert_eq!(units["distance"], "yd");
    assert_eq!(units["velocity"], "fps");
    assert_eq!(units["energy"], "ft-lb");

    let axes = &legend["axes"];
    assert_eq!(object_keys(axes), expected_legend_axes_keys());
    for axis in ["x", "y", "z"] {
        let text = axes[axis]
            .as_str()
            .unwrap_or_else(|| panic!("$.legend.axes.{axis} must be a string"));
        assert!(!text.is_empty(), "$.legend.axes.{axis} must not be empty");
    }
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
    assert_legend_shape(document);
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

/// MBA-1315: the `legend.units` labels must track `--units`, not just describe the imperial
/// default pinned by [`assert_legend_shape`].
#[test]
fn legacy_trajectory_json_legend_units_track_metric() {
    let mut command = Command::new(env!("CARGO_BIN_EXE_ballistics"));
    command.args([
        "--units",
        "metric",
        "trajectory",
        "--velocity",
        "823",
        "--bc",
        "0.475",
        "--mass",
        "10.9",
        "--diameter",
        "7.82",
        "--max-range",
        "25",
        "--ignore-ground-impact",
        "--output",
        "json",
    ]);
    let output = command
        .output()
        .expect("run metric legacy trajectory JSON command");
    assert!(
        output.status.success(),
        "metric legacy trajectory command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let document: Value = serde_json::from_slice(&output.stdout).unwrap_or_else(|error| {
        panic!(
            "metric legacy trajectory stdout is not JSON: {error}; stdout={}",
            String::from_utf8_lossy(&output.stdout)
        )
    });

    assert_eq!(document["units"], "metric");
    let legend_units = &document["legend"]["units"];
    assert_eq!(legend_units["distance"], "m");
    assert_eq!(legend_units["velocity"], "m/s");
    assert_eq!(legend_units["energy"], "J");
}

fn run_full_trajectory_with_wind(wind_speed_mph: &str, wind_direction_deg: &str) -> Value {
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
        "100",
        "--time-step",
        "0.001",
        "--wind-speed",
        wind_speed_mph,
        "--wind-direction",
        wind_direction_deg,
        "--full",
        "--output",
        "json",
    ]);
    let output = command
        .output()
        .expect("run windy legacy trajectory JSON command");
    assert!(
        output.status.success(),
        "windy legacy trajectory command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap_or_else(|error| {
        panic!(
            "windy legacy trajectory stdout is not JSON: {error}; stdout={}",
            String::from_utf8_lossy(&output.stdout)
        )
    })
}

/// MBA-1315 axis-doc-vs-behavior: `legend.axes.x` must describe the SAME sign convention the
/// engine actually produces, not an assumed one. Cross-check the documented wording against
/// controlled crosswind runs — wind FROM the left (`--wind-direction 270`, per the CLI's own
/// "270=from left" help text) must drift `x` positive, and FROM the right (`90`) must drift
/// it negative, exactly as `legend.axes.x` states.
#[test]
fn legacy_trajectory_json_axes_text_matches_observed_crosswind_sign() {
    let baseline = run_full_trajectory_with_wind("0", "0");
    let axis_x = baseline["legend"]["axes"]["x"]
        .as_str()
        .expect("legend.axes.x string");
    assert!(
        axis_x.contains("shooter's right") && axis_x.contains("positive"),
        "legend.axes.x wording changed; update this test's sign assumptions: {axis_x}"
    );

    fn last_x(document: &Value) -> f64 {
        document["trajectory"]
            .as_array()
            .expect("trajectory array")
            .last()
            .expect("at least one trajectory point")["x"]
            .as_f64()
            .expect("x is a number")
    }

    assert_eq!(
        last_x(&baseline),
        0.0,
        "a no-wind run must have zero lateral drift"
    );

    let from_left = run_full_trajectory_with_wind("10", "270");
    let from_right = run_full_trajectory_with_wind("10", "90");
    let x_from_left = last_x(&from_left);
    let x_from_right = last_x(&from_right);

    assert!(
        x_from_left > 0.0,
        "wind FROM the left (--wind-direction 270) must drift x positive per legend.axes.x, got {x_from_left}"
    );
    assert!(
        x_from_right < 0.0,
        "wind FROM the right (--wind-direction 90) must drift x negative per legend.axes.x, got {x_from_right}"
    );
}
