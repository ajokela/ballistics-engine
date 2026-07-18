// Unit tests for WASM bindings
#[cfg(test)]
#[cfg(target_arch = "wasm32")]
mod tests {
    use crate::wasm::*;
    use wasm_bindgen_test::*;

    wasm_bindgen_test_configure!(run_in_browser);

    #[wasm_bindgen_test]
    fn test_wasm_ballistics_creation() {
        let wasm = WasmBallistics::new();
        assert!(wasm.run_command("help").is_ok());
    }

    #[wasm_bindgen_test]
    fn test_help_command() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command("help").unwrap();
        assert!(result.contains("Ballistics Engine"));
        assert!(result.contains("trajectory"));
        assert!(result.contains("zero"));
        assert!(result.contains("monte-carlo"));
    }

    #[wasm_bindgen_test]
    fn test_basic_trajectory_command() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308")
            .unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
        assert!(result.contains("Range"));
        assert!(result.contains("Drop"));
        assert!(result.contains("Velocity"));
    }

    #[wasm_bindgen_test]
    fn test_trajectory_with_auto_zero() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --auto-zero 200")
            .unwrap();
        assert!(result.contains("Rifle zeroed at 200 yards"));
        assert!(result.contains("MOA"));
        assert!(result.contains("mrad"));
    }

    #[wasm_bindgen_test]
    fn test_metric_units() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("--units metric trajectory -v 823 -b 0.475 -m 10.9 -d 7.82")
            .unwrap();
        assert!(result.contains("m/s"));
        assert!(result.contains("meters"));
        assert!(!result.contains("yards"));
        assert!(!result.contains("fps"));
    }

    #[wasm_bindgen_test]
    fn test_json_output_format() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 -o json")
            .unwrap();
        assert!(result.contains("\"trajectory\""));
        assert!(result.contains("\"summary\""));
        assert!(result.contains("range_yards"));
        assert!(result.contains("drop_inches"));
        // MBA-1315: self-describing units/axes metadata, additive third top-level key
        // alongside the pre-existing "trajectory"/"summary" keys.
        assert!(result.contains("\"legend\""));
        assert!(result.contains("\"axes\""));
        assert!(result.contains("\"drift\""));
    }

    /// MBA-1315 axis-doc-vs-behavior (WASM parity with the native CLI legacy JSON test of the
    /// same name): `legend.axes.drift` must describe the SAME sign convention the formatter
    /// actually produces, not an assumed one. Wind FROM the left (`--wind-direction 270`)
    /// must drift `drift_inches` positive, and FROM the right (`90`) must drift it negative,
    /// exactly as the legend states -- identical in sign and source (`position.z`) to the
    /// native CLI legacy JSON's `x`.
    #[wasm_bindgen_test]
    fn json_legend_drift_axis_matches_observed_crosswind_sign() {
        let wasm = WasmBallistics::new();
        let out = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 -o json")
            .unwrap();
        let doc: serde_json::Value = serde_json::from_str(&out).unwrap();
        let drift_text = doc["legend"]["axes"]["drift"]
            .as_str()
            .expect("legend.axes.drift string");
        assert!(
            drift_text.contains("shooter's") && drift_text.contains("right"),
            "legend.axes.drift wording changed; update this test's sign assumptions: {drift_text}"
        );
        assert!(drift_text.contains("positive"));

        fn last_drift_inches(json: &str) -> f64 {
            let doc: serde_json::Value = serde_json::from_str(json).unwrap();
            doc["trajectory"]
                .as_array()
                .expect("trajectory array")
                .last()
                .expect("at least one trajectory point")["drift_inches"]
                .as_f64()
                .expect("drift_inches is a number")
        }

        let baseline = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --full -o json")
            .unwrap();
        let from_left = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --wind-speed 10 \
                 --wind-direction 270 --full -o json",
            )
            .unwrap();
        let from_right = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --wind-speed 10 \
                 --wind-direction 90 --full -o json",
            )
            .unwrap();

        assert_eq!(
            last_drift_inches(&baseline),
            0.0,
            "a no-wind run must have zero lateral drift"
        );
        assert!(
            last_drift_inches(&from_left) > 0.0,
            "wind FROM the left (--wind-direction 270) must drift right (positive) per legend.axes.drift"
        );
        assert!(
            last_drift_inches(&from_right) < 0.0,
            "wind FROM the right (--wind-direction 90) must drift left (negative) per legend.axes.drift"
        );
    }

    #[wasm_bindgen_test]
    fn test_csv_output_format() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 -o csv")
            .unwrap();
        assert!(result.contains("Range(yards),Drop(inches)"));
        assert!(result.contains(","));
        assert!(!result.contains("Trajectory Calculation Results")); // CSV shouldn't have the header
    }

    #[wasm_bindgen_test]
    fn test_zero_command() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("zero -v 2700 -b 0.475 -m 168 -d 0.308 --target-distance 300")
            .unwrap();
        assert!(result.contains("Zero Calculation Results"));
        assert!(result.contains("Target Distance"));
        assert!(result.contains("MOA Adjustment"));
        assert!(result.contains("Mrad Adjustment"));
    }

    #[wasm_bindgen_test]
    fn test_monte_carlo_command() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 -n 100")
            .unwrap();
        assert!(result.contains("Monte Carlo Simulation Results"));
        assert!(result.contains("Simulations Run: 100"));
        assert!(result.contains("Range Statistics"));
        assert!(result.contains("Impact Velocity Statistics"));
    }

    #[wasm_bindgen_test]
    fn test_estimate_bc_command() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("estimate-bc -v 2700 -m 168 -d 0.308 --data \"100,2.5;200,10.2;300,23.5\"")
            .unwrap();
        assert!(result.contains("BC Estimation Results"));
        assert!(result.contains("Estimated BC"));
        // Default --drag-model is "both", so a G1 and a G7 row are printed, each fit to the
        // 3-point drop series.
        assert!(result.contains("G1"));
        assert!(result.contains("G7"));
        assert!(result.contains("drop (3 pts)"));
    }

    #[wasm_bindgen_test]
    fn test_estimate_bc_velocity_and_drag_model() {
        let wasm = WasmBallistics::new();
        // G7-only, fit against a velocity-retention series.
        let result = wasm
            .run_command(
                "estimate-bc -v 2650 -m 77 -d 0.224 \
                 --velocity-data \"200,2270;400,1930;600,1610\" --drag-model g7",
            )
            .unwrap();
        assert!(result.contains("BC Estimation Results"));
        assert!(result.contains("G7"));
        assert!(result.contains("velocity (3 pts)"));
        // g7-only must NOT print a G1 row.
        assert!(!result.contains("G1"));
    }

    #[wasm_bindgen_test]
    fn test_advanced_physics_flags() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
             --enable-magnus --enable-coriolis --enable-spin-drift \
             --twist-rate 10 --latitude 45",
            )
            .unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
        // The trajectory should complete without errors when physics flags are enabled
    }

    #[wasm_bindgen_test]
    fn omitted_twist_matches_miller_default_and_explicit_override_wins() {
        let wasm = WasmBallistics::new();
        let final_drift = |output: &str, key: &str| {
            let json: serde_json::Value = serde_json::from_str(output).unwrap();
            json["trajectory"].as_array().unwrap().last().unwrap()[key]
                .as_f64()
                .unwrap()
        };
        let command = "trajectory -v 2750 -b 0.219 -m 77 -d 0.224 --drag-model G7 \
                       --enable-spin-drift --max-range 500 -o json";
        let default_twist = crate::stability::default_twist_inches(
            0.224 * 0.0254,
            77.0 * crate::constants::GRAINS_TO_KG,
            2750.0 * 0.3048,
        );

        let omitted = wasm.run_command(command).unwrap();
        let explicit_default = wasm
            .run_command(&format!("{command} --twist-rate {default_twist:.17}"))
            .unwrap();
        let explicit_twelve = wasm
            .run_command(&format!("{command} --twist-rate 12"))
            .unwrap();

        let omitted_drift = final_drift(&omitted, "drift_inches");
        let default_drift = final_drift(&explicit_default, "drift_inches");
        let twelve_drift = final_drift(&explicit_twelve, "drift_inches");
        assert!((omitted_drift - default_drift).abs() < 1e-9);
        assert!((omitted_drift - twelve_drift).abs() > 0.1);

        let metric_command = "trajectory --units metric -v 838.2 -b 0.219 -m 4.98951607 \
                              -d 5.6896 --drag-model G7 --enable-spin-drift \
                              --max-range 457.2 -o json";
        let metric_default_twist =
            crate::stability::default_twist_inches(5.6896 * 0.001, 4.98951607 * 0.001, 838.2);
        let omitted_metric = wasm.run_command(metric_command).unwrap();
        let explicit_default_metric = wasm
            .run_command(&format!(
                "{metric_command} --twist-rate {:.17}",
                metric_default_twist * 25.4
            ))
            .unwrap();
        assert!(
            (final_drift(&omitted_metric, "drift_cm")
                - final_drift(&explicit_default_metric, "drift_cm"))
            .abs()
                < 1e-9
        );
    }

    #[wasm_bindgen_test]
    fn test_environmental_conditions() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
             --temperature 32 --pressure 25.0 --humidity 80 --altitude 5000 \
             --wind-speed 10 --wind-direction 90",
            )
            .unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn test_full_trajectory_output() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --full --max-range 500")
            .unwrap();
        // With --full flag, we should see more data points
        let lines: Vec<&str> = result.lines().collect();
        let data_lines = lines.iter().filter(|l| l.contains("yd")).count();
        assert!(data_lines > 5); // Should have many data points with --full
    }

    #[wasm_bindgen_test]
    fn test_powder_temperature_sensitivity() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
             --use-powder-sensitivity --powder-temp 90 --powder-temp-sensitivity 1.5",
            )
            .unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn explicit_zero_velocity_ignores_linear_powder_adjustment() {
        let wasm = WasmBallistics::new();
        // Table format: the "Rifle zeroed at ..." banner is a table-only human header
        // (MBA-1294(a) keeps JSON/CSV pure), and this test reads the banner's adjustment line.
        let command = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
                       --auto-zero 300 --zero-velocity 2400 --zero-temperature 20 \
                       --max-range 100";
        let zero_line = |output: &str| output.lines().next().unwrap().to_string();

        let without_linear = zero_line(&wasm.run_command(command).unwrap());
        let with_linear = zero_line(
            &wasm
                .run_command(&format!(
                    "{command} --use-powder-sensitivity --powder-temp-sensitivity 5"
                ))
                .unwrap(),
        );

        assert!(without_linear.starts_with("Rifle zeroed at 300 yards"));
        assert_eq!(with_linear, without_linear);
    }

    #[wasm_bindgen_test]
    fn auto_zero_inherits_explicit_shot_day_powder_temperature() {
        let wasm = WasmBallistics::new();
        // Table format: reads the auto-zero banner's adjustment (banner is table-only
        // since MBA-1294(a); JSON/CSV stay pure machine output).
        let command = "trajectory -v 2650 -b 0.19 -m 77 -d 0.224 --drag-model g7 \
                       --temperature 32 --powder-temp 68 \
                       --powder-temp-curve 32:2650,77:2720 \
                       --auto-zero 100 --max-range 100";
        let adjustment = |output: &str| {
            output
                .lines()
                .next()
                .and_then(|line| line.split_once("(Adjustment: "))
                .map(|(_, value)| value.to_string())
                .expect("successful auto-zero adjustment")
        };

        let inherited = adjustment(&wasm.run_command(command).unwrap());
        let explicit = adjustment(
            &wasm
                .run_command(&format!("{command} --zero-powder-temp 68"))
                .unwrap(),
        );

        assert_eq!(inherited, explicit);

        let zero_air = adjustment(
            &wasm
                .run_command(&format!("{command} --zero-temperature 77"))
                .unwrap(),
        );
        let explicit_zero_air = adjustment(
            &wasm
                .run_command(&format!(
                    "{command} --zero-temperature 77 --zero-powder-temp 77"
                ))
                .unwrap(),
        );

        assert_eq!(zero_air, explicit_zero_air);
    }

    #[wasm_bindgen_test]
    fn metric_default_powder_temperature_represents_70_fahrenheit() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "trajectory --units metric -v 800 -b 0.3 -m 10 -d 7.62 \
                 --temperature 30 --use-powder-sensitivity \
                 --powder-temp-sensitivity 1 --max-range 1 -o json",
            )
            .unwrap();
        let json: serde_json::Value = serde_json::from_str(&result).unwrap();
        let actual_velocity = json["trajectory"][0]["velocity_mps"].as_f64().unwrap();
        let reference_temp_c = (70.0 - 32.0) * 5.0 / 9.0;
        let expected_velocity = 800.0 + (30.0 - reference_temp_c);

        assert!((actual_velocity - expected_velocity).abs() < 1e-9);
    }

    #[wasm_bindgen_test]
    fn test_shooting_angle() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --shooting-angle 15")
            .unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn test_drag_models() {
        let wasm = WasmBallistics::new();

        // Test G1 model
        let g1_result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --drag-model G1")
            .unwrap();
        assert!(g1_result.contains("Trajectory Calculation Results"));

        // Test G7 model
        let g7_result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --drag-model G7")
            .unwrap();
        assert!(g7_result.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn test_invalid_command() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command("invalid-command").unwrap();
        assert!(result.contains("Error: Unknown command"));
        assert!(result.contains("help"));
    }

    #[wasm_bindgen_test]
    fn test_invalid_parameters() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command("trajectory -v not-a-number -b 0.475 -m 168 -d 0.308");
        assert!(result.is_err());
    }

    #[wasm_bindgen_test]
    fn test_metric_monte_carlo() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("--units metric monte-carlo -v 823 -b 0.475 -m 10.9 -d 7.82 -n 50")
            .unwrap();
        assert!(result.contains("m/s"));
        assert!(result.contains("meters"));
        assert!(!result.contains("yards"));
        assert!(!result.contains("fps"));
    }

    #[wasm_bindgen_test]
    fn test_metric_zero_calculation() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "--units metric zero -v 823 -b 0.475 -m 10.9 -d 7.82 --target-distance 200",
            )
            .unwrap();
        assert!(result.contains("200 meters"));
        assert!(!result.contains("yards"));
    }

    #[wasm_bindgen_test]
    fn test_bc_segments() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --use-bc-segments")
            .unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn test_trajectory_sampling() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command(
            "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --sample-trajectory --sample-interval 25"
        ).unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn test_integration_methods() {
        let wasm = WasmBallistics::new();

        // Test RK4 (default)
        let rk4_result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308")
            .unwrap();
        assert!(rk4_result.contains("Trajectory Calculation Results"));

        // Test Euler
        let euler_result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --use-euler")
            .unwrap();
        assert!(euler_result.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn test_estimate_bc_no_data() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("estimate-bc -v 2700 -m 168 -d 0.308")
            .unwrap();
        assert!(result.contains("Error: No data provided"));
        assert!(result.contains("Example"));
    }

    #[wasm_bindgen_test]
    fn test_zero_with_sight_height() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "zero -v 2700 -b 0.475 -m 168 -d 0.308 --target-distance 200 --sight-height 2.5",
            )
            .unwrap();
        assert!(result.contains("Sight Height: 2.5 inches"));
    }

    #[wasm_bindgen_test]
    fn test_monte_carlo_with_variations() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 -n 50 \
             --velocity-std 15 --angle-std 0.2 --bc-std 0.02 \
             --wind-speed-std 3 --wind-dir-std 10",
            )
            .unwrap();
        assert!(result.contains("Monte Carlo Simulation Results"));
        assert!(result.contains("Std Dev"));
    }

    #[wasm_bindgen_test]
    fn test_command_parsing() {
        let wasm = WasmBallistics::new();

        // Test that ballistics prefix is handled
        let result1 = wasm
            .run_command("ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308")
            .unwrap();
        assert!(result1.contains("Trajectory Calculation Results"));

        // Test without prefix
        let result2 = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308")
            .unwrap();
        assert!(result2.contains("Trajectory Calculation Results"));

        // Test with ./ prefix
        let result3 = wasm
            .run_command("./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308")
            .unwrap();
        assert!(result3.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn test_empty_command() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command("").unwrap();
        assert!(result.contains("Ballistics Engine"));
    }

    #[wasm_bindgen_test]
    fn test_help_variants() {
        let wasm = WasmBallistics::new();

        let help1 = wasm.run_command("help").unwrap();
        let help2 = wasm.run_command("--help").unwrap();
        let help3 = wasm.run_command("-h").unwrap();

        assert!(help1.contains("Ballistics Engine"));
        assert!(help2.contains("Ballistics Engine"));
        assert!(help3.contains("Ballistics Engine"));
    }

    #[wasm_bindgen_test]
    fn ignore_ground_impact_reaches_full_range_where_default_truncates() {
        // Reads the number out of the table's trailing "Max Range: <N> yards" line —
        // more robust than the table's own "Bullet struck ground" heuristic, which only
        // fires when the last recorded point lands within 0.01m of the ground and can
        // miss a truncation that stopped a hair above that (as it does for these inputs).
        fn parse_max_range_yards(output: &str) -> f64 {
            output
                .lines()
                .find_map(|l| l.strip_prefix("Max Range: "))
                .and_then(|rest| rest.split_whitespace().next())
                .and_then(|n| n.parse::<f64>().ok())
                .expect("Max Range line present in trajectory output")
        }

        let wasm = WasmBallistics::new();
        // Flat-fire (no auto-zero) drop at 1000 yd far exceeds the default 60in bore
        // height, so the default run truncates at ground impact well short of --max-range.
        let default_result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 1000")
            .unwrap();
        let default_max_range = parse_max_range_yards(&default_result);
        assert!(
            default_max_range < 900.0,
            "expected the default run to truncate well short of 1000 yd, got {default_max_range}"
        );

        // --ignore-ground-impact (bero-feedback 0.23.0) disables that early termination,
        // so the same run reaches the full requested range instead.
        let ignore_result = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 1000 \
                 --ignore-ground-impact",
            )
            .unwrap();
        let ignore_max_range = parse_max_range_yards(&ignore_result);
        assert!(
            ignore_max_range >= 999.0,
            "expected --ignore-ground-impact to reach ~1000 yd, got {ignore_max_range}"
        );
    }

    #[wasm_bindgen_test]
    fn muzzle_height_warning_appears_only_above_the_1000m_threshold() {
        let wasm = WasmBallistics::new();
        // ~25 km "muzzle height" (defeating ground truncation the wrong way) must warn.
        let warned = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --muzzle-height 1000000")
            .unwrap();
        assert!(warned.starts_with("WARNING"));
        assert!(warned.contains("--altitude"));
        assert!(warned.contains("--ignore-ground-impact"));

        // An ordinary elevated stand (2400 in ~= 61 m) must not trigger the warning.
        let unwarned = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --muzzle-height 2400")
            .unwrap();
        assert!(!unwarned.contains("WARNING"));
    }

    #[wasm_bindgen_test]
    fn lead_output_json_parses_with_sane_intercept_range() {
        let wasm = WasmBallistics::new();
        let output = wasm
            .run_command("lead --target-speed 10 --target-angle 0 --range 500 -o json")
            .unwrap();
        let json: serde_json::Value = serde_json::from_str(&output).unwrap();
        let range = json["range"].as_f64().unwrap();
        let intercept_range = json["intercept_range"].as_f64().unwrap();
        // angle 0 = directly away (outbound): the target recedes during the time of
        // flight, so the corrected intercept range must be at or beyond the initial range.
        assert!(intercept_range >= range);
        assert_eq!(json["distance_unit"], "yd");
        assert_eq!(json["adjustment_unit"], "mil");
    }

    #[wasm_bindgen_test]
    fn test_all_physics_flags_combined() {
        let wasm = WasmBallistics::new();

        // Test all physics flags together
        let result = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
             --enable-magnus --enable-coriolis --enable-spin-drift \
             --enable-wind-shear --enable-pitch-damping --enable-precession \
             --sample-trajectory --use-bc-segments --use-powder-sensitivity \
             --twist-rate 10 --latitude 45",
            )
            .unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
    }

    /// MBA-1297 (field report, Bero): at 90 degrees of cant the vertical and
    /// lateral misses at the zero distance MUST be equal — the same rotation
    /// moves both the bore offset and the zero tilt. The WASM formatters used
    /// points[0].y + sight_height as the LOS, double-counting the canted
    /// bore's vertical rise (reported drop = true miss + one sight height).
    #[wasm_bindgen_test]
    fn canted_90_vertical_and_lateral_misses_are_symmetric() {
        let wasm = WasmBallistics::new();
        let out = wasm
            .run_command(
                "trajectory -b 0.19 -m 77 -d 0.224 --drag-model g7 --max-range 200 \
                 --sight-height 2.48 --ignore-ground-impact -v 2650 --auto-zero 100 \
                 --cant 90 --full -o csv",
            )
            .unwrap();
        // CSV rows: distance,drop_in,drift_in,... — find the row nearest 100 yd
        // and the first (muzzle) row.
        let mut muzzle_drop: Option<f64> = None;
        let mut drop100: Option<f64> = None;
        let mut drift100: Option<f64> = None;
        for line in out.lines() {
            let cols: Vec<&str> = line.split(',').collect();
            if cols.len() < 3 {
                continue;
            }
            if let (Ok(dist), Ok(drop), Ok(drift)) = (
                cols[0].trim().parse::<f64>(),
                cols[1].trim().parse::<f64>(),
                cols[2].trim().parse::<f64>(),
            ) {
                if dist.abs() < 1.0 && muzzle_drop.is_none() {
                    muzzle_drop = Some(drop);
                }
                if (dist - 100.0).abs() < 6.0 && drop100.is_none() {
                    drop100 = Some(drop);
                    drift100 = Some(drift);
                }
            }
        }
        let muzzle_drop = muzzle_drop.expect("muzzle row");
        let (d, w) = (drop100.expect("100yd drop"), drift100.expect("100yd drift"));
        // At 90 degrees the bore is BESIDE the scope: t=0 vertical offset from the
        // LOS is sh*cos(90) = 0, not the full sight height.
        assert!(
            muzzle_drop.abs() < 0.3,
            "muzzle drop must be ~0 at 90 deg cant (bore beside scope), got {muzzle_drop}"
        );
        // Symmetry: |vertical miss| == |lateral miss| within 5%.
        let ratio = d.abs() / w.abs();
        assert!(
            (0.95..=1.05).contains(&ratio),
            "90-deg cant symmetry broken: drop {d} vs drift {w} (ratio {ratio})"
        );
    }

    // -----------------------------------------------------------------------------
    // MBA-1325: mover ring (`trajectory --target-speed`) WASM parity.
    // -----------------------------------------------------------------------------

    #[wasm_bindgen_test]
    fn target_speed_adds_ring_column_to_trajectory_table() {
        let wasm = WasmBallistics::new();
        let without = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --full --max-range 500")
            .unwrap();
        assert!(!without.contains("Ring"), "Ring column must not appear without --target-speed");

        let with_ring = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --full --max-range 500 \
                 --target-speed 5",
            )
            .unwrap();
        assert!(with_ring.contains("Ring"), "table should carry a Ring column:\n{with_ring}");
        assert!(with_ring.contains("mil"), "ring values should be labeled mil");
    }

    #[wasm_bindgen_test]
    fn target_speed_adds_mover_ring_m_to_trajectory_json() {
        let wasm = WasmBallistics::new();
        let without = wasm
            .run_command(
                "trajectory --units metric -v 823 -b 0.475 -m 10.9 -d 7.82 --full \
                 --max-range 400 -o json",
            )
            .unwrap();
        assert!(
            !without.contains("mover_ring"),
            "mover_ring fields must not appear without --target-speed"
        );

        let result = wasm
            .run_command(
                "trajectory --units metric -v 823 -b 0.475 -m 10.9 -d 7.82 --full \
                 --max-range 400 --target-speed 3 -o json",
            )
            .unwrap();
        let json: serde_json::Value = serde_json::from_str(&result).unwrap();
        let points = json["trajectory"].as_array().expect("trajectory array");
        assert!(points.len() > 5);

        // mover_ring_m is EXACTLY target_speed_mps * time at every point (metric run,
        // so --target-speed 3 is already 3 m/s, matching the native CLI convention).
        let mut checked = 0;
        for p in points {
            let time = p["time_seconds"].as_f64().expect("time_seconds");
            let ring_m = p["mover_ring_m"].as_f64().expect("mover_ring_m present");
            let expected = 3.0 * time;
            assert!(
                (ring_m - expected).abs() < 1e-9,
                "ring_m {ring_m} != target_speed_mps*time {expected} at t={time}"
            );
            checked += 1;
        }
        assert!(checked > 5);

        // Muzzle point: mover_ring_m is 0 (present), mover_ring_mil is omitted (0/0
        // is degenerate) — matches the native CLI contract exactly.
        let muzzle = &points[0];
        assert_eq!(muzzle["range_meters"].as_f64().unwrap(), 0.0);
        assert_eq!(muzzle["mover_ring_m"].as_f64().unwrap(), 0.0);
        assert!(muzzle.get("mover_ring_mil").is_none());
    }

    // -----------------------------------------------------------------------------
    // MBA-1325: `lead` powder-temperature parity with `trajectory`.
    // -----------------------------------------------------------------------------

    #[wasm_bindgen_test]
    fn lead_powder_curve_matches_manually_resolved_velocity() {
        let wasm = WasmBallistics::new();
        // WASM's `lead` has no --temperature flag (it always solves at the default
        // atmosphere — a pre-existing gap unrelated to MBA-1325, which only adds the
        // POWDER flags), so drive the curve lookup directly via --powder-temp rather
        // than the ambient air temperature. 55F is the midpoint of the 40/70 curve
        // points, so this is a genuine correction, not a no-op.
        let resolved_json = wasm
            .run_command(
                "trajectory -v 2700 -b 0.19 -m 77 -d 0.224 --drag-model g7 --max-range 5 \
                 --sight-height 2.48 --powder-temp-curve 40:2620,70:2700,100:2760 \
                 --powder-temp 55 --full -o json",
            )
            .unwrap();
        let resolved_json: serde_json::Value = serde_json::from_str(&resolved_json).unwrap();
        let resolved_velocity = resolved_json["trajectory"][0]["velocity_fps"].as_f64().unwrap();
        assert!(
            (resolved_velocity - 2660.0).abs() < 0.5,
            "sanity: 55F should interpolate to ~2660 fps, got {resolved_velocity}"
        );

        let with_curve = wasm
            .run_command(
                "lead -v 2700 -b 0.19 -m 77 -d 0.224 --drag-model g7 --sight-height 2.48 \
                 --powder-temp-curve 40:2620,70:2700,100:2760 --powder-temp 55 \
                 --target-speed 5 --range 300 -o json",
            )
            .unwrap();
        let manual = wasm
            .run_command(&format!(
                "lead -v {resolved_velocity} -b 0.19 -m 77 -d 0.224 --drag-model g7 \
                 --sight-height 2.48 --target-speed 5 --range 300 -o json"
            ))
            .unwrap();

        let with_curve: serde_json::Value = serde_json::from_str(&with_curve).unwrap();
        let manual: serde_json::Value = serde_json::from_str(&manual).unwrap();
        // `resolved_velocity` was already round-tripped once through WASM's imprecise
        // fps<->m/s constants (m/s -> fps via *3.28084, not the exact reciprocal of
        // the *0.3048 used on the way in — see lead_without_powder_flags_uses_
        // velocity_verbatim above); feeding it back in via `-v` round-trips it AGAIN,
        // so the two lead solutions agree to within float noise (~1e-7), not bit-for-
        // bit. Compare numerically with a tolerance instead of assert_eq! on the raw
        // JSON values (iterations is an integer count, compared exactly).
        for field in ["lead_mil", "lead_moa", "lead", "tof_s", "intercept_range"] {
            let a = with_curve[field].as_f64().unwrap_or_else(|| panic!("{field} missing"));
            let b = manual[field].as_f64().unwrap_or_else(|| panic!("{field} missing"));
            assert!(
                (a - b).abs() < 1e-5,
                "field '{field}' differs between curve-corrected and manually-resolved lead runs: {a} vs {b}"
            );
        }
        assert_eq!(
            with_curve["iterations"], manual["iterations"],
            "iterations should match exactly (integer count)"
        );
    }

    // -----------------------------------------------------------------------------
    // MBA-1325 env-flags addendum: `lead` gains --temperature/--pressure/--humidity/
    // --altitude/--wind-speed/--wind-direction (native parity).
    // -----------------------------------------------------------------------------

    /// With NO environmental flags, `lead` output must be byte-identical to the build
    /// from before the flags existed. These literals were captured VERBATIM from the
    /// pre-env-flags build (commit c805229) via a temporary probe test run under
    /// `wasm-pack test --safari --headless` — a true golden-compare, not a
    /// same-code-both-sides tautology. JSON strings have no trailing newline
    /// (serde_json::to_string_pretty); the table format string ends with exactly one.
    #[wasm_bindgen_test]
    fn lead_no_env_flags_output_is_byte_identical_to_pre_env_build() {
        const GOLDEN_IMPERIAL_JSON: &str = r#"{
  "adjustment_unit": "mil",
  "distance_unit": "yd",
  "intercept_range": 399.9987936,
  "iterations": 0,
  "lead": 1.267051911310233,
  "lead_mil": 3.1676393318758076,
  "lead_moa": 10.890344022989026,
  "range": 400.0,
  "target_angle_deg": 90.0,
  "target_speed": 5.0,
  "target_speed_unit": "mph",
  "tof_s": 0.5183409815796777
}"#;
        const GOLDEN_METRIC_JSON: &str = r#"{
  "adjustment_unit": "mil",
  "distance_unit": "m",
  "intercept_range": 350.0,
  "iterations": 0,
  "lead": 1.4774105669002737,
  "lead_mil": 4.221173048286496,
  "lead_moa": 14.512392940008974,
  "range": 350.0,
  "target_angle_deg": 90.0,
  "target_speed": 3.0,
  "target_speed_unit": "m/s",
  "tof_s": 0.49247018896675787
}"#;
        const GOLDEN_IMPERIAL_TABLE: &str = "Moving-Target Lead\n\
===================\n\
Target: 5.0 mph at 90\u{b0} (0=away, 90=left-to-right, 180=toward, 270=right-to-left;\n\
positive lead = hold in direction of travel)\n\
\n\
Range: 400 yd\n\
Time of Flight: 0.518 s\n\
Lead: 1.27 yd (3.17 MIL / 10.89 MOA)\n\
Intercept Range: 400.0 yd\n\
Iterations: 0\n";

        let wasm = WasmBallistics::new();
        let g1 = wasm
            .run_command(
                "lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 5 --range 400 -o json",
            )
            .unwrap();
        assert_eq!(g1, GOLDEN_IMPERIAL_JSON, "imperial JSON drifted from pre-env-flags build");
        let g2 = wasm
            .run_command(
                "--units metric lead -v 823 -b 0.475 -m 10.9 -d 7.82 --target-speed 3 \
                 --range 350 -o json",
            )
            .unwrap();
        assert_eq!(g2, GOLDEN_METRIC_JSON, "metric JSON drifted from pre-env-flags build");
        let g3 = wasm
            .run_command("lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 5 --range 400")
            .unwrap();
        assert_eq!(g3, GOLDEN_IMPERIAL_TABLE, "imperial table drifted from pre-env-flags build");
    }

    /// Hot air (100F) vs the standard default (59F) must change the lead solution —
    /// lower density, less drag, shorter TOF, smaller hold — proving --temperature is
    /// actually plumbed into the solve, not just parsed. (Fails on the pre-addendum
    /// build with "Unknown flag: --temperature".)
    #[wasm_bindgen_test]
    fn lead_temperature_is_live_hot_air_shortens_lead() {
        let wasm = WasmBallistics::new();
        let base = "lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 5 --range 400 -o json";
        let std_out = wasm.run_command(base).unwrap();
        let hot_out = wasm
            .run_command(&format!("{base} --temperature 100"))
            .unwrap();
        let std_json: serde_json::Value = serde_json::from_str(&std_out).unwrap();
        let hot_json: serde_json::Value = serde_json::from_str(&hot_out).unwrap();
        let std_tof = std_json["tof_s"].as_f64().unwrap();
        let hot_tof = hot_json["tof_s"].as_f64().unwrap();
        let std_mil = std_json["lead_mil"].as_f64().unwrap();
        let hot_mil = hot_json["lead_mil"].as_f64().unwrap();
        assert!(
            hot_tof < std_tof,
            "hot air must shorten TOF: hot {hot_tof} vs standard {std_tof}"
        );
        assert!(
            hot_mil < std_mil && (std_mil - hot_mil) / std_mil > 1e-4,
            "hot air must shrink the hold by more than float noise: hot {hot_mil} vs standard {std_mil}"
        );
    }

    /// A 20 mph headwind (--wind-direction 0 = wind-FROM dead ahead) adds drag along
    /// the flight path, lengthening TOF and growing the hold vs calm air — proving the
    /// wind flags reach the WindConditions the solve actually uses.
    #[wasm_bindgen_test]
    fn lead_headwind_is_live_and_lengthens_tof() {
        let wasm = WasmBallistics::new();
        let base = "lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 5 --range 400 -o json";
        let calm_out = wasm.run_command(base).unwrap();
        let wind_out = wasm
            .run_command(&format!("{base} --wind-speed 20 --wind-direction 0"))
            .unwrap();
        let calm_json: serde_json::Value = serde_json::from_str(&calm_out).unwrap();
        let wind_json: serde_json::Value = serde_json::from_str(&wind_out).unwrap();
        let calm_tof = calm_json["tof_s"].as_f64().unwrap();
        let wind_tof = wind_json["tof_s"].as_f64().unwrap();
        assert!(
            wind_tof > calm_tof && (wind_tof - calm_tof) / calm_tof > 1e-4,
            "a 20 mph headwind must lengthen TOF beyond float noise: wind {wind_tof} vs calm {calm_tof}"
        );
    }

    #[wasm_bindgen_test]
    fn lead_without_powder_flags_uses_velocity_verbatim() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command(
                "trajectory -v 2700 -b 0.19 -m 77 -d 0.224 --drag-model g7 --max-range 5 \
                 --sight-height 2.48 --temperature 100 --full -o json",
            )
            .unwrap();
        let json: serde_json::Value = serde_json::from_str(&result).unwrap();
        let velocity = json["trajectory"][0]["velocity_fps"].as_f64().unwrap();
        // Tolerance, not exact equality: the WASM imperial round-trip (m/s -> fps via
        // a fixed *3.28084, not the exact reciprocal of the fps -> m/s *0.3048 used on
        // the way in) carries an inherent ~3e-8 relative bias (~0.0001 fps at 2700),
        // pre-existing and unrelated to MBA-1325. 1e-3 is >100x looser than that noise
        // floor while still being >1000x tighter than any real powder correction
        // (which is many fps), so this still catches a genuine verbatim-velocity
        // regression.
        assert!(
            (velocity - 2700.0).abs() < 1e-3,
            "no powder flags should keep -v 2700 verbatim regardless of --temperature, got {velocity}"
        );
    }

    // -----------------------------------------------------------------------------
    // MBA-1325 review fix pass: M1 (target-speed bound), M2 (curve parser parity),
    // M3 (ring table --adjustment-unit).
    // -----------------------------------------------------------------------------

    /// M1: --target-speed enforces the same [0, 300] bound the native flags carry via
    /// f64_range(0.0, 300.0), on BOTH WASM commands — rejected with a clear error, not
    /// silently clamped or let through.
    #[wasm_bindgen_test]
    fn target_speed_out_of_range_is_rejected_on_both_commands() {
        let wasm = WasmBallistics::new();
        let err = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 301")
            .unwrap_err();
        assert!(
            err.as_string().unwrap_or_default().contains("--target-speed must be between 0 and 300"),
            "trajectory should reject 301: {err:?}"
        );
        let err = wasm
            .run_command(
                "lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 1000000000 --range 400",
            )
            .unwrap_err();
        assert!(
            err.as_string().unwrap_or_default().contains("--target-speed must be between 0 and 300"),
            "lead should reject 1e9 mph: {err:?}"
        );
        // Boundary values stay legal (300 is lead's documented cap).
        assert!(wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 300 --max-range 200")
            .is_ok());
    }

    /// M2: duplicate temperatures in --powder-temp-curve are rejected on both WASM
    /// commands with the native parser's message (they now share one parser, so this
    /// also guards against the two copies drifting again).
    #[wasm_bindgen_test]
    fn powder_temp_curve_duplicate_temperatures_rejected() {
        let wasm = WasmBallistics::new();
        for cmd in [
            "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 5 \
             --powder-temp-curve 40:2620,40:2700,100:2760",
            "lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 5 --range 300 \
             --powder-temp-curve 40:2620,40:2700,100:2760",
        ] {
            let err = wasm.run_command(cmd).unwrap_err();
            assert!(
                err.as_string()
                    .unwrap_or_default()
                    .contains("--powder-temp-curve has duplicate temperatures"),
                "duplicate curve temps must be rejected for `{cmd}`: {err:?}"
            );
        }
    }

    /// M3: the trajectory table's Ring cell honors --adjustment-unit, and the MOA
    /// figure is exactly mil x 3.438 (the crate's MBA-724 printed-table dial constant,
    /// shared with the native table). --target-speed 300 makes ring_mil large enough
    /// (~190) that a wrong constant (exact-angle 3437.7467/1000) would shift the value
    /// by > 0.012 — past the 0.01 print quantum, so the 2dp string match below pins
    /// the constant, not just the ballpark.
    #[wasm_bindgen_test]
    fn ring_table_honors_moa_adjustment_unit() {
        let wasm = WasmBallistics::new();
        let base = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 400 \
                    --target-speed 300 --full";

        // Full-precision mil from JSON (adjustment-unit-invariant by contract).
        let json_out = wasm.run_command(&format!("{base} -o json")).unwrap();
        let json: serde_json::Value = serde_json::from_str(&json_out).unwrap();
        let mil = json["trajectory"]
            .as_array()
            .expect("points")
            .last()
            .expect("last point")["mover_ring_mil"]
            .as_f64()
            .expect("mover_ring_mil");
        assert!(mil > 160.0, "sensitivity precondition: ring_mil > 160, got {mil}");

        // Last ring cell of a table, as the printed numeric string before " <unit>".
        let last_ring_cell = |table: &str, unit: &str| -> String {
            let suffix = format!(" {unit}");
            table
                .lines()
                .filter_map(|line| {
                    let cell = line.rsplit('|').next()?.trim();
                    cell.strip_suffix(suffix.as_str()).map(str::to_string)
                })
                .last()
                .unwrap_or_else(|| panic!("no ring cell ending in '{suffix}' found:\n{table}"))
        };

        let moa_table = wasm
            .run_command(&format!("{base} --adjustment-unit moa"))
            .unwrap();
        assert!(moa_table.contains(" moa"), "moa run must label ring cells moa");
        assert!(!moa_table.contains(" mil"), "moa run must not label ring cells mil");
        assert_eq!(
            last_ring_cell(&moa_table, "moa"),
            format!("{:.2}", mil * 3.438),
            "Ring moa cell must be exactly mil x 3.438, 2dp-rounded"
        );

        let mil_table = wasm.run_command(base).unwrap();
        assert_eq!(
            last_ring_cell(&mil_table, "mil"),
            format!("{:.2}", mil),
            "default (mil) Ring cell must match JSON mover_ring_mil, 2dp-rounded"
        );
    }

    // -----------------------------------------------------------------------------
    // MBA-1328: custom Mach:Cd drag table, mirroring the BC5D bytes-loader pattern
    // (load_bc5d_table / has_bc5d_table above). loadDragTable(bytes) parses the SAME
    // CSV format native --drag-table accepts (DragTable::from_csv_str) and, once
    // loaded, is applied automatically to every trajectory/zero/lead/monte-carlo run
    // — no --use-* gate flag needed (unlike BC5D, which needs --use-bc-segments).
    // -----------------------------------------------------------------------------

    /// A deliberately non-G1/non-G7 curve: flat Cd 0.5 across the whole Mach range.
    /// 6 points, comfortably above try_new's 2-point minimum, spanning Mach 0-3.
    const FLAT_CD_CSV: &str = "mach,cd\n0.0,0.5\n0.5,0.5\n1.0,0.5\n1.5,0.5\n2.0,0.5\n3.0,0.5\n";

    #[wasm_bindgen_test]
    fn drag_table_not_loaded_by_default() {
        let wasm = WasmBallistics::new();
        assert!(!wasm.has_drag_table());
    }

    #[wasm_bindgen_test]
    fn load_drag_table_reports_point_count_and_mach_range() {
        let wasm = WasmBallistics::new();
        let summary = wasm.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();
        assert!(wasm.has_drag_table());
        assert!(summary.contains("6 points"), "got: {summary}");
        assert!(summary.contains("0.000-3.000"), "got: {summary}");
    }

    #[wasm_bindgen_test]
    fn load_drag_table_replaces_previous_table() {
        let wasm = WasmBallistics::new();
        wasm.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();
        assert!(wasm.has_drag_table());
        // A second, different (2-point) table must fully replace the first, not merge
        // with or append to it.
        let summary = wasm
            .load_drag_table(b"mach,cd\n0.5,0.30\n2.5,0.25\n")
            .unwrap();
        assert!(wasm.has_drag_table());
        assert!(summary.contains("2 points"), "got: {summary}");
    }

    /// Live-change: a table loaded via loadDragTable() must actually change the
    /// physics of a trajectory run, not just be stored inertly. Compares the flat,
    /// deliberately non-G7-shaped Cd curve above against the unloaded (G1-model) run
    /// of otherwise-identical args.
    #[wasm_bindgen_test]
    fn drag_table_live_change_alters_trajectory_output() {
        let command = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 100 -o json";

        let unloaded = WasmBallistics::new();
        let baseline = unloaded.run_command(command).unwrap();

        let loaded = WasmBallistics::new();
        loaded.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();
        let with_table = loaded.run_command(command).unwrap();

        assert_ne!(
            baseline, with_table,
            "loading a drag table must change trajectory output for identical args"
        );
    }

    /// Live-change also applies to `zero`, `lead`, and `monte-carlo` (deliberately
    /// broader than the BC5D auto-apply, which only wires into `trajectory` — see
    /// handle_trajectory_command's BC5D block vs the unconditional drag_table checks
    /// added to all four handlers).
    #[wasm_bindgen_test]
    fn drag_table_live_change_alters_zero_lead_and_monte_carlo() {
        let zero_cmd = "zero -v 2700 -b 0.475 -m 168 -d 0.308 --target-distance 300";
        let lead_cmd =
            "lead -v 2700 -b 0.475 -m 168 -d 0.308 --target-speed 10 --range 300 -o json";
        let mc_cmd = "monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 -n 50";

        for command in [zero_cmd, lead_cmd, mc_cmd] {
            let unloaded = WasmBallistics::new();
            let baseline = unloaded.run_command(command).unwrap();

            let loaded = WasmBallistics::new();
            loaded.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();
            let with_table = loaded.run_command(command).unwrap();

            assert_ne!(
                baseline, with_table,
                "loading a drag table must change `{command}` output"
            );
        }
    }

    /// Clearing reverts the physics: after load → clear, a run must be byte-identical to a
    /// never-loaded instance's run of the same args — not merely "different from the CDM
    /// run". This is exactly the property Bero needs for a G7-vs-CDM comparison on ONE
    /// instance (load → solve CDM → clear → solve G7), the inverse of the live-change test
    /// above.
    #[wasm_bindgen_test]
    fn clear_drag_table_reverts_trajectory_output_to_g7_baseline() {
        let command = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 100 -o json";

        let baseline = WasmBallistics::new().run_command(command).unwrap();

        let engine = WasmBallistics::new();
        engine.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();
        let with_table = engine.run_command(command).unwrap();
        assert_ne!(
            baseline, with_table,
            "sanity: the loaded CDM run must differ from the G7 baseline"
        );

        assert!(
            engine.clear_drag_table(),
            "clear must report it unloaded a table"
        );
        assert!(!engine.has_drag_table());
        let after_clear = engine.run_command(command).unwrap();
        assert_eq!(
            baseline, after_clear,
            "after clearDragTable the solve must revert byte-for-byte to the G7-BC baseline"
        );
    }

    /// Return value + idempotence: clear reports whether it actually unloaded a table, and
    /// a second clear (or a clear on a fresh instance) is a harmless `false` no-op.
    #[wasm_bindgen_test]
    fn clear_drag_table_reports_prior_state_and_is_idempotent() {
        let engine = WasmBallistics::new();
        assert!(
            !engine.clear_drag_table(),
            "no table loaded → clear returns false"
        );

        engine.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();
        assert!(engine.has_drag_table());
        assert!(
            engine.clear_drag_table(),
            "a loaded table → clear returns true"
        );
        assert!(
            !engine.has_drag_table(),
            "after clear, no table remains loaded"
        );
        assert!(
            !engine.clear_drag_table(),
            "a second clear is an idempotent false no-op"
        );
    }

    /// Golden-unloaded: with NO table loaded, trajectory output for a fixed command
    /// must remain byte-identical to a literal captured from the pre-MBA-1328 build
    /// (harvested via a temporary probe test run under `wasm-pack test --node`,
    /// following the MBA-1325 golden-compare methodology — see
    /// `lead_no_env_flags_output_is_byte_identical_to_pre_env_build` above). Guards
    /// against the new drag_table plumbing silently perturbing the no-table path.
    #[wasm_bindgen_test]
    fn drag_table_unloaded_output_is_byte_identical_to_pre_change_build() {
        const GOLDEN: &str = "Trajectory Calculation Results\n\
==============================\n\
\n\
Range | Drop | Drift | Velocity | Energy | Time\n\
------|------|-------|----------|--------|------\n\
0 yd   | 2.0 in | 0.0 in  | 2700 fps   | 2719 ft-lb | 0.000 s\n\
55 yd  | 2.7 in | 0.0 in  | 2595 fps   | 2512 ft-lb | 0.062 s\n\
100 yd | 4.5 in | 0.0 in  | 2510 fps   | 2350 ft-lb | 0.115 s\n\
\n\
Max Range: 100 yards\n\
Max Height: 60.0 inches\n\
Time of Flight: 0.12 seconds\n\
Impact Velocity: 2510 fps\n";

        let wasm = WasmBallistics::new();
        assert!(!wasm.has_drag_table());
        let out = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 100")
            .unwrap();
        assert_eq!(
            out, GOLDEN,
            "unloaded trajectory output drifted from the pre-MBA-1328 build"
        );
    }

    #[wasm_bindgen_test]
    fn load_drag_table_rejects_malformed_csv_with_parser_message() {
        let wasm = WasmBallistics::new();
        let err = wasm
            .load_drag_table(b"0.5,0.23\n1.0,notanumber\n")
            .unwrap_err();
        let msg = err.as_string().unwrap_or_default();
        assert!(
            msg.contains("line 2"),
            "expected the DragTable::from_csv_str parser's line-numbered message, got: {msg}"
        );
        assert!(
            !wasm.has_drag_table(),
            "a failed load must not leave a table installed"
        );
    }

    #[wasm_bindgen_test]
    fn load_drag_table_rejects_invalid_utf8_cleanly() {
        let wasm = WasmBallistics::new();
        // 0xFF is never a valid UTF-8 lead byte.
        let err = wasm
            .load_drag_table(&[0xFF, 0xFE, 0x00, 0x00])
            .unwrap_err();
        let msg = err.as_string().unwrap_or_default();
        assert!(
            msg.contains("UTF-8"),
            "expected a clean UTF-8 decode error, got: {msg}"
        );
        assert!(!wasm.has_drag_table());
    }

    /// Deliverable #4: the WASM surface never hard-requires `-b`/`--bc` (every handler
    /// defaults it, e.g. `let mut bc = 0.475;` in handle_trajectory_command) — unlike native
    /// `zero`/`monte-carlo`, which clap-require it even though it's ignored once a table is
    /// active (CLI_USAGE.md). Confirms the solve path actually succeeds end-to-end with a
    /// loaded table and the default bc, omitting `-b` entirely, across all four commands.
    #[wasm_bindgen_test]
    fn drag_table_solve_succeeds_without_explicit_bc_flag() {
        let wasm = WasmBallistics::new();
        wasm.load_drag_table(FLAT_CD_CSV.as_bytes()).unwrap();

        let trajectory = wasm
            .run_command("trajectory -v 2700 -m 168 -d 0.308 --max-range 100")
            .unwrap();
        assert!(trajectory.contains("Trajectory Calculation Results"));

        let zero = wasm
            .run_command("zero -v 2700 -m 168 -d 0.308 --target-distance 300")
            .unwrap();
        assert!(zero.contains("Zero Calculation Results"));

        let lead = wasm
            .run_command("lead -v 2700 -m 168 -d 0.308 --target-speed 10 --range 300")
            .unwrap();
        assert!(lead.contains("Moving-Target Lead"));

        let monte_carlo = wasm
            .run_command("monte-carlo -v 2700 -m 168 -d 0.308 -n 50")
            .unwrap();
        assert!(monte_carlo.contains("Monte Carlo Simulation Results"));
    }

    /// MBA-1294(a): the auto-zero "Rifle zeroed at ..." banner is a table-only human header;
    /// JSON output must stay pure machine output (a text prefix makes it unparseable).
    #[wasm_bindgen_test]
    fn test_json_output_pure_with_auto_zero() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --auto-zero 200 -o json")
            .unwrap();
        assert!(
            !result.contains("Rifle zeroed at"),
            "auto-zero banner leaked into JSON output: {result}"
        );
        // The prefix, if present, would break this parse.
        let _: serde_json::Value =
            serde_json::from_str(&result).expect("auto-zero JSON output must be pure JSON");
    }

    /// MBA-1294(a): same purity guarantee for CSV output.
    #[wasm_bindgen_test]
    fn test_csv_output_pure_with_auto_zero() {
        let wasm = WasmBallistics::new();
        let result = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --auto-zero 200 -o csv")
            .unwrap();
        assert!(
            !result.contains("Rifle zeroed at"),
            "auto-zero banner leaked into CSV output: {result}"
        );
    }

    /// MBA-1294(c): --print-bc-segments appends the BC ladder in TABLE view only; it must
    /// never contaminate a JSON payload.
    #[wasm_bindgen_test]
    fn test_print_bc_segments_table_only() {
        let wasm = WasmBallistics::new();
        let table = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --print-bc-segments")
            .unwrap();
        // A plain scalar BC yields the "none active" note; either way the report block is present.
        assert!(
            table.contains("BC segments") || table.contains("BC Segments"),
            "bc-segments report missing from table view: {table}"
        );
        let json = wasm
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --print-bc-segments -o json")
            .unwrap();
        assert!(
            !json.contains("BC Segments (active)") && !json.contains("No BC segments active"),
            "bc-segments report leaked into JSON output: {json}"
        );
        let _: serde_json::Value =
            serde_json::from_str(&json).expect("--print-bc-segments must not break JSON");
    }

    /// MBA-1294(c): the `lead` command rejects an unknown -o value instead of silently
    /// defaulting to a format; a valid value is accepted.
    #[wasm_bindgen_test]
    fn test_lead_invalid_output_rejected() {
        let wasm = WasmBallistics::new();
        assert!(
            wasm.run_command("lead -v 2700 -m 168 -d 0.308 --target-speed 10 --range 300 -o json")
                .is_ok(),
            "lead -o json should be accepted"
        );
        let err = wasm
            .run_command("lead -v 2700 -m 168 -d 0.308 --target-speed 10 --range 300 -o bogus")
            .expect_err("lead must reject an invalid -o value");
        assert!(
            format!("{err:?}").contains("Invalid --output"),
            "unexpected lead -o error: {err:?}"
        );
    }

    /// MBA-1339: the WASM CLI accepts --bore-height as an alias for --muzzle-height on the
    /// same inches/mm units (unifying with the native CLI), so the two must produce identical
    /// trajectories.
    #[wasm_bindgen_test]
    fn test_bore_height_aliases_muzzle_height() {
        let wasm = WasmBallistics::new();
        let via_muzzle = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --muzzle-height 40 --max-range 300",
            )
            .unwrap();
        let via_bore = wasm
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --bore-height 40 --max-range 300",
            )
            .unwrap();
        assert_eq!(
            via_muzzle, via_bore,
            "--bore-height must alias --muzzle-height in the WASM CLI"
        );
    }

    // --- MBA-737: `powder` command (parity with the native CLI handler) ---

    #[wasm_bindgen_test]
    fn test_trajectory_plot_appends_charts() {
        // MBA-1337 p3: --plot appends both charts after the table (table-only).
        let out = WasmBallistics::new()
            .run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 300 --plot")
            .unwrap();
        assert!(out.contains("Drop vs Range:"), "{}", out);
        assert!(out.contains("Lateral Drift vs Range:"), "{}", out);
        let ascii = WasmBallistics::new()
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 300 --plot ascii",
            )
            .unwrap();
        assert!(ascii.contains("Drop vs Range:"), "{}", ascii);
        // JSON stays pure machine output even with --plot.
        let json = WasmBallistics::new()
            .run_command(
                "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 300 --plot -o json",
            )
            .unwrap();
        assert!(!json.contains("Drop vs Range:"), "{}", json);
    }

    #[wasm_bindgen_test]
    fn test_powder_linear_model_mba1296_pin() {
        // The MBA-1296 repro: 1.0 fps/degF, 40F day, 70F reference -> 2770, never 8400.
        let out = WasmBallistics::new()
            .run_command("powder -v 2800 --temperature 40")
            .unwrap();
        assert!(out.contains("2770.0 fps"), "{}", out);
        assert!(out.contains("(-30.0)"), "{}", out);
    }

    #[wasm_bindgen_test]
    fn test_powder_curve_interpolates_at_powder_temp() {
        let out = WasmBallistics::new()
            .run_command("powder --powder-temp-curve 40:2620,70:2700,100:2760 --powder-temp 55")
            .unwrap();
        assert!(out.contains("2660.0 fps"), "{}", out);
        assert!(out.contains("measured curve, 3 points"), "{}", out);
    }

    #[wasm_bindgen_test]
    fn test_powder_sweep_csv_includes_end_row() {
        let out = WasmBallistics::new()
            .run_command("powder -v 2700 --sweep 20:110:30 -o csv")
            .unwrap();
        assert!(out.contains("2650.0,-50.0"), "{}", out);
        // The END row must not be dropped to float drift in the sweep expansion.
        assert!(out.contains("2740.0,40.0"), "{}", out);
    }

    #[wasm_bindgen_test]
    fn test_powder_json_resolved_object() {
        let out = WasmBallistics::new()
            .run_command("powder -v 2800 --temperature 40 -o json")
            .unwrap();
        assert!(out.contains("\"velocity\": 2770.0"), "{}", out);
        assert!(out.contains("\"shift\": -30.0"), "{}", out);
    }

    #[wasm_bindgen_test]
    fn test_powder_metric_resolves_in_celsius() {
        // 0.54864 m/s/degC x (5 - 21) degC = -8.778 m/s -> 814.2
        let out = WasmBallistics::new()
            .run_command("--units metric powder -v 823 --temperature 5 --powder-temp 21")
            .unwrap();
        assert!(out.contains("814.2 m/s"), "{}", out);
    }

    #[wasm_bindgen_test]
    fn test_powder_missing_inputs_error() {
        assert!(WasmBallistics::new()
            .run_command("powder --temperature 40")
            .is_err());
        assert!(WasmBallistics::new().run_command("powder -v 2800").is_err());
    }

    #[wasm_bindgen_test]
    fn test_powder_energy_with_mass() {
        // 168 gr at 2700 fps -> ~2719 ft-lb
        let out = WasmBallistics::new()
            .run_command("powder -v 2700 --temperature 70 -m 168")
            .unwrap();
        assert!(out.contains("2719 ft\u{b7}lb"), "{}", out);
    }
}
