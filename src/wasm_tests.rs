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
            77.0 * 0.00006479891,
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
        let command = "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
                       --auto-zero 300 --zero-velocity 2400 --zero-temperature 20 \
                       --max-range 100 -o json";
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
        let command = "trajectory -v 2650 -b 0.19 -m 77 -d 0.224 --drag-model g7 \
                       --temperature 32 --powder-temp 68 \
                       --powder-temp-curve 32:2650,77:2720 \
                       --auto-zero 100 --max-range 100 -o json";
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
}
