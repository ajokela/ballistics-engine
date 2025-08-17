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
        let result = wasm.run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308").unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
        assert!(result.contains("Range"));
        assert!(result.contains("Drop"));
        assert!(result.contains("Velocity"));
    }

    #[wasm_bindgen_test]
    fn test_trajectory_with_auto_zero() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --auto-zero 200").unwrap();
        assert!(result.contains("Rifle zeroed at 200 yards"));
        assert!(result.contains("MOA"));
        assert!(result.contains("mrad"));
    }

    #[wasm_bindgen_test]
    fn test_metric_units() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command("--units metric trajectory -v 823 -b 0.475 -m 10.9 -d 7.82").unwrap();
        assert!(result.contains("m/s"));
        assert!(result.contains("meters"));
        assert!(!result.contains("yards"));
        assert!(!result.contains("fps"));
    }

    #[wasm_bindgen_test]
    fn test_json_output_format() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 -o json").unwrap();
        assert!(result.contains("\"trajectory\""));
        assert!(result.contains("\"summary\""));
        assert!(result.contains("range_yards"));
        assert!(result.contains("drop_inches"));
    }

    #[wasm_bindgen_test]
    fn test_csv_output_format() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 -o csv").unwrap();
        assert!(result.contains("Range(yards),Drop(inches)"));
        assert!(result.contains(","));
        assert!(!result.contains("Trajectory Calculation Results")); // CSV shouldn't have the header
    }

    #[wasm_bindgen_test]
    fn test_zero_command() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command("zero -v 2700 -b 0.475 -m 168 -d 0.308 --target-distance 300").unwrap();
        assert!(result.contains("Zero Calculation Results"));
        assert!(result.contains("Target Distance"));
        assert!(result.contains("MOA Adjustment"));
        assert!(result.contains("Mrad Adjustment"));
    }

    #[wasm_bindgen_test]
    fn test_monte_carlo_command() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command("monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 -n 100").unwrap();
        assert!(result.contains("Monte Carlo Simulation Results"));
        assert!(result.contains("Simulations Run: 100"));
        assert!(result.contains("Range Statistics"));
        assert!(result.contains("Impact Velocity Statistics"));
    }

    #[wasm_bindgen_test]
    fn test_estimate_bc_command() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command("estimate-bc -v 2700 -m 168 -d 0.308 --data \"100,2.5;200,10.2;300,23.5\"").unwrap();
        assert!(result.contains("BC Estimation Results"));
        assert!(result.contains("Estimated BC"));
        assert!(result.contains("Based on 3 data points"));
    }

    #[wasm_bindgen_test]
    fn test_advanced_physics_flags() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command(
            "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
             --enable-magnus --enable-coriolis --enable-spin-drift \
             --twist-rate 10 --latitude 45"
        ).unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
        // The trajectory should complete without errors when physics flags are enabled
    }

    #[wasm_bindgen_test]
    fn test_environmental_conditions() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command(
            "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
             --temperature 32 --pressure 25.0 --humidity 80 --altitude 5000 \
             --wind-speed 10 --wind-direction 90"
        ).unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn test_full_trajectory_output() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --full --max-range 500").unwrap();
        // With --full flag, we should see more data points
        let lines: Vec<&str> = result.lines().collect();
        let data_lines = lines.iter().filter(|l| l.contains("yd")).count();
        assert!(data_lines > 5); // Should have many data points with --full
    }

    #[wasm_bindgen_test]
    fn test_powder_temperature_sensitivity() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command(
            "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
             --use-powder-sensitivity --powder-temp 90 --powder-temp-sensitivity 1.5"
        ).unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn test_shooting_angle() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command(
            "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --shooting-angle 15"
        ).unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn test_drag_models() {
        let wasm = WasmBallistics::new();
        
        // Test G1 model
        let g1_result = wasm.run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --drag-model G1").unwrap();
        assert!(g1_result.contains("Trajectory Calculation Results"));
        
        // Test G7 model
        let g7_result = wasm.run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --drag-model G7").unwrap();
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
        let result = wasm.run_command(
            "--units metric monte-carlo -v 823 -b 0.475 -m 10.9 -d 7.82 -n 50"
        ).unwrap();
        assert!(result.contains("m/s"));
        assert!(result.contains("meters"));
        assert!(!result.contains("yards"));
        assert!(!result.contains("fps"));
    }

    #[wasm_bindgen_test]
    fn test_metric_zero_calculation() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command(
            "--units metric zero -v 823 -b 0.475 -m 10.9 -d 7.82 --target-distance 200"
        ).unwrap();
        assert!(result.contains("200 meters"));
        assert!(!result.contains("yards"));
    }

    #[wasm_bindgen_test]
    fn test_bc_segments() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command(
            "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --use-bc-segments"
        ).unwrap();
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
        let rk4_result = wasm.run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308").unwrap();
        assert!(rk4_result.contains("Trajectory Calculation Results"));
        
        // Test Euler
        let euler_result = wasm.run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --use-euler").unwrap();
        assert!(euler_result.contains("Trajectory Calculation Results"));
    }

    #[wasm_bindgen_test]
    fn test_estimate_bc_no_data() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command("estimate-bc -v 2700 -m 168 -d 0.308").unwrap();
        assert!(result.contains("Error: No trajectory data provided"));
        assert!(result.contains("Example"));
    }

    #[wasm_bindgen_test]
    fn test_zero_with_sight_height() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command(
            "zero -v 2700 -b 0.475 -m 168 -d 0.308 --target-distance 200 --sight-height 2.5"
        ).unwrap();
        assert!(result.contains("Sight Height: 2.5 inches"));
    }

    #[wasm_bindgen_test]
    fn test_monte_carlo_with_variations() {
        let wasm = WasmBallistics::new();
        let result = wasm.run_command(
            "monte-carlo -v 2700 -b 0.475 -m 168 -d 0.308 -n 50 \
             --velocity-std 15 --angle-std 0.2 --bc-std 0.02 \
             --wind-speed-std 3 --wind-dir-std 10"
        ).unwrap();
        assert!(result.contains("Monte Carlo Simulation Results"));
        assert!(result.contains("Std Dev"));
    }

    #[wasm_bindgen_test]
    fn test_command_parsing() {
        let wasm = WasmBallistics::new();
        
        // Test that ballistics prefix is handled
        let result1 = wasm.run_command("ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308").unwrap();
        assert!(result1.contains("Trajectory Calculation Results"));
        
        // Test without prefix
        let result2 = wasm.run_command("trajectory -v 2700 -b 0.475 -m 168 -d 0.308").unwrap();
        assert!(result2.contains("Trajectory Calculation Results"));
        
        // Test with ./ prefix
        let result3 = wasm.run_command("./ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308").unwrap();
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
    fn test_all_physics_flags_combined() {
        let wasm = WasmBallistics::new();
        
        // Test all physics flags together
        let result = wasm.run_command(
            "trajectory -v 2700 -b 0.475 -m 168 -d 0.308 \
             --enable-magnus --enable-coriolis --enable-spin-drift \
             --enable-wind-shear --enable-pitch-damping --enable-precession \
             --sample-trajectory --use-bc-segments --use-powder-sensitivity \
             --twist-rate 10 --latitude 45"
        ).unwrap();
        assert!(result.contains("Trajectory Calculation Results"));
    }
}