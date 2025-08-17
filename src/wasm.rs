// WASM bindings for the ballistics engine
use wasm_bindgen::prelude::*;

use crate::cli_api::{
    TrajectorySolver, 
    BallisticInputs as InternalBallisticInputs, 
    WindConditions, AtmosphericConditions,
};
use crate::drag_model::DragModel;

#[wasm_bindgen]
pub struct WasmBallistics;

#[wasm_bindgen]
impl WasmBallistics {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        WasmBallistics
    }

    /// Run a command and return the output
    #[wasm_bindgen(js_name = runCommand)]
    pub fn run_command(&self, command: &str) -> Result<String, JsValue> {
        // Parse simple trajectory command
        let args: Vec<&str> = command.split_whitespace().collect();
        
        if args.is_empty() {
            return Ok(self.show_help());
        }

        // Skip "ballistics" if present
        let args = if !args.is_empty() && (args[0] == "ballistics" || args[0] == "./ballistics") {
            &args[1..]
        } else {
            &args[..]
        };

        if args.is_empty() || args[0] == "help" || args[0] == "--help" || args[0] == "-h" {
            return Ok(self.show_help());
        }

        if args[0] != "trajectory" {
            return Ok(format!("Error: Unknown command '{}'\n\n{}", args[0], self.show_help()));
        }

        // Parse trajectory arguments
        let mut velocity = 2700.0; // fps
        let mut bc = 0.475;
        let mut mass = 168.0; // grains
        let mut diameter = 0.308; // inches
        let mut drag_model = "G1";
        let mut max_range = 1000.0; // yards
        let mut wind_speed = 0.0; // mph
        let mut wind_direction = 90.0; // degrees

        let mut i = 1;
        while i < args.len() {
            match args[i] {
                "-v" | "--velocity" => {
                    if i + 1 < args.len() {
                        velocity = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid velocity"))?;
                        i += 1;
                    }
                }
                "-b" | "--bc" => {
                    if i + 1 < args.len() {
                        bc = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid BC"))?;
                        i += 1;
                    }
                }
                "-m" | "--mass" => {
                    if i + 1 < args.len() {
                        mass = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid mass"))?;
                        i += 1;
                    }
                }
                "-d" | "--diameter" => {
                    if i + 1 < args.len() {
                        diameter = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid diameter"))?;
                        i += 1;
                    }
                }
                "--drag-model" => {
                    if i + 1 < args.len() {
                        drag_model = args[i + 1];
                        i += 1;
                    }
                }
                "--max-range" => {
                    if i + 1 < args.len() {
                        max_range = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid max range"))?;
                        i += 1;
                    }
                }
                "--wind-speed" => {
                    if i + 1 < args.len() {
                        wind_speed = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid wind speed"))?;
                        i += 1;
                    }
                }
                "--wind-direction" => {
                    if i + 1 < args.len() {
                        wind_direction = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid wind direction"))?;
                        i += 1;
                    }
                }
                _ => {}
            }
            i += 1;
        }

        // Build inputs
        let mut inputs = InternalBallisticInputs::default();
        
        // Convert from imperial to metric
        inputs.muzzle_velocity = velocity * 0.3048; // fps to m/s
        inputs.ballistic_coefficient = bc;
        inputs.mass = mass * 0.00006479891; // grains to kg
        inputs.diameter = diameter * 0.0254; // inches to meters
        inputs.drag_model = DragModel::from_str(drag_model)
            .ok_or_else(|| JsValue::from_str("Invalid drag model"))?;
        
        // Set sight height to 2 inches (typical scope height)
        inputs.sight_height = 2.0 * 0.0254; // inches to meters
        
        // Set wind
        let mut wind = WindConditions::default();
        wind.speed = wind_speed * 0.44704; // mph to m/s
        wind.direction = wind_direction;
        
        // Default atmosphere
        let atmosphere = AtmosphericConditions::default();
        
        // Create solver and calculate
        let mut solver = TrajectorySolver::new(inputs, wind, atmosphere);
        solver.set_max_range(max_range * 0.9144); // yards to meters
        
        match solver.solve() {
            Ok(result) => {
                Ok(self.format_trajectory_output(&result))
            }
            Err(e) => Ok(format!("Error: {}", e))
        }
    }

    fn show_help(&self) -> String {
        r#"Ballistics Engine - WebAssembly Version

Usage: ballistics trajectory [OPTIONS]

Options:
  -v, --velocity <VEL>         Muzzle velocity (fps) [default: 2700]
  -b, --bc <BC>                Ballistic coefficient [default: 0.475]
  -m, --mass <MASS>            Mass (grains) [default: 168]
  -d, --diameter <DIAMETER>    Diameter (inches) [default: 0.308]
  --drag-model <MODEL>         Drag model (G1/G7) [default: G1]
  --max-range <RANGE>          Maximum range (yards) [default: 1000]
  --wind-speed <SPEED>         Wind speed (mph) [default: 0]
  --wind-direction <DIR>       Wind direction (degrees, 90=right) [default: 90]

Example:
  ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308 --max-range 1000

Try it now with default values:
  ballistics trajectory"#.to_string()
    }

    fn format_trajectory_output(&self, result: &crate::cli_api::TrajectoryResult) -> String {
        let mut output = String::new();
        output.push_str("Trajectory Calculation Results\n");
        output.push_str("==============================\n\n");
        output.push_str("Range(yd) | Drop(in) | Drift(in) | Velocity(fps) | Energy(ft-lb) | Time(s)\n");
        output.push_str("----------|----------|-----------|---------------|---------------|--------\n");

        // Sample every 100 yards, or more frequently if max range is small
        let max_range_yards = result.max_range * 1.09361;
        let sample_interval = if max_range_yards < 500.0 { 50.0 } else { 100.0 };
        let mut current_range = 0.0;
        
        // Get initial height for reference (sight height)
        let sight_height_inches = if !result.points.is_empty() {
            result.points[0].position.y * 39.3701
        } else {
            2.0
        };
        
        for point in &result.points {
            let range_yards = point.position.z * 1.09361; // m to yards
            
            if range_yards >= current_range || (result.points.last() == Some(point)) {
                // Calculate drop relative to line of sight (not absolute position)
                let drop_inches = (sight_height_inches - point.position.y * 39.3701);
                let drift_inches = point.position.x * 39.3701; // m to inches
                let velocity_fps = point.velocity_magnitude * 3.28084; // m/s to fps
                
                // Calculate energy from velocity and mass
                // KE = 0.5 * mass * velocity^2
                // Convert to ft-lbs: divide by 450240
                let mass_grains = 168.0; // Default mass, should get from result if available
                let energy_ftlb = 0.5 * mass_grains * velocity_fps * velocity_fps / 450240.0;
                
                output.push_str(&format!(
                    "{:9.0} | {:8.1} | {:9.1} | {:13.0} | {:13.0} | {:7.3}\n",
                    range_yards, drop_inches, drift_inches, velocity_fps, energy_ftlb, point.time
                ));
                
                current_range += sample_interval;
            }
        }
        
        output.push_str(&format!("\nMax Range: {:.0} yards\n", result.max_range * 1.09361));
        output.push_str(&format!("Max Height: {:.1} inches\n", result.max_height * 39.3701));
        output.push_str(&format!("Time of Flight: {:.2} seconds\n", result.time_of_flight));
        output.push_str(&format!("Impact Velocity: {:.0} fps\n", result.impact_velocity * 3.28084));

        output
    }
}