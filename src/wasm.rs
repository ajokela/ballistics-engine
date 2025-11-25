// WASM bindings for the ballistics engine with full CLI feature parity
use wasm_bindgen::prelude::*;
use serde_json;

use crate::cli_api::{
    TrajectorySolver, 
    BallisticInputs as InternalBallisticInputs, 
    WindConditions, AtmosphericConditions,
    calculate_zero_angle_with_conditions,
    run_monte_carlo,
    MonteCarloParams,
    estimate_bc_from_trajectory,
};
use crate::drag_model::DragModel;

#[wasm_bindgen]
pub struct WasmBallistics;

// Unit system for conversions
#[derive(Debug, Clone, Copy, PartialEq)]
enum UnitSystem {
    Imperial,
    Metric,
}

impl UnitSystem {
    fn from_str(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "metric" => UnitSystem::Metric,
            _ => UnitSystem::Imperial,
        }
    }
}

// Output format
#[derive(Debug, Clone, Copy, PartialEq)]
enum OutputFormat {
    Table,
    Json,
    Csv,
}

impl OutputFormat {
    fn from_str(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "json" => OutputFormat::Json,
            "csv" => OutputFormat::Csv,
            _ => OutputFormat::Table,
        }
    }
}

#[wasm_bindgen]
impl WasmBallistics {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        WasmBallistics
    }

    /// Run a command and return the output
    #[wasm_bindgen(js_name = runCommand)]
    pub fn run_command(&self, command: &str) -> Result<String, JsValue> {
        // Handle quoted arguments properly
        let mut args: Vec<String> = Vec::new();
        let mut current_arg = String::new();
        let mut in_quotes = false;
        let mut quote_char = ' ';
        
        for c in command.chars() {
            if !in_quotes && (c == '\'' || c == '"') {
                in_quotes = true;
                quote_char = c;
            } else if in_quotes && c == quote_char {
                in_quotes = false;
                quote_char = ' ';
            } else if !in_quotes && c.is_whitespace() {
                if !current_arg.is_empty() {
                    args.push(current_arg.clone());
                    current_arg.clear();
                }
            } else {
                current_arg.push(c);
            }
        }
        
        if !current_arg.is_empty() {
            args.push(current_arg);
        }
        
        let args: Vec<&str> = args.iter().map(|s| s.as_str()).collect();
        
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

        // Parse global unit system first
        let mut units = UnitSystem::Imperial;
        for i in 0..args.len() {
            if args[i] == "--units" || args[i] == "-u" {
                if i + 1 < args.len() {
                    units = UnitSystem::from_str(args[i + 1]);
                }
                break;
            }
        }

        // Route to appropriate command handler
        match args[0] {
            "version" => Ok("Ballistics Engine v0.4.2\nWASM Build\n".to_string()),
            "trajectory" => self.handle_trajectory_command(&args[1..], units),
            "zero" => self.handle_zero_command(&args[1..], units),
            "monte-carlo" | "montecarlo" => self.handle_monte_carlo_command(&args[1..], units),
            "estimate-bc" => self.handle_estimate_bc_command(&args[1..], units),
            _ => Ok(format!("Error: Unknown command '{}'\n\n{}", args[0], self.show_help())),
        }
    }

    fn handle_trajectory_command(&self, args: &[&str], units: UnitSystem) -> Result<String, JsValue> {
        // Default values based on unit system
        let (default_velocity, default_mass, default_diameter, default_temp, default_pressure) = match units {
            UnitSystem::Imperial => (2700.0, 168.0, 0.308, 59.0, 29.92),
            UnitSystem::Metric => (823.0, 10.9, 7.82, 15.0, 1013.25),
        };

        // Initialize all parameters with defaults
        let mut velocity = default_velocity;
        let mut angle = 0.0;
        let mut bc = 0.475;
        let mut mass = default_mass;
        let mut diameter = default_diameter;
        let mut drag_model = "G1";
        let mut max_range = if units == UnitSystem::Imperial { 1000.0 } else { 914.4 };
        let mut time_step = 0.001;
        let mut wind_speed = 0.0;
        let mut wind_direction = 90.0;
        let mut temperature = default_temp;
        let mut pressure = default_pressure;
        let mut humidity = 50.0;
        let mut altitude = 0.0;
        let mut output_format = OutputFormat::Table;
        let mut full = false;
        let mut auto_zero: Option<f64> = None;
        let mut sight_height = if units == UnitSystem::Imperial { 2.0 } else { 50.0 }; // inches or mm
        let mut muzzle_height = if units == UnitSystem::Imperial { 60.0 } else { 1500.0 }; // inches or mm (standing)
        let mut target_height = 0.0; // inches or mm (ground level)
        
        // Advanced physics flags
        let mut enable_magnus = false;
        let mut enable_coriolis = false;
        let mut use_euler = false;
        let mut use_rk4_fixed = false;  // Use fixed-step RK4 instead of adaptive RK45
        let mut enable_spin_drift = false;
        let mut enable_wind_shear = false;
        let mut sample_trajectory = false;
        let mut sample_interval = 10.0;
        let mut enable_pitch_damping = false;
        let mut enable_precession = false;
        let mut use_bc_segments = false;
        let mut use_powder_sensitivity = false;
        
        // Additional parameters
        let mut twist_rate: Option<f64> = None;
        let mut twist_right = true;
        let mut latitude: Option<f64> = None;
        let mut shooting_angle = 0.0;
        let mut powder_temp_sensitivity = 1.0;
        let mut powder_temp = if units == UnitSystem::Imperial { 70.0 } else { 21.0 };

        // Parse arguments
        let mut i = 0;
        while i < args.len() {
            match args[i] {
                "-v" | "--velocity" => {
                    if i + 1 < args.len() {
                        velocity = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid velocity"))?;
                        i += 1;
                    }
                }
                "-a" | "--angle" => {
                    if i + 1 < args.len() {
                        angle = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid angle"))?;
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
                "--time-step" => {
                    if i + 1 < args.len() {
                        time_step = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid time step"))?;
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
                "--temperature" => {
                    if i + 1 < args.len() {
                        temperature = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid temperature"))?;
                        i += 1;
                    }
                }
                "--pressure" => {
                    if i + 1 < args.len() {
                        pressure = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid pressure"))?;
                        i += 1;
                    }
                }
                "--humidity" => {
                    if i + 1 < args.len() {
                        humidity = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid humidity"))?;
                        i += 1;
                    }
                }
                "--altitude" => {
                    if i + 1 < args.len() {
                        altitude = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid altitude"))?;
                        i += 1;
                    }
                }
                "-o" | "--output" => {
                    if i + 1 < args.len() {
                        output_format = OutputFormat::from_str(args[i + 1]);
                        i += 1;
                    }
                }
                "--full" => full = true,
                "--auto-zero" | "-z" => {
                    if i + 1 < args.len() {
                        auto_zero = Some(args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid zero distance"))?);
                        i += 1;
                    }
                }
                "--sight-height" => {
                    if i + 1 < args.len() {
                        sight_height = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid sight height"))?;
                        i += 1;
                    }
                }
                "--muzzle-height" => {
                    if i + 1 < args.len() {
                        muzzle_height = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid muzzle height"))?;
                        i += 1;
                    }
                }
                "--target-height" => {
                    if i + 1 < args.len() {
                        target_height = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid target height"))?;
                        i += 1;
                    }
                }
                "--enable-magnus" => enable_magnus = true,
                "--enable-coriolis" => enable_coriolis = true,
                "--use-euler" => use_euler = true,
                "--use-rk4-fixed" => use_rk4_fixed = true,
                "--enable-spin-drift" => enable_spin_drift = true,
                "--enable-wind-shear" => enable_wind_shear = true,
                "--sample-trajectory" => sample_trajectory = true,
                "--sample-interval" => {
                    if i + 1 < args.len() {
                        sample_interval = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid sample interval"))?;
                        i += 1;
                    }
                }
                "--enable-pitch-damping" => enable_pitch_damping = true,
                "--enable-precession" => enable_precession = true,
                "--use-bc-segments" => use_bc_segments = true,
                "--use-powder-sensitivity" => use_powder_sensitivity = true,
                "--twist-rate" => {
                    if i + 1 < args.len() {
                        twist_rate = Some(args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid twist rate"))?);
                        i += 1;
                    }
                }
                "--twist-right" => {
                    if i + 1 < args.len() {
                        twist_right = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid twist direction"))?;
                        i += 1;
                    }
                }
                "--latitude" => {
                    if i + 1 < args.len() {
                        latitude = Some(args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid latitude"))?);
                        i += 1;
                    }
                }
                "--shooting-angle" => {
                    if i + 1 < args.len() {
                        shooting_angle = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid shooting angle"))?;
                        i += 1;
                    }
                }
                "--powder-temp-sensitivity" => {
                    if i + 1 < args.len() {
                        powder_temp_sensitivity = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid powder temp sensitivity"))?;
                        i += 1;
                    }
                }
                "--powder-temp" => {
                    if i + 1 < args.len() {
                        powder_temp = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid powder temp"))?;
                        i += 1;
                    }
                }
                _ => {}
            }
            i += 1;
        }

        // Build inputs with unit conversions
        let mut inputs = InternalBallisticInputs::default();
        
        // Convert units to metric (internal representation)
        match units {
            UnitSystem::Imperial => {
                inputs.muzzle_velocity = velocity * 0.3048; // fps to m/s
                inputs.bullet_mass = mass * 0.00006479891; // grains to kg
                inputs.bullet_diameter = diameter * 0.0254; // inches to meters
                inputs.sight_height = sight_height * 0.0254; // inches to meters
                inputs.muzzle_height = muzzle_height * 0.0254; // inches to meters
                inputs.target_height = target_height * 0.0254; // inches to meters
            }
            UnitSystem::Metric => {
                inputs.muzzle_velocity = velocity; // already m/s
                inputs.bullet_mass = mass * 0.001; // grams to kg
                inputs.bullet_diameter = diameter * 0.001; // mm to meters
                inputs.sight_height = sight_height * 0.001; // mm to meters
                inputs.muzzle_height = muzzle_height * 0.001; // mm to meters
                inputs.target_height = target_height * 0.001; // mm to meters
            }
        }
        
        inputs.bc_value = bc;
        inputs.bc_type = DragModel::from_str(drag_model)
            .ok_or_else(|| JsValue::from_str("Invalid drag model"))?;
        inputs.muzzle_angle = angle * std::f64::consts::PI / 180.0; // degrees to radians
        inputs.shooting_angle = shooting_angle * std::f64::consts::PI / 180.0;
        
        // Set advanced physics flags
        // Magnus and Coriolis are controlled by enable_advanced_effects
        if enable_magnus || enable_coriolis {
            inputs.enable_advanced_effects = true;
        }
        // Set integration method: Euler < RK4 fixed < RK45 adaptive (default)
        if use_euler {
            inputs.use_rk4 = false;
            inputs.use_adaptive_rk45 = false;
        } else if use_rk4_fixed {
            inputs.use_rk4 = true;
            inputs.use_adaptive_rk45 = false;  // Fixed-step RK4
        } else {
            inputs.use_rk4 = true;
            inputs.use_adaptive_rk45 = true;   // Default: adaptive RK45
        }
        inputs.use_enhanced_spin_drift = enable_spin_drift;
        inputs.enable_wind_shear = enable_wind_shear;
        inputs.enable_trajectory_sampling = sample_trajectory;
        inputs.sample_interval = sample_interval;
        inputs.enable_pitch_damping = enable_pitch_damping;
        inputs.enable_precession_nutation = enable_precession;
        inputs.use_bc_segments = use_bc_segments;
        inputs.use_powder_sensitivity = use_powder_sensitivity;
        
        // Set additional parameters
        if let Some(rate) = twist_rate {
            inputs.twist_rate = rate;
        }
        inputs.is_twist_right = twist_right;
        if let Some(lat) = latitude {
            inputs.latitude = Some(lat);
        }
        inputs.powder_temp_sensitivity = powder_temp_sensitivity;
        
        // Adjust velocity for powder temperature if enabled
        if use_powder_sensitivity {
            let temp_diff = match units {
                UnitSystem::Imperial => powder_temp - 70.0,
                UnitSystem::Metric => powder_temp - 21.0,
            };
            let velocity_adjustment = temp_diff * powder_temp_sensitivity;
            inputs.muzzle_velocity += velocity_adjustment * (if units == UnitSystem::Imperial { 0.3048 } else { 1.0 });
        }
        
        // Set wind conditions
        let mut wind = WindConditions::default();
        match units {
            UnitSystem::Imperial => {
                wind.speed = wind_speed * 0.44704; // mph to m/s
            }
            UnitSystem::Metric => {
                wind.speed = wind_speed; // already m/s
            }
        }
        wind.direction = wind_direction;
        
        // Set atmospheric conditions
        let mut atmosphere = AtmosphericConditions::default();
        match units {
            UnitSystem::Imperial => {
                atmosphere.temperature = (temperature - 32.0) * 5.0 / 9.0; // F to C
                atmosphere.pressure = pressure * 33.863886666667; // inHg to hPa
                atmosphere.altitude = altitude * 0.3048; // feet to meters
            }
            UnitSystem::Metric => {
                atmosphere.temperature = temperature;
                atmosphere.pressure = pressure;
                atmosphere.altitude = altitude;
            }
        }
        atmosphere.humidity = humidity;
        
        // Handle auto-zero if specified
        let mut zero_info = String::new();
        if let Some(zero_distance) = auto_zero {
            let zero_distance_m = match units {
                UnitSystem::Imperial => zero_distance * 0.9144, // yards to meters
                UnitSystem::Metric => zero_distance,
            };
            
            match calculate_zero_angle_with_conditions(
                inputs.clone(), 
                zero_distance_m,
                inputs.target_height, // Use target height from inputs
                wind.clone(), 
                atmosphere.clone()
            ) {
                Ok(zero_angle) => {
                    inputs.muzzle_angle = zero_angle;
                    let moa_adjustment = zero_angle * 180.0 / std::f64::consts::PI * 60.0;
                    let mrad_adjustment = zero_angle * 1000.0;
                    zero_info = format!("Rifle zeroed at {} {} (Adjustment: {:.2} MOA / {:.2} mrad up)\n\n", 
                                      zero_distance, 
                                      if units == UnitSystem::Imperial { "yards" } else { "meters" },
                                      moa_adjustment,
                                      mrad_adjustment);
                }
                Err(e) => {
                    return Ok(format!("Error calculating zero: {}\n\nTry a shorter zero distance or check your ballistic parameters.", e));
                }
            }
        }
        
        // Create solver and calculate
        let mut solver = TrajectorySolver::new(inputs.clone(), wind, atmosphere);
        let max_range_m = match units {
            UnitSystem::Imperial => max_range * 0.9144, // yards to meters
            UnitSystem::Metric => max_range,
        };
        solver.set_max_range(max_range_m);
        solver.set_time_step(time_step);
        
        match solver.solve() {
            Ok(result) => {
                let output = match output_format {
                    OutputFormat::Table => self.format_trajectory_table(&result, auto_zero, units, full),
                    OutputFormat::Json => self.format_trajectory_json(&result, units),
                    OutputFormat::Csv => self.format_trajectory_csv(&result, units, full),
                };
                Ok(format!("{}{}", zero_info, output))
            }
            Err(e) => Ok(format!("Error: {}", e))
        }
    }

    fn handle_zero_command(&self, args: &[&str], units: UnitSystem) -> Result<String, JsValue> {
        // Default values
        let (default_velocity, default_mass, default_diameter) = match units {
            UnitSystem::Imperial => (2700.0, 168.0, 0.308),
            UnitSystem::Metric => (823.0, 10.9, 7.82),
        };

        let mut velocity = default_velocity;
        let mut bc = 0.475;
        let mut mass = default_mass;
        let mut diameter = default_diameter;
        let mut target_distance = if units == UnitSystem::Imperial { 100.0 } else { 100.0 };
        let mut sight_height = if units == UnitSystem::Imperial { 2.0 } else { 50.0 };
        let mut muzzle_height = if units == UnitSystem::Imperial { 60.0 } else { 1500.0 }; // inches or mm
        let mut target_height = 0.0; // inches or mm
        let mut drag_model = "G1";

        // Parse arguments
        let mut i = 0;
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
                "--target-distance" => {
                    if i + 1 < args.len() {
                        target_distance = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid target distance"))?;
                        i += 1;
                    }
                }
                "--sight-height" => {
                    if i + 1 < args.len() {
                        sight_height = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid sight height"))?;
                        i += 1;
                    }
                }
                "--muzzle-height" => {
                    if i + 1 < args.len() {
                        muzzle_height = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid muzzle height"))?;
                        i += 1;
                    }
                }
                "--target-height" => {
                    if i + 1 < args.len() {
                        target_height = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid target height"))?;
                        i += 1;
                    }
                }
                "--drag-model" => {
                    if i + 1 < args.len() {
                        drag_model = args[i + 1];
                        i += 1;
                    }
                }
                _ => {}
            }
            i += 1;
        }

        // Build inputs
        let mut inputs = InternalBallisticInputs::default();
        
        // Convert units
        match units {
            UnitSystem::Imperial => {
                inputs.muzzle_velocity = velocity * 0.3048;
                inputs.bullet_mass = mass * 0.00006479891;
                inputs.bullet_diameter = diameter * 0.0254;
                inputs.sight_height = sight_height * 0.0254;
            }
            UnitSystem::Metric => {
                inputs.muzzle_velocity = velocity;
                inputs.bullet_mass = mass * 0.001;
                inputs.bullet_diameter = diameter * 0.001;
                inputs.sight_height = sight_height * 0.001;
            }
        }
        
        inputs.bc_value = bc;
        inputs.bc_type = DragModel::from_str(drag_model)
            .ok_or_else(|| JsValue::from_str("Invalid drag model"))?;

        let target_distance_m = match units {
            UnitSystem::Imperial => target_distance * 0.9144,
            UnitSystem::Metric => target_distance,
        };

        match calculate_zero_angle_with_conditions(
            inputs,
            target_distance_m,
            0.0,
            WindConditions::default(),
            AtmosphericConditions::default()
        ) {
            Ok(zero_angle) => {
                let zero_degrees = zero_angle * 180.0 / std::f64::consts::PI;
                let moa_adjustment = zero_degrees * 60.0;
                let mrad_adjustment = zero_angle * 1000.0;
                
                Ok(format!(
                    "Zero Calculation Results\n\
                     ========================\n\
                     Target Distance: {} {}\n\
                     Zero Angle: {:.4}°\n\
                     MOA Adjustment: {:.2} MOA up\n\
                     Mrad Adjustment: {:.2} mrad up\n\
                     Sight Height: {} {}\n",
                    target_distance,
                    if units == UnitSystem::Imperial { "yards" } else { "meters" },
                    zero_degrees,
                    moa_adjustment,
                    mrad_adjustment,
                    sight_height,
                    if units == UnitSystem::Imperial { "inches" } else { "mm" }
                ))
            }
            Err(e) => Ok(format!("Error calculating zero: {}", e))
        }
    }

    fn handle_monte_carlo_command(&self, args: &[&str], units: UnitSystem) -> Result<String, JsValue> {
        // Default values
        let (default_velocity, default_mass, default_diameter) = match units {
            UnitSystem::Imperial => (2700.0, 168.0, 0.308),
            UnitSystem::Metric => (823.0, 10.9, 7.82),
        };

        let mut velocity = default_velocity;
        let mut angle = 0.0;
        let mut bc = 0.475;
        let mut mass = default_mass;
        let mut diameter = default_diameter;
        let mut num_sims = 1000;
        let mut velocity_std = if units == UnitSystem::Imperial { 10.0 } else { 3.0 };
        let mut angle_std = 0.1;
        let mut bc_std = 0.01;
        let mut wind_speed_std = if units == UnitSystem::Imperial { 2.0 } else { 0.5 };
        let mut wind_dir_std = 5.0;

        // Parse arguments
        let mut i = 0;
        while i < args.len() {
            match args[i] {
                "-v" | "--velocity" => {
                    if i + 1 < args.len() {
                        velocity = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid velocity"))?;
                        i += 1;
                    }
                }
                "-a" | "--angle" => {
                    if i + 1 < args.len() {
                        angle = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid angle"))?;
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
                "-n" | "--num-sims" => {
                    if i + 1 < args.len() {
                        num_sims = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid number of simulations"))?;
                        i += 1;
                    }
                }
                "--velocity-std" => {
                    if i + 1 < args.len() {
                        velocity_std = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid velocity std"))?;
                        i += 1;
                    }
                }
                "--angle-std" => {
                    if i + 1 < args.len() {
                        angle_std = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid angle std"))?;
                        i += 1;
                    }
                }
                "--bc-std" => {
                    if i + 1 < args.len() {
                        bc_std = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid BC std"))?;
                        i += 1;
                    }
                }
                "--wind-speed-std" => {
                    if i + 1 < args.len() {
                        wind_speed_std = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid wind speed std"))?;
                        i += 1;
                    }
                }
                "--wind-dir-std" => {
                    if i + 1 < args.len() {
                        wind_dir_std = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid wind direction std"))?;
                        i += 1;
                    }
                }
                _ => {}
            }
            i += 1;
        }

        // Build inputs
        let mut inputs = InternalBallisticInputs::default();
        
        // Convert units
        match units {
            UnitSystem::Imperial => {
                inputs.muzzle_velocity = velocity * 0.3048;
                inputs.bullet_mass = mass * 0.00006479891;
                inputs.bullet_diameter = diameter * 0.0254;
            }
            UnitSystem::Metric => {
                inputs.muzzle_velocity = velocity;
                inputs.bullet_mass = mass * 0.001;
                inputs.bullet_diameter = diameter * 0.001;
            }
        }
        
        inputs.bc_value = bc;
        inputs.muzzle_angle = angle * std::f64::consts::PI / 180.0;

        // Create Monte Carlo parameters
        let params = MonteCarloParams {
            num_simulations: num_sims,
            velocity_std_dev: velocity_std * (if units == UnitSystem::Imperial { 0.3048 } else { 1.0 }),
            angle_std_dev: angle_std * std::f64::consts::PI / 180.0,
            bc_std_dev: bc_std,
            wind_speed_std_dev: wind_speed_std * (if units == UnitSystem::Imperial { 0.44704 } else { 1.0 }),
            base_wind_speed: 0.0,
            base_wind_direction: 0.0,
            azimuth_std_dev: 0.001,
            target_distance: None,
        };

        match run_monte_carlo(inputs, params) {
            Ok(results) => {
                // Calculate statistics
                let mean_range: f64 = results.ranges.iter().sum::<f64>() / results.ranges.len() as f64;
                let mean_velocity: f64 = results.impact_velocities.iter().sum::<f64>() / results.impact_velocities.len() as f64;
                
                let range_std = {
                    let variance: f64 = results.ranges.iter()
                        .map(|r| (r - mean_range).powi(2))
                        .sum::<f64>() / results.ranges.len() as f64;
                    variance.sqrt()
                };
                
                let velocity_std = {
                    let variance: f64 = results.impact_velocities.iter()
                        .map(|v| (v - mean_velocity).powi(2))
                        .sum::<f64>() / results.impact_velocities.len() as f64;
                    variance.sqrt()
                };

                // Convert back to display units
                let (range_unit, velocity_unit) = match units {
                    UnitSystem::Imperial => ("yards", "fps"),
                    UnitSystem::Metric => ("meters", "m/s"),
                };
                
                let display_mean_range = match units {
                    UnitSystem::Imperial => mean_range * 1.09361,
                    UnitSystem::Metric => mean_range,
                };
                
                let display_range_std = match units {
                    UnitSystem::Imperial => range_std * 1.09361,
                    UnitSystem::Metric => range_std,
                };
                
                let display_mean_velocity = match units {
                    UnitSystem::Imperial => mean_velocity * 3.28084,
                    UnitSystem::Metric => mean_velocity,
                };
                
                let display_velocity_std = match units {
                    UnitSystem::Imperial => velocity_std * 3.28084,
                    UnitSystem::Metric => velocity_std,
                };

                Ok(format!(
                    "Monte Carlo Simulation Results\n\
                     ==============================\n\
                     Simulations Run: {}\n\n\
                     Range Statistics:\n\
                     - Mean: {:.1} {}\n\
                     - Std Dev: {:.1} {}\n\
                     - Min: {:.1} {}\n\
                     - Max: {:.1} {}\n\n\
                     Impact Velocity Statistics:\n\
                     - Mean: {:.0} {}\n\
                     - Std Dev: {:.0} {}\n\
                     - Min: {:.0} {}\n\
                     - Max: {:.0} {}",
                    num_sims,
                    display_mean_range, range_unit,
                    display_range_std, range_unit,
                    results.ranges.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap() * 
                        (if units == UnitSystem::Imperial { 1.09361 } else { 1.0 }), range_unit,
                    results.ranges.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap() * 
                        (if units == UnitSystem::Imperial { 1.09361 } else { 1.0 }), range_unit,
                    display_mean_velocity, velocity_unit,
                    display_velocity_std, velocity_unit,
                    results.impact_velocities.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap() * 
                        (if units == UnitSystem::Imperial { 3.28084 } else { 1.0 }), velocity_unit,
                    results.impact_velocities.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap() * 
                        (if units == UnitSystem::Imperial { 3.28084 } else { 1.0 }), velocity_unit,
                ))
            }
            Err(e) => Ok(format!("Error running Monte Carlo simulation: {}", e))
        }
    }

    fn handle_estimate_bc_command(&self, args: &[&str], units: UnitSystem) -> Result<String, JsValue> {
        let (default_velocity, default_mass, default_diameter) = match units {
            UnitSystem::Imperial => (2700.0, 168.0, 0.308),
            UnitSystem::Metric => (823.0, 10.9, 7.82),
        };

        let mut velocity = default_velocity;
        let mut mass = default_mass;
        let mut diameter = default_diameter;
        let mut data_points: Vec<(f64, f64)> = Vec::new();

        // Parse arguments
        let mut i = 0;
        while i < args.len() {
            match args[i] {
                "-v" | "--velocity" => {
                    if i + 1 < args.len() {
                        velocity = args[i + 1].parse()
                            .map_err(|_| JsValue::from_str("Invalid velocity"))?;
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
                "--data" => {
                    // Parse distance,drop pairs
                    if i + 1 < args.len() {
                        // Remove quotes if present
                        let data_str = args[i + 1].trim_matches('\'').trim_matches('"');
                        let pairs: Vec<&str> = data_str.split(';').collect();
                        for pair in pairs {
                            let parts: Vec<&str> = pair.split(',').collect();
                            if parts.len() == 2 {
                                let distance: f64 = parts[0].trim().parse()
                                    .map_err(|e| JsValue::from_str(&format!("Invalid distance '{}': {}", parts[0], e)))?;
                                let drop: f64 = parts[1].trim().parse()
                                    .map_err(|e| JsValue::from_str(&format!("Invalid drop '{}': {}", parts[1], e)))?;
                                data_points.push((distance, drop));
                            }
                        }
                        i += 1;
                    }
                }
                _ => {}
            }
            i += 1;
        }

        if data_points.is_empty() {
            return Ok("Error: No trajectory data provided. Use --data with distance,drop pairs separated by semicolons.\nExample: --data \"100,2.5;200,10.2;300,23.5\"".to_string());
        }

        // Convert units
        let velocity_mps = match units {
            UnitSystem::Imperial => velocity * 0.3048,
            UnitSystem::Metric => velocity,
        };
        
        let mass_kg = match units {
            UnitSystem::Imperial => mass * 0.00006479891,
            UnitSystem::Metric => mass * 0.001,
        };
        
        let diameter_m = match units {
            UnitSystem::Imperial => diameter * 0.0254,
            UnitSystem::Metric => diameter * 0.001,
        };

        // Convert data points to meters
        let metric_points: Vec<(f64, f64)> = data_points.iter().map(|(dist, drop)| {
            match units {
                UnitSystem::Imperial => (*dist * 0.9144, *drop * 0.0254),
                UnitSystem::Metric => (*dist, *drop * 0.001),
            }
        }).collect();

        match estimate_bc_from_trajectory(velocity_mps, mass_kg, diameter_m, &metric_points) {
            Ok(estimated_bc) => {
                Ok(format!(
                    "BC Estimation Results\n\
                     ====================\n\
                     Estimated BC: {:.3}\n\
                     Based on {} data points\n\
                     Velocity: {} {}\n\
                     Mass: {} {}\n\
                     Diameter: {} {}",
                    estimated_bc,
                    data_points.len(),
                    velocity, if units == UnitSystem::Imperial { "fps" } else { "m/s" },
                    mass, if units == UnitSystem::Imperial { "grains" } else { "grams" },
                    diameter, if units == UnitSystem::Imperial { "inches" } else { "mm" }
                ))
            }
            Err(e) => Ok(format!("Error estimating BC: {}", e))
        }
    }

    fn format_trajectory_table(&self, result: &crate::cli_api::TrajectoryResult, zero_distance: Option<f64>, units: UnitSystem, full: bool) -> String {
        let mut output = String::new();
        output.push_str("Trajectory Calculation Results\n");
        output.push_str("==============================\n\n");
        output.push_str("Range | Drop | Drift | Velocity | Energy | Time\n");
        output.push_str("------|------|-------|----------|--------|------\n");

        // Determine sampling interval based on max range and full flag
        let max_range_display = match units {
            UnitSystem::Imperial => result.max_range * 1.09361, // m to yards
            UnitSystem::Metric => result.max_range,
        };
        
        let sample_interval = if full {
            if max_range_display < 100.0 { 10.0 } else { 25.0 }
        } else {
            if max_range_display < 500.0 { 50.0 } else { 100.0 }
        };
        
        let mut current_range = 0.0;
        
        // Get initial height for reference (sight height)
        let sight_height = if !result.points.is_empty() {
            result.points[0].position.y
        } else {
            0.05 // Default 2 inches
        };
        
        for (idx, point) in result.points.iter().enumerate() {
            let range_display = match units {
                UnitSystem::Imperial => point.position.x * 1.09361, // X is downrange, m to yards
                UnitSystem::Metric => point.position.x,  // X is downrange
            };
            
            let is_last_point = idx == result.points.len() - 1;
            
            // Show point if it's at the sampling interval OR if it's the last point OR if it's the zero distance
            let should_show = range_display >= current_range || 
                            is_last_point || 
                            (zero_distance.is_some() && (range_display - zero_distance.unwrap()).abs() < 1.0);
            
            if should_show {
                let drop = sight_height - point.position.y;
                let drift = point.position.z;  // Z is lateral (windage)
                let velocity = point.velocity_magnitude;
                
                // Format values based on unit system
                let (range_str, drop_str, drift_str, velocity_str, energy_str) = match units {
                    UnitSystem::Imperial => {
                        let energy_ftlb = 0.5 * (result.points[0].kinetic_energy / 0.5) * 
                                        (velocity / result.points[0].velocity_magnitude).powi(2) * 0.737562149;
                        (
                            format!("{:.0} yd", range_display),
                            format!("{:.1} in", drop * 39.3701),
                            format!("{:.1} in", drift * 39.3701),
                            format!("{:.0} fps", velocity * 3.28084),
                            format!("{:.0} ft-lb", energy_ftlb)
                        )
                    },
                    UnitSystem::Metric => {
                        (
                            format!("{:.0} m", range_display),
                            format!("{:.1} cm", drop * 100.0),
                            format!("{:.1} cm", drift * 100.0),
                            format!("{:.0} m/s", velocity),
                            format!("{:.0} J", point.kinetic_energy)
                        )
                    }
                };
                
                output.push_str(&format!(
                    "{:6} | {:6} | {:7} | {:10} | {:8} | {:.3} s\n",
                    range_str, drop_str, drift_str, velocity_str, energy_str, point.time
                ));
                
                if range_display >= current_range {
                    current_range += sample_interval;
                }
            }
        }
        
        // Add summary
        output.push_str(&format!("\nMax Range: {:.0} {}\n", 
            if units == UnitSystem::Imperial { result.max_range * 1.09361 } else { result.max_range },
            if units == UnitSystem::Imperial { "yards" } else { "meters" }
        ));
        
        // Check termination reason
        if let Some(last_point) = result.points.last() {
            // Debug info about last point
            let last_height = last_point.position.y;
            let last_range = last_point.position.x;  // X is downrange
            let last_time = last_point.time;
            
            // Ground threshold is typically around 0.01m (1cm), so if y is close to or below that, it hit ground
            if last_height <= 0.01 {
                output.push_str(&format!("Bullet struck ground at {:.0} {}\n",
                    if units == UnitSystem::Imperial { last_range * 1.09361 } else { last_range },
                    if units == UnitSystem::Imperial { "yards" } else { "meters" }
                ));
            } else if last_time >= 99.0 {
                output.push_str(&format!("Trajectory timeout at {:.0} {} (time limit reached)\n",
                    if units == UnitSystem::Imperial { last_range * 1.09361 } else { last_range },
                    if units == UnitSystem::Imperial { "yards" } else { "meters" }
                ));
            }
        }
        
        output.push_str(&format!("Max Height: {:.1} {}\n", 
            if units == UnitSystem::Imperial { result.max_height * 39.3701 } else { result.max_height * 100.0 },
            if units == UnitSystem::Imperial { "inches" } else { "cm" }
        ));
        output.push_str(&format!("Time of Flight: {:.2} seconds\n", result.time_of_flight));
        output.push_str(&format!("Impact Velocity: {:.0} {}\n", 
            if units == UnitSystem::Imperial { result.impact_velocity * 3.28084 } else { result.impact_velocity },
            if units == UnitSystem::Imperial { "fps" } else { "m/s" }
        ));

        output
    }

    fn format_trajectory_json(&self, result: &crate::cli_api::TrajectoryResult, units: UnitSystem) -> String {
        let points: Vec<serde_json::Value> = result.points.iter().map(|p| {
            match units {
                UnitSystem::Imperial => {
                    serde_json::json!({
                        "range_yards": p.position.x * 1.09361,  // X is downrange
                        "drop_inches": (result.points[0].position.y - p.position.y) * 39.3701,
                        "drift_inches": p.position.z * 39.3701,  // Z is lateral (windage)
                        "velocity_fps": p.velocity_magnitude * 3.28084,
                        "energy_ftlb": p.kinetic_energy * 0.737562149,
                        "time_seconds": p.time
                    })
                },
                UnitSystem::Metric => {
                    serde_json::json!({
                        "range_meters": p.position.x,  // X is downrange
                        "drop_cm": (result.points[0].position.y - p.position.y) * 100.0,
                        "drift_cm": p.position.z * 100.0,  // Z is lateral (windage)
                        "velocity_mps": p.velocity_magnitude,
                        "energy_joules": p.kinetic_energy,
                        "time_seconds": p.time
                    })
                }
            }
        }).collect();

        let summary = match units {
            UnitSystem::Imperial => {
                serde_json::json!({
                    "max_range_yards": result.max_range * 1.09361,
                    "max_height_inches": result.max_height * 39.3701,
                    "time_of_flight_seconds": result.time_of_flight,
                    "impact_velocity_fps": result.impact_velocity * 3.28084
                })
            },
            UnitSystem::Metric => {
                serde_json::json!({
                    "max_range_meters": result.max_range,
                    "max_height_cm": result.max_height * 100.0,
                    "time_of_flight_seconds": result.time_of_flight,
                    "impact_velocity_mps": result.impact_velocity
                })
            }
        };

        let output = serde_json::json!({
            "trajectory": points,
            "summary": summary
        });

        serde_json::to_string_pretty(&output).unwrap_or_else(|_| "Error formatting JSON".to_string())
    }

    fn format_trajectory_csv(&self, result: &crate::cli_api::TrajectoryResult, units: UnitSystem, full: bool) -> String {
        let mut output = String::new();
        
        // Header
        match units {
            UnitSystem::Imperial => {
                output.push_str("Range(yards),Drop(inches),Drift(inches),Velocity(fps),Energy(ft-lb),Time(seconds)\n");
            },
            UnitSystem::Metric => {
                output.push_str("Range(meters),Drop(cm),Drift(cm),Velocity(m/s),Energy(joules),Time(seconds)\n");
            }
        }

        // Determine sampling interval
        let max_range_display = match units {
            UnitSystem::Imperial => result.max_range * 1.09361,
            UnitSystem::Metric => result.max_range,
        };
        
        let sample_interval = if full {
            if max_range_display < 100.0 { 10.0 } else { 25.0 }
        } else {
            if max_range_display < 500.0 { 50.0 } else { 100.0 }
        };
        
        let mut current_range = 0.0;
        let sight_height = if !result.points.is_empty() {
            result.points[0].position.y
        } else {
            0.05
        };

        for (idx, point) in result.points.iter().enumerate() {
            let range_display = match units {
                UnitSystem::Imperial => point.position.x * 1.09361,  // X is downrange
                UnitSystem::Metric => point.position.x,  // X is downrange
            };

            let is_last_point = idx == result.points.len() - 1;

            if range_display >= current_range || is_last_point {
                let drop = sight_height - point.position.y;

                match units {
                    UnitSystem::Imperial => {
                        let energy_ftlb = point.kinetic_energy * 0.737562149;
                        output.push_str(&format!(
                            "{:.1},{:.2},{:.2},{:.0},{:.0},{:.3}\n",
                            range_display,
                            drop * 39.3701,
                            point.position.z * 39.3701,  // Z is lateral (windage)
                            point.velocity_magnitude * 3.28084,
                            energy_ftlb,
                            point.time
                        ));
                    },
                    UnitSystem::Metric => {
                        output.push_str(&format!(
                            "{:.1},{:.2},{:.2},{:.0},{:.0},{:.3}\n",
                            range_display,
                            drop * 100.0,
                            point.position.z * 100.0,  // Z is lateral (windage)
                            point.velocity_magnitude,
                            point.kinetic_energy,
                            point.time
                        ));
                    }
                }
                
                if range_display >= current_range {
                    current_range += sample_interval;
                }
            }
        }

        output
    }

    fn show_help(&self) -> String {
        r#"Ballistics Engine - WebAssembly Version

Usage: ballistics [OPTIONS] <COMMAND>

Commands:
  trajectory      Calculate ballistic trajectory
  zero           Calculate sight adjustment for zero
  monte-carlo    Run Monte Carlo simulation
  estimate-bc    Estimate BC from trajectory data
  help           Show this help message

Global Options:
  -u, --units <SYSTEM>  Unit system (imperial/metric) [default: imperial]

Trajectory Command:
  ballistics trajectory [OPTIONS]
  
  Basic Options:
    -v, --velocity <VEL>         Muzzle velocity (fps/m/s)
    -b, --bc <BC>                Ballistic coefficient
    -m, --mass <MASS>            Mass (grains/grams)
    -d, --diameter <DIA>         Diameter (inches/mm)
    -a, --angle <ANGLE>          Launch angle (degrees)
    --drag-model <MODEL>         Drag model (G1/G7)
    --max-range <RANGE>          Maximum range (yards/meters)
    -z, --auto-zero <DIST>       Auto-zero at distance
    -o, --output <FORMAT>        Output format (table/json/csv)
    --full                       Show all trajectory points
    
  Environmental:
    --wind-speed <SPEED>         Wind speed (mph/m/s)
    --wind-direction <DIR>       Wind direction (degrees)
    --temperature <TEMP>         Temperature (F/C)
    --pressure <PRESSURE>        Pressure (inHg/hPa)
    --humidity <HUMIDITY>        Humidity (0-100%)
    --altitude <ALT>             Altitude (feet/meters)
    
  Advanced Physics:
    --enable-magnus              Enable Magnus effect
    --enable-coriolis            Enable Coriolis effect
    --enable-spin-drift          Enable enhanced spin drift
    --enable-wind-shear          Enable altitude-dependent wind
    --enable-pitch-damping       Enable transonic stability
    --enable-precession          Enable angular motion physics
    --use-euler                  Use Euler integration (less accurate)
    --use-rk4-fixed              Use fixed-step RK4 (default: adaptive RK45)
    --sample-trajectory          Enable trajectory sampling
    --use-bc-segments            Use velocity-based BC
    --use-powder-sensitivity     Enable powder temp sensitivity
    
  Additional Parameters:
    --twist-rate <RATE>          Barrel twist (inches per turn)
    --twist-right <BOOL>         Right-hand twist (true/false)
    --latitude <LAT>             Latitude for Coriolis (degrees)
    --shooting-angle <ANGLE>     Uphill/downhill angle (degrees)
    --sight-height <HEIGHT>      Sight height above bore (inches/mm)
    --muzzle-height <HEIGHT>     Shooter height above ground (inches/mm)
    --target-height <HEIGHT>     Target height above ground (inches/mm)
    --powder-temp <TEMP>         Powder temperature
    --powder-temp-sensitivity    Velocity change per degree

Zero Command:
  ballistics zero [OPTIONS]
  
  Options:
    -v, --velocity <VEL>         Muzzle velocity
    -b, --bc <BC>                Ballistic coefficient
    -m, --mass <MASS>            Mass
    -d, --diameter <DIA>         Diameter
    --target-distance <DIST>     Target distance for zero
    --sight-height <HEIGHT>      Sight height above bore

Monte Carlo Command:
  ballistics monte-carlo [OPTIONS]
  
  Options:
    -v, --velocity <VEL>         Base velocity
    -b, --bc <BC>                Base BC
    -m, --mass <MASS>            Mass
    -d, --diameter <DIA>         Diameter
    -n, --num-sims <N>           Number of simulations
    --velocity-std <STD>         Velocity std deviation
    --angle-std <STD>            Angle std deviation
    --bc-std <STD>               BC std deviation
    --wind-speed-std <STD>       Wind speed std deviation
    --wind-dir-std <STD>         Wind direction std deviation

Examples:
  ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308
  ballistics trajectory --auto-zero 200 --enable-spin-drift
  ballistics --units metric trajectory -v 823 -b 0.475 -m 10.9
  ballistics zero --target-distance 300
  ballistics monte-carlo -n 1000 --velocity-std 10"#.to_string()
    }
}

// ============================================================================
// Object-Oriented API
// ============================================================================

/// Object-oriented calculator for programmatic use
/// Provides a type-safe, fluent API alternative to the CLI interface
#[wasm_bindgen]
pub struct Calculator {
    // Projectile properties
    bc: f64,
    velocity_fps: f64,
    mass_grains: f64,
    diameter_inches: f64,
    drag_model: String,

    // Environmental conditions
    wind_speed_mph: f64,
    wind_direction_deg: f64,
    temperature_f: f64,
    pressure_inhg: f64,
    humidity_percent: f64,
    altitude_ft: f64,

    // Shooting parameters
    sight_height_inches: f64,
    zero_range_yards: Option<f64>,
    max_range_yards: f64,

    // Advanced options
    enable_spin_drift: bool,
    enable_coriolis: bool,
    twist_rate_inches: Option<f64>,
    latitude_deg: Option<f64>,
}

#[wasm_bindgen]
impl Calculator {
    /// Create a new calculator with default values
    /// Defaults: .308 Winchester 168gr at 2700 fps, standard atmosphere
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Calculator {
            bc: 0.475,
            velocity_fps: 2700.0,
            mass_grains: 168.0,
            diameter_inches: 0.308,
            drag_model: "G7".to_string(),

            wind_speed_mph: 0.0,
            wind_direction_deg: 90.0,
            temperature_f: 59.0,
            pressure_inhg: 29.92,
            humidity_percent: 50.0,
            altitude_ft: 0.0,

            sight_height_inches: 2.0,
            zero_range_yards: None,
            max_range_yards: 1000.0,

            enable_spin_drift: false,
            enable_coriolis: false,
            twist_rate_inches: None,
            latitude_deg: None,
        }
    }

    // Projectile property setters (fluent API)

    #[wasm_bindgen(js_name = setBC)]
    pub fn set_bc(mut self, bc: f64) -> Self {
        self.bc = bc;
        self
    }

    #[wasm_bindgen(js_name = setVelocity)]
    pub fn set_velocity(mut self, velocity_fps: f64) -> Self {
        self.velocity_fps = velocity_fps;
        self
    }

    #[wasm_bindgen(js_name = setMass)]
    pub fn set_mass(mut self, mass_grains: f64) -> Self {
        self.mass_grains = mass_grains;
        self
    }

    #[wasm_bindgen(js_name = setDiameter)]
    pub fn set_diameter(mut self, diameter_inches: f64) -> Self {
        self.diameter_inches = diameter_inches;
        self
    }

    #[wasm_bindgen(js_name = setDragModel)]
    pub fn set_drag_model(mut self, model: &str) -> Self {
        self.drag_model = model.to_string();
        self
    }

    // Environmental setters

    #[wasm_bindgen(js_name = setWind)]
    pub fn set_wind(mut self, speed_mph: f64, direction_deg: f64) -> Self {
        self.wind_speed_mph = speed_mph;
        self.wind_direction_deg = direction_deg;
        self
    }

    #[wasm_bindgen(js_name = setTemperature)]
    pub fn set_temperature(mut self, temp_f: f64) -> Self {
        self.temperature_f = temp_f;
        self
    }

    #[wasm_bindgen(js_name = setPressure)]
    pub fn set_pressure(mut self, pressure_inhg: f64) -> Self {
        self.pressure_inhg = pressure_inhg;
        self
    }

    #[wasm_bindgen(js_name = setHumidity)]
    pub fn set_humidity(mut self, humidity: f64) -> Self {
        self.humidity_percent = humidity;
        self
    }

    #[wasm_bindgen(js_name = setAltitude)]
    pub fn set_altitude(mut self, altitude_ft: f64) -> Self {
        self.altitude_ft = altitude_ft;
        self
    }

    // Shooting parameter setters

    #[wasm_bindgen(js_name = setSightHeight)]
    pub fn set_sight_height(mut self, height_inches: f64) -> Self {
        self.sight_height_inches = height_inches;
        self
    }

    #[wasm_bindgen(js_name = setZeroRange)]
    pub fn set_zero_range(mut self, range_yards: f64) -> Self {
        self.zero_range_yards = Some(range_yards);
        self
    }

    #[wasm_bindgen(js_name = setMaxRange)]
    pub fn set_max_range(mut self, range_yards: f64) -> Self {
        self.max_range_yards = range_yards;
        self
    }

    // Advanced option setters

    #[wasm_bindgen(js_name = enableSpinDrift)]
    pub fn enable_spin_drift_opt(mut self, enabled: bool, twist_rate: Option<f64>) -> Self {
        self.enable_spin_drift = enabled;
        self.twist_rate_inches = twist_rate;
        self
    }

    #[wasm_bindgen(js_name = enableCoriolis)]
    pub fn enable_coriolis_opt(mut self, enabled: bool, latitude: Option<f64>) -> Self {
        self.enable_coriolis = enabled;
        self.latitude_deg = latitude;
        self
    }

    // Calculation method

    /// Calculate trajectory and return result as JavaScript object
    /// Returns: { range_yards, drop_inches, windage_inches, velocity_fps, energy_ftlb, time_sec }
    #[wasm_bindgen(js_name = calculateTrajectory)]
    pub fn calculate_trajectory(&self, range_yards: f64) -> Result<JsValue, JsValue> {
        // Build CLI command from parameters
        let mut cmd = format!(
            "ballistics trajectory -v {} -b {} -m {} -d {} --drag-model {} --max-range {}",
            self.velocity_fps,
            self.bc,
            self.mass_grains,
            self.diameter_inches,
            self.drag_model,
            range_yards
        );

        // Add environmental parameters
        if self.wind_speed_mph > 0.0 {
            cmd.push_str(&format!(" --wind-speed {} --wind-direction {}",
                self.wind_speed_mph, self.wind_direction_deg));
        }

        cmd.push_str(&format!(" --temperature {} --pressure {} --humidity {} --altitude {}",
            self.temperature_f, self.pressure_inhg, self.humidity_percent, self.altitude_ft));

        // Add shooting parameters
        cmd.push_str(&format!(" --sight-height {}", self.sight_height_inches));

        if let Some(zero) = self.zero_range_yards {
            cmd.push_str(&format!(" --auto-zero {}", zero));
        }

        // Add advanced options
        if self.enable_spin_drift {
            cmd.push_str(" --enable-spin-drift");
            if let Some(twist) = self.twist_rate_inches {
                cmd.push_str(&format!(" --twist-rate {}", twist));
            }
        }

        if self.enable_coriolis {
            cmd.push_str(" --enable-coriolis");
            if let Some(lat) = self.latitude_deg {
                cmd.push_str(&format!(" --latitude {}", lat));
            }
        }

        // Use JSON output format for easy parsing
        cmd.push_str(" -o json");

        // Execute command
        let wasm_ballistics = WasmBallistics::new();
        let result_str = wasm_ballistics.run_command(&cmd)?;

        // Strip any text before the JSON (like zero info messages)
        let json_str = if let Some(json_start) = result_str.find('{') {
            &result_str[json_start..]
        } else {
            &result_str
        };

        // Parse JSON result
        let json_result: serde_json::Value = serde_json::from_str(json_str)
            .map_err(|e| JsValue::from_str(&format!("JSON parse error: {}. Result was: {}", e, &json_str[..json_str.len().min(500)])))?;

        // Find the point closest to the requested range
        if let Some(trajectory) = json_result.get("trajectory").and_then(|t| t.as_array()) {
            let target_point = trajectory.iter()
                .min_by_key(|p| {
                    let range = p.get("range_yards").and_then(|r| r.as_f64()).unwrap_or(0.0);
                    ((range - range_yards).abs() * 100.0) as i64
                })
                .ok_or_else(|| JsValue::from_str("No trajectory points found"))?;

            // Manually construct JavaScript object to avoid Map conversion
            let result = js_sys::Object::new();
            js_sys::Reflect::set(&result, &"range_yards".into(),
                &target_point.get("range_yards").and_then(|v| v.as_f64()).unwrap_or(0.0).into())?;
            js_sys::Reflect::set(&result, &"drop_inches".into(),
                &target_point.get("drop_inches").and_then(|v| v.as_f64()).unwrap_or(0.0).into())?;
            js_sys::Reflect::set(&result, &"drift_inches".into(),
                &target_point.get("drift_inches").and_then(|v| v.as_f64()).unwrap_or(0.0).into())?;
            js_sys::Reflect::set(&result, &"velocity_fps".into(),
                &target_point.get("velocity_fps").and_then(|v| v.as_f64()).unwrap_or(0.0).into())?;
            js_sys::Reflect::set(&result, &"energy_ftlb".into(),
                &target_point.get("energy_ftlb").and_then(|v| v.as_f64()).unwrap_or(0.0).into())?;
            js_sys::Reflect::set(&result, &"time_seconds".into(),
                &target_point.get("time_seconds").and_then(|v| v.as_f64()).unwrap_or(0.0).into())?;

            Ok(result.into())
        } else {
            Err(JsValue::from_str("Invalid trajectory data"))
        }
    }

    /// Get full trajectory table as array of points
    /// Returns array of: [{ range_yards, drop_inches, windage_inches, velocity_fps, energy_ftlb, time_sec }, ...]
    #[wasm_bindgen(js_name = getFullTrajectory)]
    pub fn get_full_trajectory(&self) -> Result<JsValue, JsValue> {
        // Build CLI command (similar to calculate_trajectory but get full table)
        let mut cmd = format!(
            "ballistics trajectory -v {} -b {} -m {} -d {} --drag-model {} --max-range {} --full",
            self.velocity_fps,
            self.bc,
            self.mass_grains,
            self.diameter_inches,
            self.drag_model,
            self.max_range_yards
        );

        // Add environmental parameters
        if self.wind_speed_mph > 0.0 {
            cmd.push_str(&format!(" --wind-speed {} --wind-direction {}",
                self.wind_speed_mph, self.wind_direction_deg));
        }

        cmd.push_str(&format!(" --temperature {} --pressure {} --humidity {} --altitude {}",
            self.temperature_f, self.pressure_inhg, self.humidity_percent, self.altitude_ft));

        cmd.push_str(&format!(" --sight-height {}", self.sight_height_inches));

        if let Some(zero) = self.zero_range_yards {
            cmd.push_str(&format!(" --auto-zero {}", zero));
        }

        if self.enable_spin_drift {
            cmd.push_str(" --enable-spin-drift");
            if let Some(twist) = self.twist_rate_inches {
                cmd.push_str(&format!(" --twist-rate {}", twist));
            }
        }

        if self.enable_coriolis {
            cmd.push_str(" --enable-coriolis");
            if let Some(lat) = self.latitude_deg {
                cmd.push_str(&format!(" --latitude {}", lat));
            }
        }

        cmd.push_str(" -o json");

        // Execute command
        let wasm_ballistics = WasmBallistics::new();
        let result_str = wasm_ballistics.run_command(&cmd)?;

        // Strip any text before the JSON (like zero info messages)
        let json_str = if let Some(json_start) = result_str.find('{') {
            &result_str[json_start..]
        } else {
            &result_str
        };

        // Parse JSON result and return trajectory array
        let json_result: serde_json::Value = serde_json::from_str(json_str)
            .map_err(|e| JsValue::from_str(&format!("JSON parse error: {}. Result: {}", e, &json_str[..json_str.len().min(500)])))?;

        if let Some(trajectory) = json_result.get("trajectory") {
            // Convert trajectory array to JavaScript array of plain objects
            let js_array = js_sys::Array::new();

            if let Some(points) = trajectory.as_array() {
                for point in points {
                    let js_point = js_sys::Object::new();

                    // Extract each field and add to JavaScript object
                    if let Some(range) = point.get("range_yards").and_then(|v| v.as_f64()) {
                        js_sys::Reflect::set(&js_point, &"range_yards".into(), &range.into())?;
                    }
                    if let Some(drop) = point.get("drop_inches").and_then(|v| v.as_f64()) {
                        js_sys::Reflect::set(&js_point, &"drop_inches".into(), &drop.into())?;
                    }
                    if let Some(windage) = point.get("windage_inches").and_then(|v| v.as_f64()) {
                        js_sys::Reflect::set(&js_point, &"windage_inches".into(), &windage.into())?;
                    }
                    if let Some(velocity) = point.get("velocity_fps").and_then(|v| v.as_f64()) {
                        js_sys::Reflect::set(&js_point, &"velocity_fps".into(), &velocity.into())?;
                    }
                    if let Some(energy) = point.get("energy_ftlb").and_then(|v| v.as_f64()) {
                        js_sys::Reflect::set(&js_point, &"energy_ftlb".into(), &energy.into())?;
                    }
                    if let Some(time) = point.get("time_sec").and_then(|v| v.as_f64()) {
                        js_sys::Reflect::set(&js_point, &"time_sec".into(), &time.into())?;
                    }

                    js_array.push(&js_point);
                }
            }

            Ok(js_array.into())
        } else {
            Err(JsValue::from_str("No trajectory data found"))
        }
    }
}