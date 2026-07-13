// WASM bindings for the ballistics engine with full CLI feature parity
use serde_json;
use wasm_bindgen::prelude::*;

use crate::bc_table_5d::Bc5dTable;
use crate::cli_api::{
    calculate_zero_angle_with_conditions, estimate_bc_fit, run_monte_carlo_with_direction_std_dev,
    AtmosphericConditions, BallisticInputs as InternalBallisticInputs, BcFitMode, MonteCarloParams,
    TrajectorySolver, WindConditions,
};
use crate::drag_model::DragModel;
use crate::moving_target::calculate_lead;
use std::cell::RefCell;

#[wasm_bindgen]
pub struct WasmBallistics {
    /// Optional BC5D correction table loaded from an in-memory `.bin`. When
    /// present, trajectory runs with `--use-bc-segments` synthesize
    /// velocity-dependent BC segments from it (offline parity with the online
    /// solver's ClusterBCDegradation + BC-segment + weather corrections).
    bc5d_table: RefCell<Option<Bc5dTable>>,
}

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
        WasmBallistics {
            bc5d_table: RefCell::new(None),
        }
    }

    /// Load a BC5D correction table from the raw bytes of a `bc5d_<caliber>.bin`
    /// file. The host (browser `fetch()` or Node `fs`/`fetch`) is responsible
    /// for retrieving the file — WASM has no filesystem or network — and passes
    /// the bytes here.
    ///
    /// Once loaded, any `trajectory` run that includes `--use-bc-segments` will
    /// apply velocity-dependent BC segments synthesized from this table. Load a
    /// table matching the bullet's caliber (e.g. `bc5d_308.bin` for a .308).
    ///
    /// Returns a short human-readable summary of the loaded table. Replaces any
    /// previously loaded table.
    #[wasm_bindgen(js_name = loadBc5dTable)]
    pub fn load_bc5d_table(&self, bytes: &[u8]) -> Result<String, JsValue> {
        let table = Bc5dTable::from_bytes(bytes)
            .map_err(|e| JsValue::from_str(&format!("Failed to parse BC5D table: {}", e)))?;
        let summary = format!(
            "Loaded BC5D table: caliber {:.3}\", {} cells, api_version {}",
            table.caliber(),
            table.total_cells(),
            table.api_version()
        );
        *self.bc5d_table.borrow_mut() = Some(table);
        Ok(summary)
    }

    /// Report whether a BC5D table is currently loaded.
    #[wasm_bindgen(js_name = hasBc5dTable)]
    pub fn has_bc5d_table(&self) -> bool {
        self.bc5d_table.borrow().is_some()
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

        // MBA-1289: `--units <system>` may appear anywhere — including BEFORE the command
        // word, the order the help text advertises. Units were already parsed from the
        // full list above, so strip the flag and its value here; otherwise the dispatch
        // below would read `--units` as the command and fail with "Unknown command".
        let mut stripped: Vec<&str> = Vec::with_capacity(args.len());
        let mut i = 0;
        while i < args.len() {
            if args[i] == "--units" || args[i] == "-u" {
                i += 2; // skip the flag and its value (a dangling flag skips harmlessly)
            } else {
                stripped.push(args[i]);
                i += 1;
            }
        }
        let args = stripped;
        if args.is_empty() {
            return Ok(self.show_help());
        }

        // Route to appropriate command handler
        match args[0] {
            "version" => Ok(format!(
                "Ballistics Engine v{}\nWASM Build\n",
                env!("CARGO_PKG_VERSION")
            )),
            "trajectory" => self.handle_trajectory_command(&args[1..], units),
            "zero" => self.handle_zero_command(&args[1..], units),
            "monte-carlo" | "montecarlo" => self.handle_monte_carlo_command(&args[1..], units),
            "estimate-bc" => self.handle_estimate_bc_command(&args[1..], units),
            "lead" => self.handle_lead_command(&args[1..], units),
            _ => Ok(format!(
                "Error: Unknown command '{}'\n\n{}",
                args[0],
                self.show_help()
            )),
        }
    }

    fn handle_trajectory_command(
        &self,
        args: &[&str],
        units: UnitSystem,
    ) -> Result<String, JsValue> {
        // Default values based on unit system
        let (default_velocity, default_mass, default_diameter, default_temp, default_pressure) =
            match units {
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
        let mut max_range = if units == UnitSystem::Imperial {
            1000.0
        } else {
            914.4
        };
        let mut time_step = 0.001;
        let mut wind_speed = 0.0;
        // f64 is required: `wind_direction.to_radians()` below needs a known receiver
        // type (method resolution can't infer it from the field assignment alone).
        let mut wind_direction: f64 = 0.0;
        // Raw "SPEED:ANGLE:UNTIL" strings; every --wind-segment occurrence is collected
        // (the parse loop visits all args, so repeats accumulate here).
        let mut wind_segment_strs: Vec<String> = Vec::new();
        // Raw "VMIN:VMAX:BC" strings from every --bc-segment occurrence (manual
        // velocity-keyed BC segments; take precedence over a loaded BC5D table).
        let mut bc_segment_strs: Vec<String> = Vec::new();
        let mut temperature = default_temp;
        let mut pressure = default_pressure;
        let mut humidity = 50.0;
        let mut altitude = 0.0;
        let mut output_format = OutputFormat::Table;
        let mut full = false;
        let mut auto_zero: Option<f64> = None;
        let mut sight_height = if units == UnitSystem::Imperial {
            2.0
        } else {
            50.0
        }; // inches or mm
        let mut muzzle_height = if units == UnitSystem::Imperial {
            60.0
        } else {
            1500.0
        }; // inches or mm (standing)
        let mut target_height = 0.0; // inches or mm (ground level)

        // Advanced physics flags
        let mut enable_magnus = false;
        let mut enable_coriolis = false;
        let mut use_euler = false;
        let mut use_rk4_fixed = false; // Use fixed-step RK4 instead of adaptive RK45
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
        let mut shot_direction: Option<f64> = None; // compass bearing, degrees, 0=N (Coriolis)
        let mut shooting_angle = 0.0;
        let mut cant_angle_deg = 0.0;
        let mut powder_temp_sensitivity = if units == UnitSystem::Imperial {
            1.0
        } else {
            0.3048 / (5.0 / 9.0)
        };
        let mut powder_temp = if units == UnitSystem::Imperial {
            crate::constants::DEFAULT_POWDER_REFERENCE_TEMP_F
        } else {
            crate::constants::DEFAULT_POWDER_REFERENCE_TEMP_C
        };
        // Optional measured powder-temperature -> velocity curve ("TEMP:VEL,..."),
        // parsed after unit resolution. Supersedes the linear sensitivity model.
        let mut powder_temp_curve_str: Option<String> = None;
        // Track whether --powder-temp was explicitly given. When a curve is present it
        // becomes the powder temperature the curve is interpolated at (decoupled from the
        // air temperature); when not given, the curve falls back to the air temperature.
        let mut powder_temp_provided = false;

        // Zero-day condition overrides. When --auto-zero is used, these let the zero
        // ANGLE be solved under the conditions the rifle was actually zeroed in (a
        // different day: air temperature, pressure, humidity, altitude, and — via the
        // caller's own powder-temp/velocity table — muzzle velocity), while the
        // trajectory itself runs under the current shot-day conditions. Omitting all zero-day
        // flags reuses the shot-day values exactly; coupled powder fallbacks are resolved below.
        let mut zero_velocity: Option<f64> = None;
        let mut zero_temperature: Option<f64> = None;
        let mut zero_pressure: Option<f64> = None;
        let mut zero_humidity: Option<f64> = None;
        let mut zero_altitude: Option<f64> = None;
        let mut zero_powder_temp: Option<f64> = None;

        // Parse arguments
        let mut i = 0;
        while i < args.len() {
            match args[i] {
                "-v" | "--velocity" => {
                    if i + 1 < args.len() {
                        velocity = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid velocity"))?;
                        i += 1;
                    }
                }
                "-a" | "--angle" => {
                    if i + 1 < args.len() {
                        angle = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid angle"))?;
                        i += 1;
                    }
                }
                "-b" | "--bc" => {
                    if i + 1 < args.len() {
                        bc = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid BC"))?;
                        i += 1;
                    }
                }
                "-m" | "--mass" => {
                    if i + 1 < args.len() {
                        mass = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid mass"))?;
                        i += 1;
                    }
                }
                "-d" | "--diameter" => {
                    if i + 1 < args.len() {
                        diameter = args[i + 1]
                            .parse()
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
                        max_range = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid max range"))?;
                        i += 1;
                    }
                }
                "--time-step" => {
                    if i + 1 < args.len() {
                        time_step = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid time step"))?;
                        i += 1;
                    }
                }
                "--wind-speed" => {
                    if i + 1 < args.len() {
                        wind_speed = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid wind speed"))?;
                        i += 1;
                    }
                }
                "--wind-direction" => {
                    if i + 1 < args.len() {
                        wind_direction = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid wind direction"))?;
                        i += 1;
                    }
                }
                "--wind-segment" => {
                    if i + 1 < args.len() {
                        wind_segment_strs.push(args[i + 1].to_string());
                        i += 1;
                    }
                }
                "--bc-segment" => {
                    if i + 1 < args.len() {
                        bc_segment_strs.push(args[i + 1].to_string());
                        i += 1;
                    }
                }
                "--temperature" => {
                    if i + 1 < args.len() {
                        temperature = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid temperature"))?;
                        i += 1;
                    }
                }
                "--pressure" => {
                    if i + 1 < args.len() {
                        pressure = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid pressure"))?;
                        i += 1;
                    }
                }
                "--humidity" => {
                    if i + 1 < args.len() {
                        humidity = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid humidity"))?;
                        i += 1;
                    }
                }
                "--altitude" => {
                    if i + 1 < args.len() {
                        altitude = args[i + 1]
                            .parse()
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
                        auto_zero = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid zero distance"))?,
                        );
                        i += 1;
                    }
                }
                "--sight-height" => {
                    if i + 1 < args.len() {
                        sight_height = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid sight height"))?;
                        i += 1;
                    }
                }
                "--muzzle-height" => {
                    if i + 1 < args.len() {
                        muzzle_height = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid muzzle height"))?;
                        i += 1;
                    }
                }
                "--target-height" => {
                    if i + 1 < args.len() {
                        target_height = args[i + 1]
                            .parse()
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
                        sample_interval = args[i + 1]
                            .parse()
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
                        twist_rate = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid twist rate"))?,
                        );
                        i += 1;
                    }
                }
                "--twist-right" => {
                    if i + 1 < args.len() {
                        twist_right = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid twist direction"))?;
                        i += 1;
                    }
                }
                "--latitude" => {
                    if i + 1 < args.len() {
                        latitude = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid latitude"))?,
                        );
                        i += 1;
                    }
                }
                "--shot-direction" => {
                    if i + 1 < args.len() {
                        shot_direction = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid shot-direction"))?,
                        );
                        i += 1;
                    }
                }
                "--shooting-angle" => {
                    if i + 1 < args.len() {
                        shooting_angle = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid shooting angle"))?;
                        i += 1;
                    }
                }
                "--cant" | "--cant-angle" => {
                    if i + 1 < args.len() {
                        cant_angle_deg = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid cant angle"))?;
                        i += 1;
                    }
                }
                "--powder-temp-sensitivity" => {
                    if i + 1 < args.len() {
                        powder_temp_sensitivity = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid powder temp sensitivity"))?;
                        i += 1;
                    }
                }
                "--powder-temp" => {
                    if i + 1 < args.len() {
                        powder_temp = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid powder temp"))?;
                        powder_temp_provided = true;
                        i += 1;
                    }
                }
                "--powder-temp-curve" => {
                    if i + 1 < args.len() {
                        powder_temp_curve_str = Some(args[i + 1].to_string());
                        i += 1;
                    }
                }
                "--zero-velocity" => {
                    if i + 1 < args.len() {
                        zero_velocity = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid zero velocity"))?,
                        );
                        i += 1;
                    }
                }
                "--zero-temperature" => {
                    if i + 1 < args.len() {
                        zero_temperature = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid zero temperature"))?,
                        );
                        i += 1;
                    }
                }
                "--zero-pressure" => {
                    if i + 1 < args.len() {
                        zero_pressure = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid zero pressure"))?,
                        );
                        i += 1;
                    }
                }
                "--zero-humidity" => {
                    if i + 1 < args.len() {
                        zero_humidity = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid zero humidity"))?,
                        );
                        i += 1;
                    }
                }
                "--zero-altitude" => {
                    if i + 1 < args.len() {
                        zero_altitude = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid zero altitude"))?,
                        );
                        i += 1;
                    }
                }
                "--zero-powder-temp" => {
                    if i + 1 < args.len() {
                        zero_powder_temp = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid zero powder temp"))?,
                        );
                        i += 1;
                    }
                }
                // --units/-u (+ its value) is consumed globally in run_command, which
                // pre-scans it to set the unit system before dispatch. Skip it here so
                // it isn't rejected as an unknown flag (this is what blocked metric input).
                "--units" | "-u" => {
                    i += 1;
                }
                // Reject unrecognized flags instead of silently ignoring them, so a
                // typo or a flag that isn't wired into this WASM surface is caught
                // immediately rather than looking like a no-op. (The native CLI's clap
                // parser already does this; the hand-rolled WASM parser did not.)
                other if other.starts_with('-') => {
                    return Err(JsValue::from_str(&format!("Unknown flag: {}", other)));
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
        // MBA-1135: mass-based length estimate (mirrors CLI/FFI), replacing the mass-blind
        // 4.5-caliber heuristic. WASM otherwise left it at the struct default regardless of the
        // supplied caliber/weight, skewing the Miller Sg / Litz spin drift / Magnus.
        inputs.bullet_length =
            crate::stability::estimate_bullet_length_m(inputs.bullet_diameter, inputs.bullet_mass);
        if inputs.bullet_length <= 0.0 {
            inputs.bullet_length = inputs.bullet_diameter * 4.5;
        }

        inputs.bc_value = bc;
        inputs.bc_type = DragModel::from_str(drag_model)
            .ok_or_else(|| JsValue::from_str("Invalid drag model"))?;
        inputs.muzzle_angle = angle * std::f64::consts::PI / 180.0; // degrees to radians
        inputs.shooting_angle = shooting_angle * std::f64::consts::PI / 180.0;
        inputs.cant_angle = cant_angle_deg * std::f64::consts::PI / 180.0;
        inputs.ground_threshold = 0.0;

        // Set advanced physics flags. enable_advanced_effects remains the umbrella
        // flag, but Magnus and Coriolis are now gated independently so enabling one
        // does not silently enable the other.
        if enable_magnus || enable_coriolis {
            inputs.enable_advanced_effects = true;
        }
        inputs.enable_magnus = enable_magnus;
        inputs.enable_coriolis = enable_coriolis;
        // Set integration method: Euler < RK4 fixed < RK45 adaptive (default)
        if use_euler {
            inputs.use_rk4 = false;
            inputs.use_adaptive_rk45 = false;
        } else if use_rk4_fixed {
            inputs.use_rk4 = true;
            inputs.use_adaptive_rk45 = false; // Fixed-step RK4
        } else {
            inputs.use_rk4 = true;
            inputs.use_adaptive_rk45 = true; // Default: adaptive RK45
        }
        inputs.use_enhanced_spin_drift = enable_spin_drift;
        inputs.enable_wind_shear = enable_wind_shear;
        inputs.enable_trajectory_sampling = sample_trajectory;
        inputs.sample_interval = sample_interval;
        inputs.enable_pitch_damping = enable_pitch_damping;
        inputs.enable_precession_nutation = enable_precession;
        inputs.use_bc_segments = use_bc_segments;
        inputs.use_powder_sensitivity = use_powder_sensitivity;

        // Velocity-keyed BC segments, in priority order:
        //   1. manual --bc-segment "VMIN:VMAX:BC" pairs (explicit user input)
        //   2. a BC5D table loaded via loadBc5dTable() + --use-bc-segments
        // Velocities are in the command's display units (fps imperial, m/s metric);
        // the solver compares against fps, so convert.
        let vel_to_fps = match units {
            UnitSystem::Imperial => 1.0,
            UnitSystem::Metric => 3.280_839_895,
        };
        let mut manual_bc_segments: Vec<crate::BCSegmentData> = Vec::new();
        for s in &bc_segment_strs {
            let parts: Vec<&str> = s.split(':').collect();
            if parts.len() != 3 {
                return Err(JsValue::from_str(&format!(
                    "--bc-segment expects VMIN:VMAX:BC (e.g. 1500:1800:0.243), got '{}'",
                    s
                )));
            }
            let vmin: f64 = parts[0].trim().parse().map_err(|_| {
                JsValue::from_str(&format!("--bc-segment: invalid VMIN in '{}'", s))
            })?;
            let vmax: f64 = parts[1].trim().parse().map_err(|_| {
                JsValue::from_str(&format!("--bc-segment: invalid VMAX in '{}'", s))
            })?;
            let bcv: f64 = parts[2].trim().parse().map_err(|_| {
                JsValue::from_str(&format!("--bc-segment: invalid BC in '{}'", s))
            })?;
            if !(vmin < vmax) {
                return Err(JsValue::from_str(&format!(
                    "--bc-segment: VMIN must be < VMAX in '{}'",
                    s
                )));
            }
            if bcv <= 0.0 {
                return Err(JsValue::from_str(&format!(
                    "--bc-segment: BC must be > 0 in '{}'",
                    s
                )));
            }
            manual_bc_segments.push(crate::BCSegmentData {
                velocity_min: vmin * vel_to_fps,
                velocity_max: vmax * vel_to_fps,
                bc_value: bcv,
            });
        }

        if !manual_bc_segments.is_empty() {
            // Manual segments win; imply --use-bc-segments so they're applied.
            inputs.bc_segments_data = Some(manual_bc_segments);
            inputs.use_bc_segments = true;
        } else if use_bc_segments {
            // BC5D offline correction: synthesize velocity-dependent BC segments
            // from a loaded table. The table's native units are grains + fps.
            if let Some(table) = self.bc5d_table.borrow().as_ref() {
                let (weight_grains, muzzle_fps) = match units {
                    UnitSystem::Imperial => (mass, velocity),
                    UnitSystem::Metric => (mass * 15.4323584, velocity * 3.280839895),
                };
                if let Some(schedule) =
                    table.generate_segment_schedule(bc, drag_model, weight_grains, muzzle_fps)
                {
                    inputs.bc_segments_data = Some(schedule.segments);
                    inputs.bc_value = schedule.fallback_bc;
                }
            }
        }

        // Set additional parameters
        let explicit_twist_inches = twist_rate.map(|rate| match units {
            UnitSystem::Imperial => rate,
            UnitSystem::Metric => rate / 25.4,
        });
        inputs.twist_rate = crate::stability::resolve_twist_inches(
            explicit_twist_inches,
            inputs.bullet_diameter,
            inputs.bullet_mass,
            inputs.muzzle_velocity,
        );
        inputs.is_twist_right = twist_right;
        if let Some(lat) = latitude {
            inputs.latitude = Some(lat);
        }
        inputs.shot_azimuth = shot_direction.map(|d| d.to_radians()).unwrap_or(0.0);
        inputs.powder_temp_sensitivity = match units {
            UnitSystem::Imperial => powder_temp_sensitivity * 0.3048 / (5.0 / 9.0),
            UnitSystem::Metric => powder_temp_sensitivity,
        };
        inputs.powder_temp = match units {
            UnitSystem::Imperial => (powder_temp - 32.0) * 5.0 / 9.0,
            UnitSystem::Metric => powder_temp,
        };
        // When --powder-temp was explicitly given, it becomes the POWDER temperature the
        // curve is interpolated at (decoupled from --temperature / air density). When not
        // given, powder_curve_temp_c stays None so the curve falls back to the air temp.
        inputs.powder_curve_temp_c = if powder_temp_provided {
            Some(inputs.powder_temp)
        } else {
            None
        };
        // Parse the optional powder-temperature -> velocity curve into SI points.
        if let Some(curve_str) = &powder_temp_curve_str {
            let mut pts: Vec<(f64, f64)> = Vec::new();
            for raw in curve_str.split(',') {
                let part = raw.trim();
                if part.is_empty() {
                    continue;
                }
                let (t_str, v_str) = part.split_once(':').ok_or_else(|| {
                    JsValue::from_str("--powder-temp-curve point must be TEMP:VELOCITY")
                })?;
                let t: f64 = t_str
                    .trim()
                    .parse()
                    .map_err(|_| JsValue::from_str("Invalid temperature in --powder-temp-curve"))?;
                let v: f64 = v_str
                    .trim()
                    .parse()
                    .map_err(|_| JsValue::from_str("Invalid velocity in --powder-temp-curve"))?;
                if !(v > 0.0) {
                    return Err(JsValue::from_str(
                        "--powder-temp-curve velocity must be positive",
                    ));
                }
                let (t_c, v_ms) = match units {
                    UnitSystem::Imperial => ((t - 32.0) * 5.0 / 9.0, v * 0.3048),
                    UnitSystem::Metric => (t, v),
                };
                pts.push((t_c, v_ms));
            }
            if pts.len() < 2 {
                return Err(JsValue::from_str(
                    "--powder-temp-curve needs at least 2 TEMP:VELOCITY points",
                ));
            }
            pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            inputs.powder_temp_curve = Some(pts);
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
        // WindConditions.direction is RADIANS (0=North, PI/2=East); --wind-direction is degrees.
        // Convert (matches native CLI); previously a 90-degree crosswind was fed as 90 radians.
        wind.direction = wind_direction.to_radians();

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
        inputs.temperature = atmosphere.temperature;
        inputs.pressure = atmosphere.pressure;
        inputs.humidity = (humidity / 100.0).clamp(0.0, 1.0);
        inputs.altitude = atmosphere.altitude;

        // Handle auto-zero if specified
        let mut zero_info = String::new();
        if let Some(zero_distance) = auto_zero {
            let zero_distance_m = match units {
                UnitSystem::Imperial => zero_distance * 0.9144, // yards to meters
                UnitSystem::Metric => zero_distance,
            };

            // Build the condition set the zero ANGLE is solved under. It starts from the
            // shot-day inputs/atmosphere and is overridden only by whichever --zero-*
            // flags the caller supplied, so a rifle zeroed on a different day (different
            // air density and/or muzzle velocity) is modeled correctly while the
            // trajectory below still runs under the current shot-day conditions.
            let mut zero_inputs = inputs.clone();
            let mut zero_atmosphere = atmosphere.clone();
            let mut zero_day_overridden = false;
            if let Some(zv) = zero_velocity {
                zero_inputs.muzzle_velocity = match units {
                    UnitSystem::Imperial => zv * 0.3048, // fps to m/s
                    UnitSystem::Metric => zv,
                };
                // An explicit zero-day velocity takes precedence: disable both velocity
                // adjustment models for the zero solve so neither changes the supplied value.
                // (zero_inputs is a clone of inputs, which may carry the shot-day models.)
                zero_inputs.powder_temp_curve = None;
                zero_inputs.use_powder_sensitivity = false;
                zero_day_overridden = true;
            }
            if let Some(zt) = zero_temperature {
                let t_c = match units {
                    UnitSystem::Imperial => (zt - 32.0) * 5.0 / 9.0, // F to C
                    UnitSystem::Metric => zt,
                };
                zero_atmosphere.temperature = t_c;
                zero_inputs.temperature = t_c;
                zero_day_overridden = true;
            }
            if let Some(zp) = zero_pressure {
                let p_hpa = match units {
                    UnitSystem::Imperial => zp * 33.863886666667, // inHg to hPa
                    UnitSystem::Metric => zp,
                };
                zero_atmosphere.pressure = p_hpa;
                zero_inputs.pressure = p_hpa;
                zero_day_overridden = true;
            }
            if let Some(zh) = zero_humidity {
                zero_atmosphere.humidity = zh;
                zero_inputs.humidity = (zh / 100.0).clamp(0.0, 1.0);
                zero_day_overridden = true;
            }
            if let Some(za) = zero_altitude {
                let a_m = match units {
                    UnitSystem::Imperial => za * 0.3048, // feet to meters
                    UnitSystem::Metric => za,
                };
                zero_atmosphere.altitude = a_m;
                zero_inputs.altitude = a_m;
                zero_day_overridden = true;
            }
            // An explicit zero-day powder temperature wins. Otherwise an explicit zero-day air
            // temperature retains the established "powder follows zero-day air" behavior. With
            // neither override, keep the cloned shot-day powder temperature so a no-override
            // zero solve reproduces the flight conditions exactly.
            if let Some(zpt) = zero_powder_temp {
                zero_inputs.powder_curve_temp_c = Some(match units {
                    UnitSystem::Imperial => (zpt - 32.0) * 5.0 / 9.0,
                    UnitSystem::Metric => zpt,
                });
                zero_day_overridden = true;
            } else if zero_temperature.is_some() {
                zero_inputs.powder_curve_temp_c = None;
            }

            match calculate_zero_angle_with_conditions(
                zero_inputs.clone(),
                zero_distance_m,
                zero_inputs.muzzle_height + zero_inputs.sight_height, // Zero crosses the line of sight (matches CLI)
                wind.clone(),
                zero_atmosphere.clone(),
            ) {
                Ok(zero_angle) => {
                    inputs.muzzle_angle = zero_angle;
                    let moa_adjustment = zero_angle * 180.0 / std::f64::consts::PI * 60.0;
                    let mrad_adjustment = zero_angle * 1000.0;
                    zero_info = format!(
                        "Rifle zeroed at {} {}{} (Adjustment: {:.2} MOA / {:.2} mrad up)\n\n",
                        zero_distance,
                        if units == UnitSystem::Imperial {
                            "yards"
                        } else {
                            "meters"
                        },
                        if zero_day_overridden {
                            " under supplied zero-day conditions"
                        } else {
                            ""
                        },
                        moa_adjustment,
                        mrad_adjustment
                    );
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

        // Downrange-segmented wind overrides the scalar wind when present.
        if !wind_segment_strs.is_empty() {
            if enable_wind_shear {
                return Err(JsValue::from_str(
                    "--wind-segment cannot be combined with --enable-wind-shear \
                     (downrange segments + altitude shear is not yet a defined model)",
                ));
            }
            let imperial = matches!(units, UnitSystem::Imperial);
            let mut segments = Vec::with_capacity(wind_segment_strs.len());
            for s in &wind_segment_strs {
                let seg = crate::wind::parse_wind_segment_str(s, imperial)
                    .map_err(|err| JsValue::from_str(&err))?;
                segments.push(seg);
            }
            solver.set_wind_segments(segments);
        }

        match solver.solve() {
            Ok(result) => {
                let output = match output_format {
                    OutputFormat::Table => self.format_trajectory_table(
                        &result,
                        auto_zero,
                        units,
                        full,
                        inputs.sight_height,
                    ),
                    OutputFormat::Json => {
                        self.format_trajectory_json(&result, units, inputs.sight_height)
                    }
                    OutputFormat::Csv => {
                        self.format_trajectory_csv(&result, units, full, inputs.sight_height)
                    }
                };
                Ok(format!("{}{}", zero_info, output))
            }
            Err(e) => Ok(format!("Error: {}", e)),
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
        let mut target_distance = if units == UnitSystem::Imperial {
            100.0
        } else {
            100.0
        };
        let mut sight_height = if units == UnitSystem::Imperial {
            2.0
        } else {
            50.0
        };
        // Heights above GROUND (--muzzle-height / --target-height) do NOT change the zero ANGLE —
        // for a same-elevation target they cancel — so they are intentionally not parsed here
        // (silently ignored by the catch-all arm below). The SIGHT height IS honored: the zero
        // targets the line-of-sight height at the zero distance (see the calculate_zero call).
        let mut drag_model = "G1";

        // Parse arguments
        let mut i = 0;
        while i < args.len() {
            match args[i] {
                "-v" | "--velocity" => {
                    if i + 1 < args.len() {
                        velocity = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid velocity"))?;
                        i += 1;
                    }
                }
                "-b" | "--bc" => {
                    if i + 1 < args.len() {
                        bc = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid BC"))?;
                        i += 1;
                    }
                }
                "-m" | "--mass" => {
                    if i + 1 < args.len() {
                        mass = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid mass"))?;
                        i += 1;
                    }
                }
                "-d" | "--diameter" => {
                    if i + 1 < args.len() {
                        diameter = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid diameter"))?;
                        i += 1;
                    }
                }
                "--target-distance" => {
                    if i + 1 < args.len() {
                        target_distance = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid target distance"))?;
                        i += 1;
                    }
                }
                "--sight-height" => {
                    if i + 1 < args.len() {
                        sight_height = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid sight height"))?;
                        i += 1;
                    }
                }
                "--drag-model" => {
                    if i + 1 < args.len() {
                        drag_model = args[i + 1];
                        i += 1;
                    }
                }
                // --units/-u (+ its value) is consumed globally in run_command, which
                // pre-scans it to set the unit system before dispatch. Skip it here so
                // it isn't rejected as an unknown flag (this is what blocked metric input).
                "--units" | "-u" => {
                    i += 1;
                }
                // Reject unrecognized flags instead of silently ignoring them, so a
                // typo or a flag that isn't wired into this WASM surface is caught
                // immediately rather than looking like a no-op. (The native CLI's clap
                // parser already does this; the hand-rolled WASM parser did not.)
                other if other.starts_with('-') => {
                    return Err(JsValue::from_str(&format!("Unknown flag: {}", other)));
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
        // MBA-1135: mass-based length estimate (mirrors CLI/FFI); replaces the mass-blind
        // 4.5-caliber heuristic. WASM otherwise left it at the struct default.
        inputs.bullet_length =
            crate::stability::estimate_bullet_length_m(inputs.bullet_diameter, inputs.bullet_mass);
        if inputs.bullet_length <= 0.0 {
            inputs.bullet_length = inputs.bullet_diameter * 4.5;
        }

        inputs.bc_value = bc;
        inputs.bc_type = DragModel::from_str(drag_model)
            .ok_or_else(|| JsValue::from_str("Invalid drag model"))?;

        let target_distance_m = match units {
            UnitSystem::Imperial => target_distance * 0.9144,
            UnitSystem::Metric => target_distance,
        };

        // MBA-951: target the line-of-sight height at the zero distance (= sight_height), matching
        // the CLI convention in every zero call. Previously 0.0, which solved a BORE-line zero and
        // ignored sight height entirely (off by the sight-height angle — ~2 MOA at 100 yd).
        let los_height = inputs.sight_height;
        match calculate_zero_angle_with_conditions(
            inputs,
            target_distance_m,
            los_height,
            WindConditions::default(),
            AtmosphericConditions::default(),
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
                    if units == UnitSystem::Imperial {
                        "yards"
                    } else {
                        "meters"
                    },
                    zero_degrees,
                    moa_adjustment,
                    mrad_adjustment,
                    sight_height,
                    if units == UnitSystem::Imperial {
                        "inches"
                    } else {
                        "mm"
                    }
                ))
            }
            Err(e) => Ok(format!("Error calculating zero: {}", e)),
        }
    }

    fn handle_lead_command(&self, args: &[&str], units: UnitSystem) -> Result<String, JsValue> {
        // Default values
        let (default_velocity, default_mass, default_diameter) = match units {
            UnitSystem::Imperial => (2700.0, 168.0, 0.308),
            UnitSystem::Metric => (823.0, 10.9, 7.82),
        };

        let mut velocity = default_velocity;
        let mut bc = 0.475;
        let mut mass = default_mass;
        let mut diameter = default_diameter;
        let mut drag_model = "G1";
        let mut sight_height = if units == UnitSystem::Imperial {
            2.0
        } else {
            50.0
        };
        let mut target_speed: Option<f64> = None;
        let mut target_angle = 90.0;
        let mut range = 500.0;
        let mut adjustment_unit = "mil";

        // Parse arguments
        let mut i = 0;
        while i < args.len() {
            match args[i] {
                "-v" | "--velocity" => {
                    if i + 1 < args.len() {
                        velocity = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid velocity"))?;
                        i += 1;
                    }
                }
                "-b" | "--bc" => {
                    if i + 1 < args.len() {
                        bc = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid BC"))?;
                        i += 1;
                    }
                }
                "-m" | "--mass" => {
                    if i + 1 < args.len() {
                        mass = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid mass"))?;
                        i += 1;
                    }
                }
                "-d" | "--diameter" => {
                    if i + 1 < args.len() {
                        diameter = args[i + 1]
                            .parse()
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
                "--sight-height" => {
                    if i + 1 < args.len() {
                        sight_height = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid sight height"))?;
                        i += 1;
                    }
                }
                "--target-speed" => {
                    if i + 1 < args.len() {
                        target_speed = Some(
                            args[i + 1]
                                .parse()
                                .map_err(|_| JsValue::from_str("Invalid target speed"))?,
                        );
                        i += 1;
                    }
                }
                "--target-angle" => {
                    if i + 1 < args.len() {
                        target_angle = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid target angle"))?;
                        i += 1;
                    }
                }
                "--range" => {
                    if i + 1 < args.len() {
                        range = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid range"))?;
                        i += 1;
                    }
                }
                "--adjustment-unit" => {
                    if i + 1 < args.len() {
                        adjustment_unit = args[i + 1];
                        i += 1;
                    }
                }
                // --units/-u (+ its value) is consumed globally in run_command, which
                // pre-scans it to set the unit system before dispatch. Skip it here so
                // it isn't rejected as an unknown flag (this is what blocked metric input).
                "--units" | "-u" => {
                    i += 1;
                }
                // Reject unrecognized flags instead of silently ignoring them, so a
                // typo or a flag that isn't wired into this WASM surface is caught
                // immediately rather than looking like a no-op. (The native CLI's clap
                // parser already does this; the hand-rolled WASM parser did not.)
                other if other.starts_with('-') => {
                    return Err(JsValue::from_str(&format!("Unknown flag: {}", other)));
                }
                _ => {}
            }
            i += 1;
        }

        let target_speed = target_speed
            .ok_or_else(|| JsValue::from_str("--target-speed is required"))?;

        let adjustment_unit_lower = adjustment_unit.to_lowercase();
        if adjustment_unit_lower != "mil" && adjustment_unit_lower != "moa" {
            return Err(JsValue::from_str(&format!(
                "Invalid --adjustment-unit '{}' (expected mil or moa)",
                adjustment_unit
            )));
        }

        // Build inputs (mirrors handle_zero_command's unit conversions)
        let mut inputs = InternalBallisticInputs::default();
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
        // MBA-1135: mass-based length estimate (mirrors CLI/FFI); replaces the mass-blind
        // 4.5-caliber heuristic. WASM otherwise left it at the struct default.
        inputs.bullet_length =
            crate::stability::estimate_bullet_length_m(inputs.bullet_diameter, inputs.bullet_mass);
        if inputs.bullet_length <= 0.0 {
            inputs.bullet_length = inputs.bullet_diameter * 4.5;
        }

        inputs.bc_value = bc;
        inputs.bc_type = DragModel::from_str(drag_model)
            .ok_or_else(|| JsValue::from_str("Invalid drag model"))?;

        // --target-speed uses the same mph (imperial) / m/s (metric) convention as
        // --wind-speed in handle_trajectory_command.
        let target_speed_mps = match units {
            UnitSystem::Imperial => target_speed * 0.44704, // mph to m/s
            UnitSystem::Metric => target_speed,
        };
        let range_m = match units {
            UnitSystem::Imperial => range * 0.9144, // yards to meters
            UnitSystem::Metric => range,
        };

        match calculate_lead(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
            target_speed_mps,
            target_angle,
            range_m,
        ) {
            Ok(sol) => {
                let (dist_unit, speed_unit) = match units {
                    UnitSystem::Imperial => ("yd", "mph"),
                    UnitSystem::Metric => ("m", "m/s"),
                };
                let lead_disp = match units {
                    UnitSystem::Imperial => sol.lead_m * 1.09361, // m to yards
                    UnitSystem::Metric => sol.lead_m,
                };
                let intercept_disp = match units {
                    UnitSystem::Imperial => sol.corrected_range_m * 1.09361,
                    UnitSystem::Metric => sol.corrected_range_m,
                };
                // Requested --adjustment-unit is listed first; both MIL and MOA are always shown.
                let lead_adj_line = if adjustment_unit_lower == "moa" {
                    format!(
                        "{:.2} MOA / {:.2} MIL",
                        sol.lead_moa, sol.lead_mil
                    )
                } else {
                    format!(
                        "{:.2} MIL / {:.2} MOA",
                        sol.lead_mil, sol.lead_moa
                    )
                };

                Ok(format!(
                    "Moving-Target Lead\n\
                     ===================\n\
                     Target: {:.1} {} at {:.0}\u{b0} \
                     (0=away, 90=left-to-right, 180=toward, 270=right-to-left;\n\
                     positive lead = hold in direction of travel)\n\n\
                     Range: {:.0} {}\n\
                     Time of Flight: {:.3} s\n\
                     Lead: {:.2} {} ({})\n\
                     Intercept Range: {:.1} {}\n\
                     Iterations: {}\n",
                    target_speed,
                    speed_unit,
                    target_angle,
                    range,
                    dist_unit,
                    sol.time_of_flight_s,
                    lead_disp,
                    dist_unit,
                    lead_adj_line,
                    intercept_disp,
                    dist_unit,
                    sol.iterations
                ))
            }
            Err(e) => Err(JsValue::from_str(&format!("Error calculating lead: {}", e))),
        }
    }

    fn handle_monte_carlo_command(
        &self,
        args: &[&str],
        units: UnitSystem,
    ) -> Result<String, JsValue> {
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
        let mut velocity_std = if units == UnitSystem::Imperial {
            10.0
        } else {
            3.0
        };
        let mut angle_std = 0.1;
        let mut bc_std = 0.01;
        let mut wind_speed_std = if units == UnitSystem::Imperial {
            2.0
        } else {
            0.5
        };
        let mut wind_direction_std = 0.0;
        let mut drag_model = "G1";

        // Parse arguments
        let mut i = 0;
        while i < args.len() {
            match args[i] {
                "-v" | "--velocity" => {
                    if i + 1 < args.len() {
                        velocity = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid velocity"))?;
                        i += 1;
                    }
                }
                "-a" | "--angle" => {
                    if i + 1 < args.len() {
                        angle = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid angle"))?;
                        i += 1;
                    }
                }
                "-b" | "--bc" => {
                    if i + 1 < args.len() {
                        bc = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid BC"))?;
                        i += 1;
                    }
                }
                "-m" | "--mass" => {
                    if i + 1 < args.len() {
                        mass = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid mass"))?;
                        i += 1;
                    }
                }
                "-d" | "--diameter" => {
                    if i + 1 < args.len() {
                        diameter = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid diameter"))?;
                        i += 1;
                    }
                }
                "-n" | "--num-sims" => {
                    if i + 1 < args.len() {
                        num_sims = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid number of simulations"))?;
                        i += 1;
                    }
                }
                "--velocity-std" => {
                    if i + 1 < args.len() {
                        velocity_std = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid velocity std"))?;
                        i += 1;
                    }
                }
                "--angle-std" => {
                    if i + 1 < args.len() {
                        angle_std = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid angle std"))?;
                        i += 1;
                    }
                }
                "--bc-std" => {
                    if i + 1 < args.len() {
                        bc_std = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid BC std"))?;
                        i += 1;
                    }
                }
                "--wind-speed-std" => {
                    if i + 1 < args.len() {
                        wind_speed_std = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid wind speed std"))?;
                        i += 1;
                    }
                }
                "--wind-direction-std" | "--wind-dir-std" => {
                    if i + 1 < args.len() {
                        wind_direction_std = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid wind direction std"))?;
                        i += 1;
                    }
                }
                "--drag-model" => {
                    if i + 1 < args.len() {
                        drag_model = args[i + 1];
                        i += 1;
                    }
                }
                // --units/-u (+ its value) is consumed globally in run_command, which
                // pre-scans it to set the unit system before dispatch. Skip it here so
                // it isn't rejected as an unknown flag (this is what blocked metric input).
                "--units" | "-u" => {
                    i += 1;
                }
                // Reject unrecognized flags instead of silently ignoring them, so a
                // typo or a flag that isn't wired into this WASM surface is caught
                // immediately rather than looking like a no-op. (The native CLI's clap
                // parser already does this; the hand-rolled WASM parser did not.)
                other if other.starts_with('-') => {
                    return Err(JsValue::from_str(&format!("Unknown flag: {}", other)));
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

        // MBA-1135: mass-based length estimate (mirrors CLI/FFI); replaces the mass-blind
        // 4.5-caliber heuristic. WASM otherwise left it at the struct default.
        inputs.bullet_length =
            crate::stability::estimate_bullet_length_m(inputs.bullet_diameter, inputs.bullet_mass);
        if inputs.bullet_length <= 0.0 {
            inputs.bullet_length = inputs.bullet_diameter * 4.5;
        }

        inputs.bc_value = bc;
        // Honor --drag-model (mirrors the trajectory/zero handlers); previously the Monte
        // Carlo path silently always used the G1 default even when G7 was intended.
        inputs.bc_type = DragModel::from_str(drag_model)
            .ok_or_else(|| JsValue::from_str("Invalid drag model (expected G1 or G7)"))?;
        inputs.muzzle_angle = angle * std::f64::consts::PI / 180.0;
        inputs.muzzle_height = 1.5;
        inputs.ground_threshold = 0.0;

        // Create Monte Carlo parameters
        let params = MonteCarloParams {
            num_simulations: num_sims,
            velocity_std_dev: velocity_std
                * (if units == UnitSystem::Imperial {
                    0.3048
                } else {
                    1.0
                }),
            angle_std_dev: angle_std * std::f64::consts::PI / 180.0,
            bc_std_dev: bc_std,
            wind_speed_std_dev: wind_speed_std
                * (if units == UnitSystem::Imperial {
                    0.44704
                } else {
                    1.0
                }),
            base_wind_speed: 0.0,
            base_wind_direction: 0.0,
            azimuth_std_dev: 0.001,
            target_distance: None,
        };

        match run_monte_carlo_with_direction_std_dev(
            inputs,
            params,
            wind_direction_std * std::f64::consts::PI / 180.0,
        ) {
            Ok(results) => {
                // Calculate statistics
                let mean_range: f64 =
                    results.ranges.iter().sum::<f64>() / results.ranges.len() as f64;
                let mean_velocity: f64 = results.impact_velocities.iter().sum::<f64>()
                    / results.impact_velocities.len() as f64;

                let range_std = {
                    let variance: f64 = results
                        .ranges
                        .iter()
                        .map(|r| (r - mean_range).powi(2))
                        .sum::<f64>()
                        / results.ranges.len() as f64;
                    variance.sqrt()
                };

                let velocity_std = {
                    let variance: f64 = results
                        .impact_velocities
                        .iter()
                        .map(|v| (v - mean_velocity).powi(2))
                        .sum::<f64>()
                        / results.impact_velocities.len() as f64;
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
                    display_mean_range,
                    range_unit,
                    display_range_std,
                    range_unit,
                    results
                        .ranges
                        .iter()
                        .min_by(|a, b| a.partial_cmp(b).unwrap())
                        .unwrap()
                        * (if units == UnitSystem::Imperial {
                            1.09361
                        } else {
                            1.0
                        }),
                    range_unit,
                    results
                        .ranges
                        .iter()
                        .max_by(|a, b| a.partial_cmp(b).unwrap())
                        .unwrap()
                        * (if units == UnitSystem::Imperial {
                            1.09361
                        } else {
                            1.0
                        }),
                    range_unit,
                    display_mean_velocity,
                    velocity_unit,
                    display_velocity_std,
                    velocity_unit,
                    results
                        .impact_velocities
                        .iter()
                        .min_by(|a, b| a.partial_cmp(b).unwrap())
                        .unwrap()
                        * (if units == UnitSystem::Imperial {
                            3.28084
                        } else {
                            1.0
                        }),
                    velocity_unit,
                    results
                        .impact_velocities
                        .iter()
                        .max_by(|a, b| a.partial_cmp(b).unwrap())
                        .unwrap()
                        * (if units == UnitSystem::Imperial {
                            3.28084
                        } else {
                            1.0
                        }),
                    velocity_unit,
                ))
            }
            Err(e) => Ok(format!("Error running Monte Carlo simulation: {}", e)),
        }
    }

    fn handle_estimate_bc_command(
        &self,
        args: &[&str],
        units: UnitSystem,
    ) -> Result<String, JsValue> {
        let (default_velocity, default_mass, default_diameter) = match units {
            UnitSystem::Imperial => (2700.0, 168.0, 0.308),
            UnitSystem::Metric => (823.0, 10.9, 7.82),
        };

        let mut velocity = default_velocity;
        let mut mass = default_mass;
        let mut diameter = default_diameter;
        let mut data_points: Vec<(f64, f64)> = Vec::new();
        let mut vel_points: Vec<(f64, f64)> = Vec::new();
        let mut drag_model_str = String::from("both");
        let mut zero_range: Option<f64> = None;
        let mut sight_height: Option<f64> = None;
        let mut temperature: Option<f64> = None;
        let mut pressure: Option<f64> = None;
        let mut humidity: f64 = 50.0;
        let mut altitude: f64 = 0.0;

        // Parse "d,v;d,v;..." pairs, tolerating surrounding quotes/whitespace.
        fn parse_pairs(raw: &str) -> Result<Vec<(f64, f64)>, JsValue> {
            let s = raw.trim().trim_matches('\'').trim_matches('"');
            let mut out = Vec::new();
            for pair in s.split(';') {
                let pair = pair.trim();
                if pair.is_empty() {
                    continue;
                }
                let parts: Vec<&str> = pair.split(',').collect();
                if parts.len() != 2 {
                    return Err(JsValue::from_str(&format!(
                        "Malformed data pair '{}': expected \"distance,value\"",
                        pair
                    )));
                }
                let d: f64 = parts[0]
                    .trim()
                    .parse()
                    .map_err(|_| JsValue::from_str(&format!("Invalid distance '{}'", parts[0].trim())))?;
                let v: f64 = parts[1]
                    .trim()
                    .parse()
                    .map_err(|_| JsValue::from_str(&format!("Invalid value '{}'", parts[1].trim())))?;
                out.push((d, v));
            }
            Ok(out)
        }

        // Parse arguments
        let mut i = 0;
        while i < args.len() {
            match args[i] {
                "-v" | "--velocity" => {
                    if i + 1 < args.len() {
                        velocity = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid velocity"))?;
                        i += 1;
                    }
                }
                "-m" | "--mass" => {
                    if i + 1 < args.len() {
                        mass = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid mass"))?;
                        i += 1;
                    }
                }
                "-d" | "--diameter" => {
                    if i + 1 < args.len() {
                        diameter = args[i + 1]
                            .parse()
                            .map_err(|_| JsValue::from_str("Invalid diameter"))?;
                        i += 1;
                    }
                }
                "--data" => {
                    // Drop data: distance,drop pairs.
                    if i + 1 < args.len() {
                        data_points = parse_pairs(args[i + 1])?;
                        i += 1;
                    }
                }
                "--velocity-data" => {
                    // Velocity-retention data: distance,velocity pairs.
                    if i + 1 < args.len() {
                        vel_points = parse_pairs(args[i + 1])?;
                        i += 1;
                    }
                }
                "--drag-model" => {
                    if i + 1 < args.len() {
                        drag_model_str = args[i + 1].to_lowercase();
                        i += 1;
                    }
                }
                "--zero-range" => {
                    if i + 1 < args.len() {
                        zero_range = Some(args[i + 1].parse().map_err(|_| JsValue::from_str("Invalid zero-range"))?);
                        i += 1;
                    }
                }
                "--sight-height" => {
                    if i + 1 < args.len() {
                        sight_height = Some(args[i + 1].parse().map_err(|_| JsValue::from_str("Invalid sight-height"))?);
                        i += 1;
                    }
                }
                "--temperature" => {
                    if i + 1 < args.len() {
                        temperature = Some(args[i + 1].parse().map_err(|_| JsValue::from_str("Invalid temperature"))?);
                        i += 1;
                    }
                }
                "--pressure" => {
                    if i + 1 < args.len() {
                        pressure = Some(args[i + 1].parse().map_err(|_| JsValue::from_str("Invalid pressure"))?);
                        i += 1;
                    }
                }
                "--humidity" => {
                    if i + 1 < args.len() {
                        humidity = args[i + 1].parse().map_err(|_| JsValue::from_str("Invalid humidity"))?;
                        i += 1;
                    }
                }
                "--altitude" => {
                    if i + 1 < args.len() {
                        altitude = args[i + 1].parse().map_err(|_| JsValue::from_str("Invalid altitude"))?;
                        i += 1;
                    }
                }
                // --units/-u (+ its value) is consumed globally in run_command, which
                // pre-scans it to set the unit system before dispatch. Skip it here so
                // it isn't rejected as an unknown flag (this is what blocked metric input).
                "--units" | "-u" => {
                    i += 1;
                }
                // Reject unrecognized flags instead of silently ignoring them, so a
                // typo or a flag that isn't wired into this WASM surface is caught
                // immediately rather than looking like a no-op. (The native CLI's clap
                // parser already does this; the hand-rolled WASM parser did not.)
                other if other.starts_with('-') => {
                    return Err(JsValue::from_str(&format!("Unknown flag: {}", other)));
                }
                _ => {}
            }
            i += 1;
        }

        if data_points.is_empty() && vel_points.is_empty() {
            return Ok("Error: No data provided. Use --data \"dist,drop;...\" and/or \
                 --velocity-data \"dist,vel;...\".\nExample: --data \"300,29.0;500,89.9;700,204.6\""
                .to_string());
        }

        // Select drag models to estimate.
        let models: Vec<DragModel> = match drag_model_str.as_str() {
            "g1" => vec![DragModel::G1],
            "g7" => vec![DragModel::G7],
            "both" | "all" | "g1,g7" | "g1g7" => vec![DragModel::G1, DragModel::G7],
            other => {
                return Err(JsValue::from_str(&format!(
                    "Unknown --drag-model '{}'; use g1, g7, or both.",
                    other
                )))
            }
        };

        // Convert scalar inputs to metric.
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

        // Atmosphere the data was measured at (defaults = ICAO standard). BC only means
        // something relative to air density, so this must match the dope card's conditions.
        let atmosphere = AtmosphericConditions {
            temperature: temperature
                .map(|t| match units {
                    UnitSystem::Imperial => (t - 32.0) * 5.0 / 9.0,
                    UnitSystem::Metric => t,
                })
                .unwrap_or(15.0),
            pressure: pressure
                .map(|p| match units {
                    UnitSystem::Imperial => p * 33.8639,
                    UnitSystem::Metric => p,
                })
                .unwrap_or(1013.25),
            humidity,
            altitude: match units {
                UnitSystem::Imperial => altitude * 0.3048,
                UnitSystem::Metric => altitude,
            },
        };
        let zero_m = zero_range.map(|z| match units {
            UnitSystem::Imperial => z * 0.9144,
            UnitSystem::Metric => z,
        });
        let sight_m = sight_height
            .map(|s| match units {
                UnitSystem::Imperial => s * 0.0254,
                UnitSystem::Metric => s / 1000.0,
            })
            .unwrap_or(0.05);

        // Convert both series to metric (drop -> m, velocity -> m/s; metric drop input is mm).
        let drop_metric: Vec<(f64, f64)> = data_points
            .iter()
            .map(|(dist, drop)| match units {
                UnitSystem::Imperial => (*dist * 0.9144, *drop * 0.0254),
                UnitSystem::Metric => (*dist, *drop * 0.001),
            })
            .collect();
        let vel_metric: Vec<(f64, f64)> = vel_points
            .iter()
            .map(|(dist, v)| match units {
                UnitSystem::Imperial => (*dist * 0.9144, *v * 0.3048),
                UnitSystem::Metric => (*dist, *v),
            })
            .collect();

        // Data bases in a stable order (drop first, then velocity).
        let mut bases: Vec<(BcFitMode, &[(f64, f64)])> = Vec::new();
        if !drop_metric.is_empty() {
            bases.push((BcFitMode::Drop, &drop_metric));
        }
        if !vel_metric.is_empty() {
            bases.push((BcFitMode::Velocity, &vel_metric));
        }

        let (vu, mu, du) = match units {
            UnitSystem::Imperial => ("fps", "grains", "inches"),
            UnitSystem::Metric => ("m/s", "grams", "mm"),
        };
        let mut lines = vec![
            "BC Estimation Results".to_string(),
            "=====================".to_string(),
            format!(
                "Inputs: v={} {}, m={} {}, d={} {}",
                velocity, vu, mass, mu, diameter, du
            ),
            String::new(),
            format!(
                "  {:<6} {:<20} {:>12}   {}",
                "Model", "Fit basis", "Estimated BC", "Fit RMS"
            ),
        ];
        // Guard the common mistake: a zeroed dope card (a point with ~0 drop) fed without
        // --zero-range, which makes the drop fit bore-referenced and returns a wrong BC.
        if zero_m.is_none()
            && drop_metric.iter().any(|(_, dr)| dr.abs() < 0.05)
            && drop_metric.iter().any(|(_, dr)| dr.abs() > 0.25)
        {
            let zd = drop_metric
                .iter()
                .min_by(|a, b| a.1.abs().partial_cmp(&b.1.abs()).unwrap())
                .map(|(d, _)| match units {
                    UnitSystem::Imperial => d / 0.9144,
                    UnitSystem::Metric => *d,
                })
                .unwrap_or(0.0);
            let du2 = if units == UnitSystem::Imperial { "yd" } else { "m" };
            lines.insert(0, format!(
                "⚠ Data looks zeroed near {zd:.0} {du2} but --zero-range not given; drop is being"
            ));
            lines.insert(1, format!(
                "  treated as bore-referenced. For a dope card, pass --zero-range {zd:.0}."
            ));
            lines.insert(2, String::new());
        }

        let mut any_unreliable = false;
        for &model in &models {
            for &(mode, pts) in &bases {
                // Zero range only shapes a drop fit; a velocity fit is frame-independent.
                let zr = match mode {
                    BcFitMode::Drop => zero_m,
                    BcFitMode::Velocity => None,
                };
                let est = estimate_bc_fit(
                    velocity_mps, mass_kg, diameter_m, pts, model, mode,
                    atmosphere.clone(), zr, sight_m,
                )
                .map_err(|e| JsValue::from_str(&format!("Error estimating BC: {}", e)))?;
                if est.at_bound {
                    any_unreliable = true;
                }
                let (rms_user, unit) = match mode {
                    BcFitMode::Drop => match units {
                        UnitSystem::Imperial => (est.rms_error / 0.0254, "in"),
                        UnitSystem::Metric => (est.rms_error * 1000.0, "mm"),
                    },
                    BcFitMode::Velocity => match units {
                        UnitSystem::Imperial => (est.rms_error / 0.3048, "fps"),
                        UnitSystem::Metric => (est.rms_error, "m/s"),
                    },
                };
                let model_name = match model {
                    DragModel::G7 => "G7",
                    DragModel::G1 => "G1",
                    _ => "G?",
                };
                let basis = match mode {
                    BcFitMode::Drop => "drop",
                    BcFitMode::Velocity => "velocity",
                };
                lines.push(format!(
                    "  {:<6} {:<20} {:>12.3}   {:>6.2} {:<4}{}",
                    model_name,
                    format!("{} ({} pts)", basis, pts.len()),
                    est.bc,
                    rms_user,
                    unit,
                    if est.at_bound { " ⚠ UNRELIABLE (hit BC limit)" } else { "" }
                ));
            }
        }
        if any_unreliable {
            lines.push(String::new());
            lines.push("⚠ A fit ran to the BC search limit — the data did not determine a real".to_string());
            lines.push("  value. Add more/longer-range points and check --zero-range / --temperature.".to_string());
        }
        Ok(lines.join("\n"))
    }

    fn format_trajectory_table(
        &self,
        result: &crate::cli_api::TrajectoryResult,
        zero_distance: Option<f64>,
        units: UnitSystem,
        full: bool,
        sight_height_above_bore_m: f64,
    ) -> String {
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
            if max_range_display < 100.0 {
                10.0
            } else {
                25.0
            }
        } else {
            if max_range_display < 500.0 {
                50.0
            } else {
                100.0
            }
        };

        let mut current_range = 0.0;

        // Get initial height for reference (sight height)
        let los_height = if !result.points.is_empty() {
            result.points[0].position.y + sight_height_above_bore_m
        } else {
            0.05 // Default 2 inches
        };

        for (idx, point) in result.points.iter().enumerate() {
            // McCoy coordinate system: X=downrange, Y=vertical, Z=lateral
            let range_display = match units {
                UnitSystem::Imperial => point.position.x * 1.09361, // X is downrange (McCoy), m to yards
                UnitSystem::Metric => point.position.x,             // X is downrange (McCoy)
            };

            let is_last_point = idx == result.points.len() - 1;

            // Show point if it's at the sampling interval OR if it's the last point OR if it's the zero distance
            let should_show = range_display >= current_range
                || is_last_point
                || (zero_distance.is_some()
                    && (range_display - zero_distance.unwrap()).abs() < 1.0);

            if should_show {
                let drop = los_height - point.position.y;
                let drift = point.position.z; // Z is lateral (windage, McCoy)
                let velocity = point.velocity_magnitude;

                // Format values based on unit system
                let (range_str, drop_str, drift_str, velocity_str, energy_str) = match units {
                    UnitSystem::Imperial => {
                        let energy_ftlb = 0.5
                            * (result.points[0].kinetic_energy / 0.5)
                            * (velocity / result.points[0].velocity_magnitude).powi(2)
                            * 0.737562149;
                        (
                            format!("{:.0} yd", range_display),
                            format!("{:.1} in", drop * 39.3701),
                            format!("{:.1} in", drift * 39.3701),
                            format!("{:.0} fps", velocity * 3.28084),
                            format!("{:.0} ft-lb", energy_ftlb),
                        )
                    }
                    UnitSystem::Metric => (
                        format!("{:.0} m", range_display),
                        format!("{:.1} cm", drop * 100.0),
                        format!("{:.1} cm", drift * 100.0),
                        format!("{:.0} m/s", velocity),
                        format!("{:.0} J", point.kinetic_energy),
                    ),
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
        output.push_str(&format!(
            "\nMax Range: {:.0} {}\n",
            if units == UnitSystem::Imperial {
                result.max_range * 1.09361
            } else {
                result.max_range
            },
            if units == UnitSystem::Imperial {
                "yards"
            } else {
                "meters"
            }
        ));

        // Check termination reason
        if let Some(last_point) = result.points.last() {
            // Debug info about last point
            // McCoy coordinate system: X=downrange, Y=vertical, Z=lateral
            let last_height = last_point.position.y;
            let last_range = last_point.position.x; // X is downrange (McCoy)
            let last_time = last_point.time;

            // Ground threshold is typically around 0.01m (1cm), so if y is close to or below that, it hit ground
            if last_height <= 0.01 {
                output.push_str(&format!(
                    "Bullet struck ground at {:.0} {}\n",
                    if units == UnitSystem::Imperial {
                        last_range * 1.09361
                    } else {
                        last_range
                    },
                    if units == UnitSystem::Imperial {
                        "yards"
                    } else {
                        "meters"
                    }
                ));
            } else if last_time >= 99.0 {
                output.push_str(&format!(
                    "Trajectory timeout at {:.0} {} (time limit reached)\n",
                    if units == UnitSystem::Imperial {
                        last_range * 1.09361
                    } else {
                        last_range
                    },
                    if units == UnitSystem::Imperial {
                        "yards"
                    } else {
                        "meters"
                    }
                ));
            }
        }

        output.push_str(&format!(
            "Max Height: {:.1} {}\n",
            if units == UnitSystem::Imperial {
                result.max_height * 39.3701
            } else {
                result.max_height * 100.0
            },
            if units == UnitSystem::Imperial {
                "inches"
            } else {
                "cm"
            }
        ));
        output.push_str(&format!(
            "Time of Flight: {:.2} seconds\n",
            result.time_of_flight
        ));
        output.push_str(&format!(
            "Impact Velocity: {:.0} {}\n",
            if units == UnitSystem::Imperial {
                result.impact_velocity * 3.28084
            } else {
                result.impact_velocity
            },
            if units == UnitSystem::Imperial {
                "fps"
            } else {
                "m/s"
            }
        ));

        output
    }

    fn format_trajectory_json(
        &self,
        result: &crate::cli_api::TrajectoryResult,
        units: UnitSystem,
        sight_height_above_bore_m: f64,
    ) -> String {
        let los_height = result
            .points
            .first()
            .map(|p| p.position.y + sight_height_above_bore_m)
            .unwrap_or(0.05);
        // McCoy coordinate system: X=downrange, Y=vertical, Z=lateral
        let points: Vec<serde_json::Value> = result
            .points
            .iter()
            .map(|p| {
                match units {
                    UnitSystem::Imperial => {
                        serde_json::json!({
                            "range_yards": p.position.x * 1.09361,  // X is downrange (McCoy)
                            "drop_inches": (los_height - p.position.y) * 39.3701,
                            "drift_inches": p.position.z * 39.3701,  // Z is lateral (windage, McCoy)
                            "velocity_fps": p.velocity_magnitude * 3.28084,
                            "energy_ftlb": p.kinetic_energy * 0.737562149,
                            "time_seconds": p.time
                        })
                    }
                    UnitSystem::Metric => {
                        serde_json::json!({
                            "range_meters": p.position.x,  // X is downrange (McCoy)
                            "drop_cm": (los_height - p.position.y) * 100.0,
                            "drift_cm": p.position.z * 100.0,  // Z is lateral (windage, McCoy)
                            "velocity_mps": p.velocity_magnitude,
                            "energy_joules": p.kinetic_energy,
                            "time_seconds": p.time
                        })
                    }
                }
            })
            .collect();

        let summary = match units {
            UnitSystem::Imperial => {
                serde_json::json!({
                    "max_range_yards": result.max_range * 1.09361,
                    "max_height_inches": result.max_height * 39.3701,
                    "time_of_flight_seconds": result.time_of_flight,
                    "impact_velocity_fps": result.impact_velocity * 3.28084
                })
            }
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

        serde_json::to_string_pretty(&output)
            .unwrap_or_else(|_| "Error formatting JSON".to_string())
    }

    fn format_trajectory_csv(
        &self,
        result: &crate::cli_api::TrajectoryResult,
        units: UnitSystem,
        full: bool,
        sight_height_above_bore_m: f64,
    ) -> String {
        let mut output = String::new();

        // Header
        match units {
            UnitSystem::Imperial => {
                output.push_str("Range(yards),Drop(inches),Drift(inches),Velocity(fps),Energy(ft-lb),Time(seconds)\n");
            }
            UnitSystem::Metric => {
                output.push_str(
                    "Range(meters),Drop(cm),Drift(cm),Velocity(m/s),Energy(joules),Time(seconds)\n",
                );
            }
        }

        // Determine sampling interval
        let max_range_display = match units {
            UnitSystem::Imperial => result.max_range * 1.09361,
            UnitSystem::Metric => result.max_range,
        };

        let sample_interval = if full {
            if max_range_display < 100.0 {
                10.0
            } else {
                25.0
            }
        } else {
            if max_range_display < 500.0 {
                50.0
            } else {
                100.0
            }
        };

        let mut current_range = 0.0;
        let los_height = if !result.points.is_empty() {
            result.points[0].position.y + sight_height_above_bore_m
        } else {
            0.05
        };

        for (idx, point) in result.points.iter().enumerate() {
            // McCoy coordinate system: X=downrange, Y=vertical, Z=lateral
            let range_display = match units {
                UnitSystem::Imperial => point.position.x * 1.09361, // X is downrange (McCoy)
                UnitSystem::Metric => point.position.x,             // X is downrange (McCoy)
            };

            let is_last_point = idx == result.points.len() - 1;

            if range_display >= current_range || is_last_point {
                let drop = los_height - point.position.y;

                match units {
                    UnitSystem::Imperial => {
                        let energy_ftlb = point.kinetic_energy * 0.737562149;
                        output.push_str(&format!(
                            "{:.1},{:.2},{:.2},{:.0},{:.0},{:.3}\n",
                            range_display,
                            drop * 39.3701,
                            point.position.z * 39.3701, // Z is lateral (windage, McCoy)
                            point.velocity_magnitude * 3.28084,
                            energy_ftlb,
                            point.time
                        ));
                    }
                    UnitSystem::Metric => {
                        output.push_str(&format!(
                            "{:.1},{:.2},{:.2},{:.0},{:.0},{:.3}\n",
                            range_display,
                            drop * 100.0,
                            point.position.z * 100.0, // Z is lateral (windage, McCoy)
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
  lead           Calculate moving-target lead (hold)
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
    --wind-direction <DIR>       Wind direction (deg; 0=headwind, 90=from right)
    --wind-segment <S:A:D>       Downrange wind seg speed:angle:until-dist (repeatable)
    --temperature <TEMP>         Temperature (F/C)
    --pressure <PRESSURE>        Pressure (inHg/hPa)
    --humidity <HUMIDITY>        Humidity (0-100%)
    --altitude <ALT>             Altitude (feet/meters)
    
  Advanced Physics:
    --enable-magnus              Enable Magnus effect
    --enable-coriolis            Enable Coriolis effect
    --enable-spin-drift          Enable empirical Litz spin drift
    --enable-wind-shear          Enable altitude-dependent wind
    --enable-pitch-damping       Enable transonic stability
    --enable-precession          Enable angular motion physics
    --use-euler                  Use Euler integration (less accurate)
    --use-rk4-fixed              Use fixed-step RK4 (default: adaptive RK45)
    --time-step <SECONDS>        Integration time step (seconds) [default: 0.001]
    --sample-trajectory          Enable trajectory sampling
    --sample-interval <DIST>     Trajectory sampling interval (yards/meters) [default: 10]
    --use-bc-segments            Use velocity-based BC (from a loaded BC5D table)
    --bc-segment <VMIN:VMAX:BC>  Manual velocity-keyed BC segment (repeatable; fps/m/s per --units)
    --use-powder-sensitivity     Enable powder temp sensitivity
    
  Additional Parameters:
    --twist-rate <RATE>          Barrel twist (inches/turn imperial, mm/turn metric)
    --twist-right <BOOL>         Right-hand twist (true/false)
    --latitude <LAT>             Latitude for Coriolis (degrees)
    --shot-direction <DEG>       Compass bearing of the shot for Coriolis (0=N, 90=E)
    --shooting-angle <ANGLE>     Uphill/downhill angle (degrees)
    --cant <DEGREES>             Rifle cant angle (degrees); positive = clockwise from the
                                 shooter, moving point of impact right and low
    --sight-height <HEIGHT>      Sight height above bore (inches/mm)
    --muzzle-height <HEIGHT>     Shooter height above ground (inches/mm)
    --target-height <HEIGHT>     Target height above ground (inches/mm)
    --powder-temp <TEMP>         Powder temperature
    --powder-temp-sensitivity <SENS>  Velocity change per degree

Zero Command:
  ballistics zero [OPTIONS]
  
  Options:
    -v, --velocity <VEL>         Muzzle velocity
    -b, --bc <BC>                Ballistic coefficient
    -m, --mass <MASS>            Mass
    -d, --diameter <DIA>         Diameter
    --drag-model <MODEL>         Drag model (G1/G7)
    --target-distance <DIST>     Target distance for zero
    --sight-height <HEIGHT>      Sight height above bore

Monte Carlo Command:
  ballistics monte-carlo [OPTIONS]
  
  Options:
    -v, --velocity <VEL>         Base velocity
    -b, --bc <BC>                Base BC
    -m, --mass <MASS>            Mass
    -d, --diameter <DIA>         Diameter
    --drag-model <MODEL>         Drag model (G1/G7)
    -n, --num-sims <N>           Number of simulations
    --velocity-std <STD>         Velocity std deviation
    --angle-std <STD>            Angle std deviation
    --bc-std <STD>               BC std deviation
    --wind-speed-std <STD>       Wind speed std deviation
    --wind-direction-std <STD>   Wind direction std deviation (degrees)

Estimate BC Command:
  ballistics estimate-bc [OPTIONS]

  Options:
    -v, --velocity <VEL>         Muzzle velocity (fps/m/s)
    -m, --mass <MASS>            Mass (grains/grams)
    -d, --diameter <DIA>         Diameter (inches/mm)
    --data <PAIRS>               Drop data: "dist,drop;..." (yd,in / m,mm)
    --velocity-data <PAIRS>      Velocity data: "dist,vel;..." (yd,fps / m,m/s)
    --drag-model <MODEL>         g1, g7, or both [default: both]
    --zero-range <DIST>          Zero range of the drop data (yd/m). Dope cards are
                                 zeroed — pass this so drop is fit below line of sight.
    --sight-height <H>           Sight height above bore (in/mm) for the zeroed fit
    --temperature <T>            Air temp the data was measured at (°F/°C) [59/15]
    --pressure <P>               Pressure the data was measured at (inHg/hPa) [29.92/1013.25]
    --humidity <H>               Relative humidity (percent) [50]
    --altitude <A>               Altitude the data was measured at (ft/m) [0]

  Prints a BC for each drag model x data basis you supply. For a DOPE CARD, pass
  --zero-range (drop is below line of sight) and set the atmosphere it was made at.
  Without --zero-range, drop is treated as bore-referenced (flat-fire). --velocity-data
  gives a velocity-retention fit (immune to zero/angle). A fit that can't be pinned
  down is flagged UNRELIABLE.

Lead Command:
  ballistics lead --target-speed <SPEED> [OPTIONS]

  Options:
    -v, --velocity <VEL>          Muzzle velocity (fps/m/s)
    -b, --bc <BC>                 Ballistic coefficient
    -m, --mass <MASS>             Mass (grains/grams)
    -d, --diameter <DIA>          Diameter (inches/mm)
    --drag-model <MODEL>          Drag model (G1/G7)
    --sight-height <HEIGHT>       Sight height above bore (inches/mm)
    --target-speed <SPEED>        Target speed (mph/m/s) [required]
    --target-angle <DEG>          Direction of target travel, degrees [default: 90]
                                  0=away, 90=left-to-right, 180=toward, 270=right-to-left;
                                  positive lead = hold in direction of travel
    --range <DIST>                Range to target (yards/meters) [default: 500]
    --adjustment-unit <UNIT>      mil or moa [default: mil]

  Note: this command assumes calm air and standard atmosphere; for a wind-aware
  lead (time of flight under wind), use the native CLI's lead subcommand.

Examples:
  ballistics trajectory -v 2700 -b 0.475 -m 168 -d 0.308
  ballistics trajectory --auto-zero 200 --enable-spin-drift
  ballistics --units metric trajectory -v 823 -b 0.475 -m 10.9
  ballistics zero --target-distance 300
  ballistics estimate-bc -v 2700 -m 168 -d 0.308 --data "100,2.1;200,9.4;300,22.8"
  ballistics estimate-bc -v 2650 -m 77 -d 0.224 --data "300,29;500,89.9" \
    --velocity-data "300,1980;500,1560" --drag-model both
  ballistics monte-carlo -n 1000 --velocity-std 10"#
            .to_string()
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
    // Downrange wind segments as (speed_mph, angle_deg, until_yards); when non-empty
    // these are emitted as --wind-segment flags and override the scalar wind above.
    wind_segments: Vec<(f64, f64, f64)>,
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
            drag_model: "G1".to_string(), // G1 matches the G1-scale default BC (0.475) and the CLI default

            wind_speed_mph: 0.0,
            wind_direction_deg: 0.0, // wind-FROM convention: 0 = headwind, matching the CLI/trajectory defaults
            wind_segments: Vec::new(),
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

    /// Add a downrange wind segment: `speed_mph` from `direction_deg`
    /// (wind-FROM: 0 = headwind, 90 = from the right) out to `until_yards`.
    /// Each segment applies from the previous boundary to `until_yards`; wind is
    /// zero beyond the last segment. When any segment is added it overrides the
    /// scalar `setWind` value. Repeatable.
    #[wasm_bindgen(js_name = addWindSegment)]
    pub fn add_wind_segment(
        mut self,
        speed_mph: f64,
        direction_deg: f64,
        until_yards: f64,
    ) -> Self {
        self.wind_segments
            .push((speed_mph, direction_deg, until_yards));
        self
    }

    /// Remove all downrange wind segments (reverts to the scalar `setWind`).
    #[wasm_bindgen(js_name = clearWindSegments)]
    pub fn clear_wind_segments(mut self) -> Self {
        self.wind_segments.clear();
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
            cmd.push_str(&format!(
                " --wind-speed {} --wind-direction {}",
                self.wind_speed_mph, self.wind_direction_deg
            ));
        }
        for (speed_mph, direction_deg, until_yards) in &self.wind_segments {
            cmd.push_str(&format!(
                " --wind-segment {}:{}:{}",
                speed_mph, direction_deg, until_yards
            ));
        }

        cmd.push_str(&format!(
            " --temperature {} --pressure {} --humidity {} --altitude {}",
            self.temperature_f, self.pressure_inhg, self.humidity_percent, self.altitude_ft
        ));

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
        let json_result: serde_json::Value = serde_json::from_str(json_str).map_err(|e| {
            JsValue::from_str(&format!(
                "JSON parse error: {}. Result was: {}",
                e,
                &json_str[..json_str.len().min(500)]
            ))
        })?;

        // Find the point closest to the requested range
        if let Some(trajectory) = json_result.get("trajectory").and_then(|t| t.as_array()) {
            let target_point = trajectory
                .iter()
                .min_by_key(|p| {
                    let range = p.get("range_yards").and_then(|r| r.as_f64()).unwrap_or(0.0);
                    ((range - range_yards).abs() * 100.0) as i64
                })
                .ok_or_else(|| JsValue::from_str("No trajectory points found"))?;

            // Manually construct JavaScript object to avoid Map conversion
            let result = js_sys::Object::new();
            js_sys::Reflect::set(
                &result,
                &"range_yards".into(),
                &target_point
                    .get("range_yards")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0)
                    .into(),
            )?;
            js_sys::Reflect::set(
                &result,
                &"drop_inches".into(),
                &target_point
                    .get("drop_inches")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0)
                    .into(),
            )?;
            js_sys::Reflect::set(
                &result,
                &"drift_inches".into(),
                &target_point
                    .get("drift_inches")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0)
                    .into(),
            )?;
            js_sys::Reflect::set(
                &result,
                &"velocity_fps".into(),
                &target_point
                    .get("velocity_fps")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0)
                    .into(),
            )?;
            js_sys::Reflect::set(
                &result,
                &"energy_ftlb".into(),
                &target_point
                    .get("energy_ftlb")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0)
                    .into(),
            )?;
            js_sys::Reflect::set(
                &result,
                &"time_seconds".into(),
                &target_point
                    .get("time_seconds")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0)
                    .into(),
            )?;

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
            cmd.push_str(&format!(
                " --wind-speed {} --wind-direction {}",
                self.wind_speed_mph, self.wind_direction_deg
            ));
        }
        for (speed_mph, direction_deg, until_yards) in &self.wind_segments {
            cmd.push_str(&format!(
                " --wind-segment {}:{}:{}",
                speed_mph, direction_deg, until_yards
            ));
        }

        cmd.push_str(&format!(
            " --temperature {} --pressure {} --humidity {} --altitude {}",
            self.temperature_f, self.pressure_inhg, self.humidity_percent, self.altitude_ft
        ));

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
        let json_result: serde_json::Value = serde_json::from_str(json_str).map_err(|e| {
            JsValue::from_str(&format!(
                "JSON parse error: {}. Result: {}",
                e,
                &json_str[..json_str.len().min(500)]
            ))
        })?;

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
                    // The JSON producer emits "drift_inches"/"time_seconds"; read those (the old
                    // "windage_inches"/"time_sec" lookups always missed, dropping both fields).
                    // Keep the public output keys (windage_inches/time_sec) unchanged.
                    if let Some(windage) = point.get("drift_inches").and_then(|v| v.as_f64()) {
                        js_sys::Reflect::set(&js_point, &"windage_inches".into(), &windage.into())?;
                    }
                    if let Some(velocity) = point.get("velocity_fps").and_then(|v| v.as_f64()) {
                        js_sys::Reflect::set(&js_point, &"velocity_fps".into(), &velocity.into())?;
                    }
                    if let Some(energy) = point.get("energy_ftlb").and_then(|v| v.as_f64()) {
                        js_sys::Reflect::set(&js_point, &"energy_ftlb".into(), &energy.into())?;
                    }
                    if let Some(time) = point.get("time_seconds").and_then(|v| v.as_f64()) {
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
