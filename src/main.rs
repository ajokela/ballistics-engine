use ballistics_engine::{
    BallisticInputs, TrajectorySolver, MonteCarloParams, 
    WindConditions, AtmosphericConditions, DragModel
};
use clap::{Parser, Subcommand, ValueEnum};
use serde::{Serialize, Deserialize};
use std::error::Error;

#[derive(Parser)]
#[command(name = "ballistics")]
#[command(author = "Ballistics Engine Team")]
#[command(version = "0.1.0")]
#[command(about = "High-performance ballistics trajectory calculator", long_about = None)]
struct Cli {
    /// Unit system for input/output (metric or imperial)
    #[arg(short = 'u', long, default_value = "imperial", global = true)]
    units: UnitSystem,
    
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Calculate a single trajectory
    Trajectory {
        /// Initial velocity (fps or m/s based on --units)
        #[arg(short = 'v', long)]
        velocity: f64,
        
        /// Launch angle (degrees)
        #[arg(short = 'a', long, default_value = "0.0")]
        angle: f64,
        
        /// Ballistic coefficient
        #[arg(short = 'b', long)]
        bc: f64,
        
        /// Mass (grains or kg based on --units)
        #[arg(short = 'm', long)]
        mass: f64,
        
        /// Diameter (inches or meters based on --units)
        #[arg(short = 'd', long)]
        diameter: f64,
        
        /// Drag model (G1, G7, Custom)
        #[arg(long, default_value = "g1")]
        drag_model: DragModelArg,
        
        /// Maximum range (yards or meters based on --units)
        #[arg(long, default_value = "1000.0")]
        max_range: f64,
        
        /// Time step (seconds)
        #[arg(long, default_value = "0.001")]
        time_step: f64,
        
        /// Wind speed (mph or m/s based on --units)
        #[arg(long, default_value = "0.0")]
        wind_speed: f64,
        
        /// Wind direction (degrees, 0=North, 90=East)
        #[arg(long, default_value = "0.0")]
        wind_direction: f64,
        
        /// Temperature (Fahrenheit or Celsius based on --units)
        #[arg(long, default_value = "59.0")]
        temperature: f64,
        
        /// Pressure (inHg or hPa based on --units)
        #[arg(long, default_value = "29.92")]
        pressure: f64,
        
        /// Humidity (0-100%)
        #[arg(long, default_value = "50.0")]
        humidity: f64,
        
        /// Altitude (feet or meters based on --units)
        #[arg(long, default_value = "0.0")]
        altitude: f64,
        
        /// Output format
        #[arg(short = 'o', long, default_value = "table")]
        output: OutputFormat,
        
        /// Show all trajectory points
        #[arg(long)]
        full: bool,
        
        /// Automatically zero to target distance (overrides --angle)
        #[arg(long)]
        auto_zero: Option<f64>,
        
        /// Sight height above bore for auto-zero
        #[arg(long, default_value = "0.05")]
        sight_height: f64,
        
        /// Enable velocity-based BC segmentation
        #[arg(long)]
        use_bc_segments: bool,
        
        /// Enable cluster-based BC degradation (overrides BC segments)
        #[arg(long)]
        use_cluster_bc: bool,
        
        /// Bullet cluster ID (0-3) for cluster BC
        #[arg(long)]
        bullet_cluster: Option<usize>,
        
        // Advanced Physics Parameters
        
        /// Enable Magnus effect (requires twist-rate)
        #[arg(long)]
        enable_magnus: bool,
        
        /// Enable Coriolis effect (requires latitude)
        #[arg(long)]
        enable_coriolis: bool,
        
        /// Enable enhanced spin drift calculations
        #[arg(long)]
        enable_spin_drift: bool,
        
        /// Enable wind shear (altitude-dependent wind)
        #[arg(long)]
        enable_wind_shear: bool,
        
        /// Barrel twist rate (inches per turn, e.g., 10 for 1:10)
        #[arg(long)]
        twist_rate: Option<f64>,
        
        /// Right-hand twist (true) or left-hand (false)
        #[arg(long, default_value = "true")]
        twist_right: bool,
        
        /// Latitude for Coriolis effect (degrees, -90 to 90)
        #[arg(long)]
        latitude: Option<f64>,
        
        /// Shooting angle (degrees, positive = uphill, negative = downhill)
        #[arg(long, default_value = "0.0")]
        shooting_angle: f64,
        
        /// Enable powder temperature sensitivity
        #[arg(long)]
        use_powder_sensitivity: bool,
        
        /// Powder temperature sensitivity (fps per degree)
        #[arg(long, default_value = "1.0")]
        powder_temp_sensitivity: f64,
        
        /// Powder temperature
        #[arg(long, default_value = "70.0")]
        powder_temp: f64,
    },
    
    /// Run Monte Carlo simulation
    MonteCarlo {
        /// Base velocity (m/s)
        #[arg(short = 'v', long)]
        velocity: f64,
        
        /// Launch angle (degrees)
        #[arg(short = 'a', long, default_value = "0.0")]
        angle: f64,
        
        /// Ballistic coefficient
        #[arg(short = 'b', long)]
        bc: f64,
        
        /// Mass (kg)
        #[arg(short = 'm', long)]
        mass: f64,
        
        /// Diameter (meters)
        #[arg(short = 'd', long)]
        diameter: f64,
        
        /// Number of simulations
        #[arg(short = 'n', long, default_value = "1000")]
        num_sims: usize,
        
        /// Velocity standard deviation (m/s)
        #[arg(long, default_value = "1.0")]
        velocity_std: f64,
        
        /// Angle standard deviation (degrees)
        #[arg(long, default_value = "0.1")]
        angle_std: f64,
        
        /// BC standard deviation
        #[arg(long, default_value = "0.01")]
        bc_std: f64,
        
        /// Wind speed standard deviation (m/s)
        #[arg(long, default_value = "1.0")]
        wind_std: f64,
        
        /// Target distance (meters)
        #[arg(long)]
        target_distance: Option<f64>,
        
        /// Output format
        #[arg(short = 'o', long, default_value = "summary")]
        output: MonteCarloOutput,
    },
    
    /// Calculate zero angle for a target
    Zero {
        /// Initial velocity (m/s)
        #[arg(short = 'v', long)]
        velocity: f64,
        
        /// Ballistic coefficient
        #[arg(short = 'b', long)]
        bc: f64,
        
        /// Mass (kg)
        #[arg(short = 'm', long)]
        mass: f64,
        
        /// Diameter (meters)
        #[arg(short = 'd', long)]
        diameter: f64,
        
        /// Target distance (meters)
        #[arg(long)]
        target_distance: f64,
        
        /// Target height (meters, negative for below)
        #[arg(long, default_value = "0.0")]
        target_height: f64,
        
        /// Sight height above bore (meters)
        #[arg(long, default_value = "0.05")]
        sight_height: f64,
        
        /// Output format
        #[arg(short = 'o', long, default_value = "table")]
        output: OutputFormat,
    },
    
    /// Estimate BC from trajectory data
    EstimateBC {
        /// Initial velocity (m/s)
        #[arg(short = 'v', long)]
        velocity: f64,
        
        /// Mass (kg)
        #[arg(short = 'm', long)]
        mass: f64,
        
        /// Diameter (meters)
        #[arg(short = 'd', long)]
        diameter: f64,
        
        /// Distance 1 (meters)
        #[arg(long)]
        distance1: f64,
        
        /// Drop at distance 1 (meters)
        #[arg(long)]
        drop1: f64,
        
        /// Distance 2 (meters)
        #[arg(long)]
        distance2: f64,
        
        /// Drop at distance 2 (meters)
        #[arg(long)]
        drop2: f64,
        
        /// Output format
        #[arg(short = 'o', long, default_value = "table")]
        output: OutputFormat,
    },
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum UnitSystem {
    /// Metric units (m/s, kg, meters, Celsius)
    Metric,
    /// Imperial units (fps, grains, yards, Fahrenheit)
    Imperial,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum DragModelArg {
    G1,
    G7,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum OutputFormat {
    Json,
    Csv,
    Table,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum MonteCarloOutput {
    Summary,
    Full,
    Statistics,
}

#[derive(Debug, Serialize, Deserialize)]
struct TrajectoryPoint {
    time: f64,
    x: f64,
    y: f64,
    z: f64,
    velocity: f64,
    energy: f64,
}

#[derive(Debug, Serialize, Deserialize)]
struct TrajectoryResult {
    max_range: f64,
    max_height: f64,
    time_of_flight: f64,
    impact_velocity: f64,
    impact_energy: f64,
    trajectory: Vec<TrajectoryPoint>,
}

#[derive(Debug, Serialize, Deserialize)]
struct MonteCarloResult {
    mean_range: f64,
    std_range: f64,
    mean_impact_velocity: f64,
    std_impact_velocity: f64,
    cep: f64,  // Circular Error Probable
    hit_probability: Option<f64>,
}

// Unit conversion functions
struct UnitConverter;

impl UnitConverter {
    // Input conversions (to metric)
    fn velocity_to_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val * 0.3048, // fps to m/s
        }
    }
    
    fn mass_to_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val * 0.00006479891, // grains to kg
        }
    }
    
    fn distance_to_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val * 0.9144, // yards to meters
        }
    }
    
    fn diameter_to_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val * 0.0254, // inches to meters
        }
    }
    
    fn wind_to_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val * 0.44704, // mph to m/s
        }
    }
    
    fn temperature_to_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => (val - 32.0) * 5.0 / 9.0, // Fahrenheit to Celsius
        }
    }
    
    fn pressure_to_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val * 33.8639, // inHg to hPa
        }
    }
    
    fn altitude_to_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val * 0.3048, // feet to meters
        }
    }
    
    // Output conversions (from metric)
    fn velocity_from_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val / 0.3048, // m/s to fps
        }
    }
    
    fn distance_from_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val / 0.9144, // meters to yards
        }
    }
    
    fn energy_from_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val * 0.737562, // Joules to ft-lbs
        }
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();
    
    match cli.command {
        Commands::Trajectory { 
            velocity, angle, bc, mass, diameter, drag_model, 
            max_range, time_step, wind_speed, wind_direction,
            temperature, pressure, humidity, altitude,
            output, full, auto_zero, sight_height,
            use_bc_segments, use_cluster_bc, bullet_cluster,
            enable_magnus, enable_coriolis, enable_spin_drift,
            enable_wind_shear, twist_rate, twist_right, latitude,
            shooting_angle, use_powder_sensitivity, 
            powder_temp_sensitivity, powder_temp
        } => {
            // Convert inputs to metric
            let velocity_metric = UnitConverter::velocity_to_metric(velocity, cli.units);
            let mass_metric = UnitConverter::mass_to_metric(mass, cli.units);
            let diameter_metric = UnitConverter::diameter_to_metric(diameter, cli.units);
            let max_range_metric = UnitConverter::distance_to_metric(max_range, cli.units);
            let wind_speed_metric = UnitConverter::wind_to_metric(wind_speed, cli.units);
            let temperature_metric = UnitConverter::temperature_to_metric(temperature, cli.units);
            let pressure_metric = UnitConverter::pressure_to_metric(pressure, cli.units);
            let altitude_metric = UnitConverter::altitude_to_metric(altitude, cli.units);
            let sight_height_metric = UnitConverter::distance_to_metric(sight_height, cli.units);
            
            // Calculate zero angle if auto-zero is specified
            let launch_angle = if let Some(zero_distance) = auto_zero {
                let zero_distance_metric = UnitConverter::distance_to_metric(zero_distance, cli.units);
                
                // Create inputs for zero calculation
                let zero_inputs = BallisticInputs {
                    muzzle_velocity: velocity_metric,
                    ballistic_coefficient: bc,
                    mass: mass_metric,
                    diameter: diameter_metric,
                    sight_height: sight_height_metric,
                    ..Default::default()
                };
                
                // Calculate zero angle
                let zero_angle = ballistics_engine::calculate_zero_angle(
                    zero_inputs,
                    zero_distance_metric,
                    0.0  // target height at zero distance
                )?;
                
                // Convert to degrees for the trajectory function
                zero_angle.to_degrees()
            } else {
                angle
            };
            
            run_trajectory(
                velocity_metric, launch_angle, bc, mass_metric, diameter_metric, drag_model,
                max_range_metric, time_step, wind_speed_metric, wind_direction,
                temperature_metric, pressure_metric, humidity, altitude_metric,
                output, full, cli.units, sight_height_metric,
                use_bc_segments, use_cluster_bc, bullet_cluster,
                enable_magnus, enable_coriolis, enable_spin_drift,
                enable_wind_shear, twist_rate, twist_right, latitude,
                shooting_angle, use_powder_sensitivity,
                powder_temp_sensitivity, powder_temp
            )?;
        },
        
        Commands::MonteCarlo {
            velocity, angle, bc, mass, diameter,
            num_sims, velocity_std, angle_std, bc_std, wind_std,
            target_distance, output
        } => {
            run_monte_carlo(
                velocity, angle, bc, mass, diameter,
                num_sims, velocity_std, angle_std, bc_std, wind_std,
                target_distance, output
            )?;
        },
        
        Commands::Zero {
            velocity, bc, mass, diameter,
            target_distance, target_height, sight_height,
            output
        } => {
            // Convert inputs to metric
            let velocity_metric = UnitConverter::velocity_to_metric(velocity, cli.units);
            let mass_metric = UnitConverter::mass_to_metric(mass, cli.units);
            let diameter_metric = UnitConverter::diameter_to_metric(diameter, cli.units);
            let target_distance_metric = UnitConverter::distance_to_metric(target_distance, cli.units);
            let target_height_metric = UnitConverter::distance_to_metric(target_height, cli.units);
            let sight_height_metric = UnitConverter::distance_to_metric(sight_height, cli.units);
            
            run_zero_calculation(
                velocity_metric, bc, mass_metric, diameter_metric,
                target_distance_metric, target_height_metric, sight_height_metric,
                output, cli.units
            )?;
        },
        
        Commands::EstimateBC {
            velocity, mass, diameter,
            distance1, drop1, distance2, drop2,
            output
        } => {
            run_bc_estimation(
                velocity, mass, diameter,
                distance1, drop1, distance2, drop2,
                output
            )?;
        },
    }
    
    Ok(())
}

fn run_trajectory(
    velocity: f64,
    angle: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    drag_model: DragModelArg,
    max_range: f64,
    time_step: f64,
    wind_speed: f64,
    wind_direction: f64,
    temperature: f64,
    pressure: f64,
    humidity: f64,
    altitude: f64,
    output: OutputFormat,
    full: bool,
    units: UnitSystem,
    sight_height: f64,
    use_bc_segments: bool,
    use_cluster_bc: bool,
    bullet_cluster: Option<usize>,
    enable_magnus: bool,
    enable_coriolis: bool,
    enable_spin_drift: bool,
    enable_wind_shear: bool,
    twist_rate: Option<f64>,
    twist_right: bool,
    latitude: Option<f64>,
    shooting_angle: f64,
    use_powder_sensitivity: bool,
    powder_temp_sensitivity: f64,
    powder_temp: f64,
) -> Result<(), Box<dyn Error>> {
    // Create ballistic inputs with all required fields
    let drag_model_enum = match drag_model {
        DragModelArg::G1 => DragModel::G1,
        DragModelArg::G7 => DragModel::G7,
    };
    
    let inputs = BallisticInputs {
        // Core fields
        muzzle_velocity: velocity,
        launch_angle: angle.to_radians(),
        ballistic_coefficient: bc,
        mass,
        diameter,
        drag_model: drag_model_enum,
        sight_height,  // Use provided sight height
        
        // Duplicate fields for internal compatibility
        altitude,
        bc_type: drag_model_enum,
        bc_value: bc,
        caliber_inches: diameter / 0.0254,  // Convert meters to inches
        weight_grains: mass / 0.00006479891,  // Convert kg to grains
        bullet_diameter: diameter,  // Keep in meters
        bullet_mass: mass,  // Keep in kg
        bullet_length: diameter * 4.0,  // Approximate
        muzzle_angle: angle.to_radians(),
        target_distance: max_range,
        temperature,
        twist_rate: twist_rate.unwrap_or(12.0),  // Default 1:12" twist if not specified
        is_twist_right: twist_right,
        shooting_angle: shooting_angle.to_radians(),
        latitude,
        ground_threshold: -10.0,
        
        // Advanced effects - now separately controlled
        enable_advanced_effects: enable_magnus || enable_coriolis,  // Either one enables the system
        use_powder_sensitivity,
        powder_temp_sensitivity: if use_powder_sensitivity { 
            UnitConverter::velocity_to_metric(powder_temp_sensitivity, units) / 
            UnitConverter::temperature_to_metric(1.0, units) 
        } else { 0.0 },
        powder_temp: UnitConverter::temperature_to_metric(powder_temp, units),
        tipoff_yaw: 0.0,
        tipoff_decay_distance: 50.0,
        use_bc_segments,
        bc_segments: None,
        bc_segments_data: None,
        use_enhanced_spin_drift: enable_spin_drift,
        use_form_factor: false,
        use_cluster_bc,
        enable_wind_shear,
        wind_shear_model: if enable_wind_shear { "exponential".to_string() } else { "none".to_string() },
        
        // Optional data
        bc_type_str: None,
        bullet_model: None,
        bullet_id: None,
        bullet_cluster: bullet_cluster.map(|id| id.to_string()),
    };
    
    // Set up wind conditions
    let wind = WindConditions {
        speed: wind_speed,
        direction: wind_direction.to_radians(),
        ..Default::default()
    };
    
    // Set up atmospheric conditions
    let atmosphere = AtmosphericConditions {
        temperature,
        pressure,
        humidity,
        altitude,
        ..Default::default()
    };
    
    // Create solver
    let mut solver = TrajectorySolver::new(inputs, wind, atmosphere);
    solver.set_max_range(max_range);
    solver.set_time_step(time_step);
    
    // Solve trajectory
    let result = solver.solve()?;
    
    // Format output
    match output {
        OutputFormat::Json => {
            let trajectory_result = TrajectoryResult {
                max_range: result.max_range,
                max_height: result.max_height,
                time_of_flight: result.time_of_flight,
                impact_velocity: result.impact_velocity,
                impact_energy: result.impact_energy,
                trajectory: if full {
                    result.points.into_iter().map(|p| TrajectoryPoint {
                        time: p.time,
                        x: p.position.x,
                        y: p.position.y,
                        z: p.position.z,
                        velocity: p.velocity_magnitude,
                        energy: p.kinetic_energy,
                    }).collect()
                } else {
                    vec![]
                },
            };
            println!("{}", serde_json::to_string_pretty(&trajectory_result)?);
        },
        
        OutputFormat::Csv => {
            if full {
                println!("time,x,y,z,velocity,energy");
                for p in result.points {
                    println!("{:.4},{:.2},{:.2},{:.2},{:.2},{:.2}",
                        p.time, p.position.x, p.position.y, p.position.z,
                        p.velocity_magnitude, p.kinetic_energy);
                }
            } else {
                println!("metric,value");
                println!("max_range,{:.2}", result.max_range);
                println!("max_height,{:.2}", result.max_height);
                println!("time_of_flight,{:.4}", result.time_of_flight);
                println!("impact_velocity,{:.2}", result.impact_velocity);
                println!("impact_energy,{:.2}", result.impact_energy);
            }
        },
        
        OutputFormat::Table => {
            // Convert outputs to user's units
            let range_display = UnitConverter::distance_from_metric(result.max_range, units);
            let height_display = UnitConverter::distance_from_metric(result.max_height, units);
            let velocity_display = UnitConverter::velocity_from_metric(result.impact_velocity, units);
            let energy_display = UnitConverter::energy_from_metric(result.impact_energy, units);
            
            let (range_unit, velocity_unit, energy_unit) = match units {
                UnitSystem::Metric => ("m", "m/s", "J"),
                UnitSystem::Imperial => ("yd", "fps", "ft-lb"),
            };
            
            println!("╔════════════════════════════════════════╗");
            println!("║         TRAJECTORY RESULTS             ║");
            println!("╠════════════════════════════════════════╣");
            println!("║ Max Range:         {:>8.2} {:3}       ║", range_display, range_unit);
            println!("║ Max Height:        {:>8.2} {:3}       ║", height_display, range_unit);
            println!("║ Time of Flight:    {:>8.3} s          ║", result.time_of_flight);
            println!("║ Impact Velocity:   {:>8.2} {:3}       ║", velocity_display, velocity_unit);
            println!("║ Impact Energy:     {:>8.2} {:5}     ║", energy_display, energy_unit);
            println!("╚════════════════════════════════════════╝");
            
            if full && !result.points.is_empty() {
                println!();
                println!("Trajectory Points:");
                let (dist_hdr, vel_hdr, energy_hdr) = match units {
                    UnitSystem::Metric => ("(m)", "(m/s)", "(J)"),
                    UnitSystem::Imperial => ("(yd)", "(fps)", "(ft-lb)"),
                };
                println!("┌──────────┬──────────┬──────────┬──────────┬──────────┐");
                println!("│ Time (s) │  X {:5} │  Y {:5} │ Vel{:5} │Energy{:5}│", dist_hdr, dist_hdr, vel_hdr, energy_hdr);
                println!("├──────────┼──────────┼──────────┼──────────┼──────────┤");
                
                let step = if result.points.len() > 20 {
                    result.points.len() / 20
                } else {
                    1
                };
                
                for (i, p) in result.points.iter().enumerate() {
                    if i % step == 0 || i == result.points.len() - 1 {
                        let x_display = UnitConverter::distance_from_metric(p.position.x, units);
                        let y_display = UnitConverter::distance_from_metric(p.position.y, units);
                        let vel_display = UnitConverter::velocity_from_metric(p.velocity_magnitude, units);
                        let energy_display = UnitConverter::energy_from_metric(p.kinetic_energy, units);
                        
                        println!("│ {:>8.3} │ {:>8.2} │ {:>8.2} │ {:>8.2} │ {:>8.2} │",
                            p.time, x_display, y_display, vel_display, energy_display);
                    }
                }
                println!("└──────────┴──────────┴──────────┴──────────┴──────────┘");
            }
        },
    }
    
    Ok(())
}

fn run_monte_carlo(
    velocity: f64,
    angle: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    num_sims: usize,
    velocity_std: f64,
    angle_std: f64,
    bc_std: f64,
    wind_std: f64,
    target_distance: Option<f64>,
    output: MonteCarloOutput,
) -> Result<(), Box<dyn Error>> {
    // Create base inputs
    let base_inputs = BallisticInputs {
        muzzle_velocity: velocity,
        launch_angle: angle.to_radians(),
        ballistic_coefficient: bc,
        mass,
        diameter,
        ..Default::default()
    };
    
    // Set up Monte Carlo parameters
    let mc_params = MonteCarloParams {
        num_simulations: num_sims,
        velocity_std_dev: velocity_std,
        angle_std_dev: angle_std.to_radians(),
        bc_std_dev: bc_std,
        wind_speed_std_dev: wind_std,
        target_distance,
        ..Default::default()
    };
    
    // Run Monte Carlo simulation
    let results = ballistics_engine::run_monte_carlo(base_inputs, mc_params)?;
    
    // Calculate statistics
    let mean_range = results.ranges.iter().sum::<f64>() / results.ranges.len() as f64;
    let variance_range: f64 = results.ranges.iter()
        .map(|r| (r - mean_range).powi(2))
        .sum::<f64>() / results.ranges.len() as f64;
    let std_range = variance_range.sqrt();
    
    let mean_velocity = results.impact_velocities.iter().sum::<f64>() / results.impact_velocities.len() as f64;
    let variance_velocity: f64 = results.impact_velocities.iter()
        .map(|v| (v - mean_velocity).powi(2))
        .sum::<f64>() / results.impact_velocities.len() as f64;
    let std_velocity = variance_velocity.sqrt();
    
    // Calculate CEP (simplified - actual CEP calculation would need lateral dispersion)
    let cep = std_range * 1.1774; // Approximation for circular normal distribution
    
    // Calculate hit probability if target distance specified
    let hit_probability = target_distance.map(|target| {
        let hits = results.ranges.iter()
            .filter(|r| (*r - target).abs() < 1.0) // Within 1m of target
            .count();
        hits as f64 / results.ranges.len() as f64
    });
    
    match output {
        MonteCarloOutput::Summary => {
            println!("╔════════════════════════════════════════╗");
            println!("║      MONTE CARLO RESULTS               ║");
            println!("║      {} simulations                   ║", num_sims);
            println!("╠════════════════════════════════════════╣");
            println!("║ Mean Range:        {:>8.2} m          ║", mean_range);
            println!("║ Std Dev Range:     {:>8.2} m          ║", std_range);
            println!("║ Mean Impact Vel:   {:>8.2} m/s        ║", mean_velocity);
            println!("║ Std Dev Velocity:  {:>8.2} m/s        ║", std_velocity);
            println!("║ CEP (approx):      {:>8.2} m          ║", cep);
            if let Some(prob) = hit_probability {
                println!("║ Hit Probability:   {:>8.1} %          ║", prob * 100.0);
            }
            println!("╚════════════════════════════════════════╝");
        },
        
        MonteCarloOutput::Full => {
            let mc_result = MonteCarloResult {
                mean_range,
                std_range,
                mean_impact_velocity: mean_velocity,
                std_impact_velocity: std_velocity,
                cep,
                hit_probability,
            };
            println!("{}", serde_json::to_string_pretty(&mc_result)?);
        },
        
        MonteCarloOutput::Statistics => {
            println!("range_min,range_max,range_mean,range_std,vel_min,vel_max,vel_mean,vel_std");
            let range_min = results.ranges.iter().fold(f64::INFINITY, |a, &b| a.min(b));
            let range_max = results.ranges.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
            let vel_min = results.impact_velocities.iter().fold(f64::INFINITY, |a, &b| a.min(b));
            let vel_max = results.impact_velocities.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
            
            println!("{:.2},{:.2},{:.2},{:.2},{:.2},{:.2},{:.2},{:.2}",
                range_min, range_max, mean_range, std_range,
                vel_min, vel_max, mean_velocity, std_velocity);
        },
    }
    
    Ok(())
}

fn run_zero_calculation(
    velocity: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    target_distance: f64,
    target_height: f64,
    sight_height: f64,
    output: OutputFormat,
    units: UnitSystem,
) -> Result<(), Box<dyn Error>> {
    // Create ballistic inputs
    let inputs = BallisticInputs {
        muzzle_velocity: velocity,
        ballistic_coefficient: bc,
        mass,
        diameter,
        sight_height,
        ..Default::default()
    };
    
    // Calculate zero angle
    let zero_angle = ballistics_engine::calculate_zero_angle(
        inputs.clone(),
        target_distance,
        target_height
    )?;
    
    // Calculate trajectory at zero angle to get additional info
    let mut zeroed_inputs = inputs;
    zeroed_inputs.launch_angle = zero_angle;
    
    let solver = TrajectorySolver::new(zeroed_inputs, Default::default(), Default::default());
    let trajectory = solver.solve()?;
    
    match output {
        OutputFormat::Json => {
            let result = serde_json::json!({
                "zero_angle_degrees": zero_angle.to_degrees(),
                "zero_angle_moa": zero_angle.to_degrees() * 60.0,
                "zero_angle_mrad": zero_angle * 1000.0,
                "sight_adjustment_moa": (zero_angle.to_degrees() * 60.0) - (sight_height / target_distance * 3437.75),
                "max_ordinate": trajectory.max_height,
                "point_blank_range": trajectory.points.iter()
                    .find(|p| p.position.y < -0.05)
                    .map(|p| p.position.x)
                    .unwrap_or(trajectory.max_range),
            });
            println!("{}", serde_json::to_string_pretty(&result)?);
        },
        
        OutputFormat::Csv => {
            println!("metric,value,unit");
            println!("zero_angle,{:.4},degrees", zero_angle.to_degrees());
            println!("zero_angle_moa,{:.2},MOA", zero_angle.to_degrees() * 60.0);
            println!("zero_angle_mrad,{:.2},mrad", zero_angle * 1000.0);
            println!("max_ordinate,{:.3},meters", trajectory.max_height);
        },
        
        OutputFormat::Table => {
            // Convert distances back to display units
            let target_dist_display = UnitConverter::distance_from_metric(target_distance, units);
            let target_height_display = UnitConverter::distance_from_metric(target_height, units);
            let sight_height_display = UnitConverter::distance_from_metric(sight_height, units);
            let max_ordinate_display = UnitConverter::distance_from_metric(trajectory.max_height, units);
            
            let dist_unit = match units {
                UnitSystem::Metric => "m",
                UnitSystem::Imperial => "yd",
            };
            
            println!("╔════════════════════════════════════════╗");
            println!("║          ZERO CALCULATION              ║");
            println!("╠════════════════════════════════════════╣");
            println!("║ Target Distance:   {:>8.1} {:3}       ║", target_dist_display, dist_unit);
            println!("║ Target Height:     {:>8.2} {:3}       ║", target_height_display, dist_unit);
            println!("║ Sight Height:      {:>8.3} {:3}       ║", sight_height_display, dist_unit);
            println!("╠════════════════════════════════════════╣");
            println!("║ Zero Angle:        {:>8.4}°          ║", zero_angle.to_degrees());
            println!("║ Zero Angle (MOA):  {:>8.2} MOA        ║", zero_angle.to_degrees() * 60.0);
            println!("║ Zero Angle (mrad): {:>8.2} mrad       ║", zero_angle * 1000.0);
            println!("║ Max Ordinate:      {:>8.3} {:3}       ║", max_ordinate_display, dist_unit);
            println!("╚════════════════════════════════════════╝");
        },
    }
    
    Ok(())
}

fn run_bc_estimation(
    velocity: f64,
    mass: f64,
    diameter: f64,
    distance1: f64,
    drop1: f64,
    distance2: f64,
    drop2: f64,
    output: OutputFormat,
) -> Result<(), Box<dyn Error>> {
    // Create trajectory points for BC estimation
    let points = vec![
        (distance1, drop1),
        (distance2, drop2),
    ];
    
    // Estimate BC
    let estimated_bc = ballistics_engine::estimate_bc_from_trajectory(
        velocity,
        mass,
        diameter,
        &points,
    )?;
    
    // Verify the estimation by running a trajectory
    let inputs = BallisticInputs {
        muzzle_velocity: velocity,
        ballistic_coefficient: estimated_bc,
        mass,
        diameter,
        ..Default::default()
    };
    
    let solver = TrajectorySolver::new(inputs, Default::default(), Default::default());
    let trajectory = solver.solve()?;
    
    // Find drops at the specified distances
    let calc_drop1 = trajectory.points.iter()
        .find(|p| p.position.x >= distance1)
        .map(|p| -p.position.y)
        .unwrap_or(0.0);
    
    let calc_drop2 = trajectory.points.iter()
        .find(|p| p.position.x >= distance2)
        .map(|p| -p.position.y)
        .unwrap_or(0.0);
    
    let error1 = ((calc_drop1 - drop1) / drop1 * 100.0).abs();
    let error2 = ((calc_drop2 - drop2) / drop2 * 100.0).abs();
    
    match output {
        OutputFormat::Json => {
            let result = serde_json::json!({
                "estimated_bc": estimated_bc,
                "verification": {
                    "distance1_m": distance1,
                    "actual_drop1_m": drop1,
                    "calculated_drop1_m": calc_drop1,
                    "error1_percent": error1,
                    "distance2_m": distance2,
                    "actual_drop2_m": drop2,
                    "calculated_drop2_m": calc_drop2,
                    "error2_percent": error2,
                }
            });
            println!("{}", serde_json::to_string_pretty(&result)?);
        },
        
        OutputFormat::Csv => {
            println!("metric,value");
            println!("estimated_bc,{:.4}", estimated_bc);
            println!("error_at_{}m_percent,{:.2}", distance1, error1);
            println!("error_at_{}m_percent,{:.2}", distance2, error2);
        },
        
        OutputFormat::Table => {
            println!("╔════════════════════════════════════════╗");
            println!("║         BC ESTIMATION RESULT           ║");
            println!("╠════════════════════════════════════════╣");
            println!("║ Estimated BC:      {:>8.4}            ║", estimated_bc);
            println!("╠════════════════════════════════════════╣");
            println!("║ Verification:                          ║");
            println!("║ At {:.0}m:                             ║", distance1);
            println!("║   Actual drop:     {:>8.3} m          ║", drop1);
            println!("║   Calculated:      {:>8.3} m          ║", calc_drop1);
            println!("║   Error:           {:>8.2} %          ║", error1);
            println!("║ At {:.0}m:                             ║", distance2);
            println!("║   Actual drop:     {:>8.3} m          ║", drop2);
            println!("║   Calculated:      {:>8.3} m          ║", calc_drop2);
            println!("║   Error:           {:>8.2} %          ║", error2);
            println!("╚════════════════════════════════════════╝");
        },
    }
    
    Ok(())
}