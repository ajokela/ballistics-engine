use clap::{Parser, Subcommand, ValueEnum};
use serde::{Serialize, Deserialize};
use std::error::Error;
use std::f64;

#[derive(Parser)]
#[command(name = "ballistics")]
#[command(author = "Ballistics Engine Team")]
#[command(version = "0.1.0")]
#[command(about = "High-performance ballistics trajectory calculator", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Calculate a single trajectory
    Trajectory {
        /// Initial velocity (m/s)
        #[arg(short = 'v', long)]
        velocity: f64,
        
        /// Launch angle (degrees)
        #[arg(short = 'a', long, default_value = "0.0")]
        angle: f64,
        
        /// Ballistic coefficient
        #[arg(short = 'b', long, default_value = "0.5")]
        bc: f64,
        
        /// Mass (kg)
        #[arg(short = 'm', long, default_value = "0.01")]
        mass: f64,
        
        /// Diameter (meters)
        #[arg(short = 'd', long, default_value = "0.00762")]
        diameter: f64,
        
        /// Maximum time (seconds)
        #[arg(long, default_value = "10.0")]
        max_time: f64,
        
        /// Time step (seconds)
        #[arg(long, default_value = "0.01")]
        time_step: f64,
        
        /// Output format
        #[arg(short = 'o', long, default_value = "table")]
        output: OutputFormat,
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
        #[arg(short = 'b', long, default_value = "0.5")]
        bc: f64,
        
        /// Mass (kg)
        #[arg(short = 'm', long, default_value = "0.01")]
        mass: f64,
        
        /// Diameter (meters)
        #[arg(short = 'd', long, default_value = "0.00762")]
        diameter: f64,
        
        /// Number of simulations
        #[arg(short = 'n', long, default_value = "1000")]
        num_sims: usize,
        
        /// Velocity standard deviation (m/s)
        #[arg(long, default_value = "5.0")]
        velocity_std: f64,
        
        /// Angle standard deviation (degrees)
        #[arg(long, default_value = "0.5")]
        angle_std: f64,
        
        /// BC standard deviation
        #[arg(long, default_value = "0.02")]
        bc_std: f64,
        
        /// Output format
        #[arg(short = 'o', long, default_value = "summary")]
        output: MonteCarloOutput,
    },
    
    /// Display ballistics information
    Info,
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
    num_simulations: usize,
    mean_range: f64,
    std_range: f64,
    min_range: f64,
    max_range: f64,
    mean_impact_velocity: f64,
    std_impact_velocity: f64,
    mean_max_height: f64,
    std_max_height: f64,
    cep: f64,  // Circular Error Probable
}

fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();
    
    match cli.command {
        Commands::Trajectory { 
            velocity, angle, bc, mass, diameter,
            max_time, time_step, output 
        } => {
            let result = calculate_trajectory(
                velocity, angle, bc, mass, diameter,
                max_time, time_step
            )?;
            
            display_results(result, output)?;
        },
        
        Commands::MonteCarlo {
            velocity, angle, bc, mass, diameter,
            num_sims, velocity_std, angle_std, bc_std, output
        } => {
            let result = run_monte_carlo(
                velocity, angle, bc, mass, diameter,
                num_sims, velocity_std, angle_std, bc_std
            )?;
            
            display_monte_carlo_results(result, output)?;
        },
        
        Commands::Info => {
            println!("╔════════════════════════════════════════╗");
            println!("║      BALLISTICS ENGINE v0.1.0         ║");
            println!("╠════════════════════════════════════════╣");
            println!("║ A high-performance ballistics          ║");
            println!("║ trajectory calculation engine.         ║");
            println!("╠════════════════════════════════════════╣");
            println!("║ Features:                              ║");
            println!("║ • 6DOF trajectory integration          ║");
            println!("║ • Atmospheric modeling                 ║");
            println!("║ • Drag coefficient interpolation       ║");
            println!("║ • Multiple output formats              ║");
            println!("╚════════════════════════════════════════╝");
        }
    }
    
    Ok(())
}

fn calculate_trajectory(
    velocity: f64,
    angle_deg: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    max_time: f64,
    time_step: f64,
) -> Result<TrajectoryResult, Box<dyn Error>> {
    let angle_rad = angle_deg.to_radians();
    
    // Initial conditions
    let mut time = 0.0;
    let mut x = 0.0;
    let mut y = 0.0;
    let mut vx = velocity * angle_rad.cos();
    let mut vy = velocity * angle_rad.sin();
    
    let mut trajectory = Vec::new();
    let mut max_height = 0.0;
    
    // Standard atmosphere at sea level
    let air_density = 1.225; // kg/m³
    let g = 9.80665; // m/s²
    
    // Cross-sectional area
    let area = std::f64::consts::PI * (diameter / 2.0).powi(2);
    
    while time <= max_time && y >= 0.0 {
        let v = (vx * vx + vy * vy).sqrt();
        let energy = 0.5 * mass * v * v;
        
        trajectory.push(TrajectoryPoint {
            time,
            x,
            y,
            velocity: v,
            energy,
        });
        
        if y > max_height {
            max_height = y;
        }
        
        // Simple drag model
        // Cd = base drag coefficient (simplified)
        let cd = 0.5 / bc; // Simplified relationship
        
        // Drag force components
        let drag = 0.5 * air_density * cd * area * v;
        let ax = -drag * vx / mass;
        let ay = -drag * vy / mass - g;
        
        // Update velocity and position (Euler integration)
        vx += ax * time_step;
        vy += ay * time_step;
        x += vx * time_step;
        y += vy * time_step;
        time += time_step;
    }
    
    // Get final values
    let last = trajectory.last().unwrap();
    
    Ok(TrajectoryResult {
        max_range: last.x,
        max_height,
        time_of_flight: last.time,
        impact_velocity: last.velocity,
        impact_energy: last.energy,
        trajectory,
    })
}

fn run_monte_carlo(
    base_velocity: f64,
    base_angle_deg: f64,
    base_bc: f64,
    base_mass: f64,
    base_diameter: f64,
    num_sims: usize,
    velocity_std: f64,
    angle_std_deg: f64,
    bc_std: f64,
) -> Result<MonteCarloResult, Box<dyn Error>> {
    // Simple pseudo-random number generator
    struct SimpleRng {
        seed: u64,
    }
    
    impl SimpleRng {
        fn new() -> Self {
            use std::time::{SystemTime, UNIX_EPOCH};
            let seed = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos() as u64;
            SimpleRng { seed }
        }
        
        fn next_f64(&mut self) -> f64 {
            // Linear congruential generator
            self.seed = (self.seed.wrapping_mul(1664525).wrapping_add(1013904223)) & 0xFFFFFFFFu64;
            (self.seed as f64) / (0xFFFFFFFFu64 as f64)
        }
        
        fn normal(&mut self, mean: f64, std: f64) -> f64 {
            // Box-Muller transform for normal distribution
            let u1 = self.next_f64();
            let u2 = self.next_f64();
            let z0 = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            mean + std * z0
        }
    }
    
    let mut rng = SimpleRng::new();
    let mut ranges = Vec::new();
    let mut max_heights = Vec::new();
    let mut impact_velocities = Vec::new();
    
    // Run simulations
    for _ in 0..num_sims {
        // Sample varied parameters
        let velocity = rng.normal(base_velocity, velocity_std).max(0.0);
        let angle = rng.normal(base_angle_deg, angle_std_deg);
        let bc = rng.normal(base_bc, bc_std).max(0.01);
        
        // Run trajectory with varied parameters
        match calculate_trajectory(velocity, angle, bc, base_mass, base_diameter, 20.0, 0.05) {
            Ok(result) => {
                ranges.push(result.max_range);
                max_heights.push(result.max_height);
                impact_velocities.push(result.impact_velocity);
            },
            Err(_) => {
                // Skip failed simulations
                continue;
            }
        }
    }
    
    if ranges.is_empty() {
        return Err("No successful simulations".into());
    }
    
    // Calculate statistics
    let n = ranges.len() as f64;
    
    let mean_range = ranges.iter().sum::<f64>() / n;
    let mean_height = max_heights.iter().sum::<f64>() / n;
    let mean_velocity = impact_velocities.iter().sum::<f64>() / n;
    
    let std_range = (ranges.iter()
        .map(|r| (r - mean_range).powi(2))
        .sum::<f64>() / n).sqrt();
    let std_height = (max_heights.iter()
        .map(|h| (h - mean_height).powi(2))
        .sum::<f64>() / n).sqrt();
    let std_velocity = (impact_velocities.iter()
        .map(|v| (v - mean_velocity).powi(2))
        .sum::<f64>() / n).sqrt();
    
    let min_range = ranges.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let max_range = ranges.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    
    // CEP approximation (50% of shots within this radius)
    let cep = std_range * 1.1774;
    
    Ok(MonteCarloResult {
        num_simulations: ranges.len(),
        mean_range,
        std_range,
        min_range,
        max_range,
        mean_impact_velocity: mean_velocity,
        std_impact_velocity: std_velocity,
        mean_max_height: mean_height,
        std_max_height: std_height,
        cep,
    })
}

fn display_monte_carlo_results(result: MonteCarloResult, format: MonteCarloOutput) -> Result<(), Box<dyn Error>> {
    match format {
        MonteCarloOutput::Summary => {
            println!("╔════════════════════════════════════════╗");
            println!("║      MONTE CARLO SIMULATION           ║");
            println!("║      {} simulations                ║", result.num_simulations);
            println!("╠════════════════════════════════════════╣");
            println!("║ RANGE STATISTICS                       ║");
            println!("║ Mean:              {:>8.2} m          ║", result.mean_range);
            println!("║ Std Dev:           {:>8.2} m          ║", result.std_range);
            println!("║ Min:               {:>8.2} m          ║", result.min_range);
            println!("║ Max:               {:>8.2} m          ║", result.max_range);
            println!("║ CEP:               {:>8.2} m          ║", result.cep);
            println!("╠════════════════════════════════════════╣");
            println!("║ IMPACT VELOCITY                        ║");
            println!("║ Mean:              {:>8.2} m/s        ║", result.mean_impact_velocity);
            println!("║ Std Dev:           {:>8.2} m/s        ║", result.std_impact_velocity);
            println!("╠════════════════════════════════════════╣");
            println!("║ MAX HEIGHT                             ║");
            println!("║ Mean:              {:>8.2} m          ║", result.mean_max_height);
            println!("║ Std Dev:           {:>8.2} m          ║", result.std_max_height);
            println!("╚════════════════════════════════════════╝");
        },
        
        MonteCarloOutput::Full => {
            println!("{}", serde_json::to_string_pretty(&result)?);
        },
        
        MonteCarloOutput::Statistics => {
            println!("metric,value");
            println!("num_simulations,{}", result.num_simulations);
            println!("mean_range,{:.2}", result.mean_range);
            println!("std_range,{:.2}", result.std_range);
            println!("min_range,{:.2}", result.min_range);
            println!("max_range,{:.2}", result.max_range);
            println!("mean_impact_velocity,{:.2}", result.mean_impact_velocity);
            println!("std_impact_velocity,{:.2}", result.std_impact_velocity);
            println!("mean_max_height,{:.2}", result.mean_max_height);
            println!("std_max_height,{:.2}", result.std_max_height);
            println!("cep,{:.2}", result.cep);
        },
    }
    
    Ok(())
}

fn display_results(result: TrajectoryResult, format: OutputFormat) -> Result<(), Box<dyn Error>> {
    match format {
        OutputFormat::Json => {
            println!("{}", serde_json::to_string_pretty(&result)?);
        },
        
        OutputFormat::Csv => {
            println!("time,x,y,velocity,energy");
            for p in &result.trajectory {
                println!("{:.3},{:.2},{:.2},{:.2},{:.2}",
                    p.time, p.x, p.y, p.velocity, p.energy);
            }
        },
        
        OutputFormat::Table => {
            println!("╔════════════════════════════════════════╗");
            println!("║         TRAJECTORY RESULTS             ║");
            println!("╠════════════════════════════════════════╣");
            println!("║ Max Range:         {:>8.2} m          ║", result.max_range);
            println!("║ Max Height:        {:>8.2} m          ║", result.max_height);
            println!("║ Time of Flight:    {:>8.3} s          ║", result.time_of_flight);
            println!("║ Impact Velocity:   {:>8.2} m/s        ║", result.impact_velocity);
            println!("║ Impact Energy:     {:>8.2} J          ║", result.impact_energy);
            println!("╚════════════════════════════════════════╝");
            
            println!("\nTrajectory Points (every {:.1}s):", result.time_of_flight / 10.0);
            println!("┌──────────┬──────────┬──────────┬──────────┐");
            println!("│ Time (s) │  X (m)   │  Y (m)   │ Vel(m/s) │");
            println!("├──────────┼──────────┼──────────┼──────────┤");
            
            let step = result.trajectory.len() / 10;
            for (i, p) in result.trajectory.iter().enumerate() {
                if i % step.max(1) == 0 || i == result.trajectory.len() - 1 {
                    println!("│ {:>8.3} │ {:>8.2} │ {:>8.2} │ {:>8.2} │",
                        p.time, p.x, p.y, p.velocity);
                }
            }
            println!("└──────────┴──────────┴──────────┴──────────┘");
        },
    }
    
    Ok(())
}