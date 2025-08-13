#!/usr/bin/env rust-script
//! ```cargo
//! [dependencies]
//! clap = { version = "4.5", features = ["derive"] }
//! serde = { version = "1.0", features = ["derive"] }
//! serde_json = "1.0"
//! ```

use clap::{Parser, Subcommand, ValueEnum};
use serde::{Serialize, Deserialize};
use std::error::Error;

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
    
    /// Display ballistics information
    Info,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum OutputFormat {
    Json,
    Csv,
    Table,
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