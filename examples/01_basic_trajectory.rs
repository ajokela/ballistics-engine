/// Basic Trajectory Calculation Example
/// 
/// This example demonstrates how to calculate a simple ballistic trajectory
/// using the ballistics engine with basic parameters.

use std::f64::consts::PI;

fn main() {
    println!("=== Basic Trajectory Example ===\n");
    
    // Define projectile parameters
    let velocity = 800.0;      // m/s
    let angle_deg = 45.0;      // degrees
    let mass = 0.02;           // kg (20 grams)
    let diameter = 0.00762;    // meters (7.62mm)
    let bc = 0.5;              // ballistic coefficient
    
    println!("Initial Parameters:");
    println!("  Velocity: {} m/s", velocity);
    println!("  Angle: {}°", angle_deg);
    println!("  Mass: {} kg", mass);
    println!("  Diameter: {} m", diameter);
    println!("  BC: {}", bc);
    println!();
    
    // Calculate trajectory
    let trajectory = calculate_trajectory(velocity, angle_deg, mass, diameter, bc);
    
    // Display results
    println!("Trajectory Results:");
    println!("  Max Range: {:.2} m", trajectory.max_range);
    println!("  Max Height: {:.2} m", trajectory.max_height);
    println!("  Time of Flight: {:.2} s", trajectory.time_of_flight);
    println!("  Impact Velocity: {:.2} m/s", trajectory.impact_velocity);
    println!("  Impact Energy: {:.2} J", trajectory.impact_energy);
    println!();
    
    // Show key trajectory points
    println!("Key Trajectory Points:");
    println!("  Time (s) |   X (m)  |   Y (m)  | Velocity (m/s)");
    println!("  ---------|----------|----------|---------------");
    
    for point in trajectory.key_points {
        println!("  {:8.2} | {:8.2} | {:8.2} | {:8.2}", 
            point.time, point.x, point.y, point.velocity);
    }
}

struct TrajectoryResult {
    max_range: f64,
    max_height: f64,
    time_of_flight: f64,
    impact_velocity: f64,
    impact_energy: f64,
    key_points: Vec<TrajectoryPoint>,
}

struct TrajectoryPoint {
    time: f64,
    x: f64,
    y: f64,
    velocity: f64,
}

fn calculate_trajectory(
    velocity: f64, 
    angle_deg: f64, 
    mass: f64, 
    diameter: f64, 
    bc: f64
) -> TrajectoryResult {
    let angle_rad = angle_deg * PI / 180.0;
    let dt = 0.01; // time step
    
    // Initial conditions
    let mut t = 0.0;
    let mut x = 0.0;
    let mut y = 0.0;
    let mut vx = velocity * angle_rad.cos();
    let mut vy = velocity * angle_rad.sin();
    
    let mut max_height = 0.0;
    let mut key_points = Vec::new();
    let point_interval = 1.0; // Record point every second
    let mut next_point_time = 0.0;
    
    // Physics constants
    let g = 9.80665;
    let air_density = 1.225;
    let area = PI * (diameter / 2.0).powi(2);
    let cd = 0.47 / bc; // Simplified drag coefficient
    
    // Simulate trajectory
    while y >= 0.0 && t < 100.0 {
        let v = (vx * vx + vy * vy).sqrt();
        
        // Record key points
        if t >= next_point_time {
            key_points.push(TrajectoryPoint {
                time: t,
                x,
                y,
                velocity: v,
            });
            next_point_time += point_interval;
        }
        
        // Track max height
        if y > max_height {
            max_height = y;
        }
        
        // Calculate drag
        let drag = 0.5 * air_density * cd * area * v;
        let ax = -drag * vx / mass;
        let ay = -drag * vy / mass - g;
        
        // Update state (Euler integration)
        vx += ax * dt;
        vy += ay * dt;
        x += vx * dt;
        y += vy * dt;
        t += dt;
    }
    
    // Final values
    let impact_velocity = (vx * vx + vy * vy).sqrt();
    let impact_energy = 0.5 * mass * impact_velocity * impact_velocity;
    
    // Add final point
    key_points.push(TrajectoryPoint {
        time: t,
        x,
        y: 0.0,
        velocity: impact_velocity,
    });
    
    TrajectoryResult {
        max_range: x,
        max_height,
        time_of_flight: t,
        impact_velocity,
        impact_energy,
        key_points,
    }
}