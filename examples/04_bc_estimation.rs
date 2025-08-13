/// Ballistic Coefficient Estimation Example
/// 
/// This example demonstrates how to estimate the ballistic coefficient (BC)
/// of a projectile from observed trajectory data points.

use std::f64::consts::PI;

fn main() {
    println!("=== Ballistic Coefficient Estimation ===\n");
    
    // Known projectile parameters
    let velocity = 762.0;       // m/s (2500 fps)
    let mass = 0.01166;         // kg (180 grains)
    let diameter = 0.00762;     // meters (.308 caliber)
    
    println!("Known Parameters:");
    println!("  Muzzle Velocity: {} m/s ({:.0} fps)", velocity, velocity / 0.3048);
    println!("  Mass: {} kg ({:.0} grains)", mass, mass * 15432.358);
    println!("  Diameter: {} mm", diameter * 1000.0);
    println!();
    
    // Observed trajectory data (distance, drop) in meters
    // These might come from doppler radar, high-speed camera, or range testing
    let observed_data = vec![
        (100.0, 0.038),   // 100m, 38mm drop
        (200.0, 0.178),   // 200m, 178mm drop
        (300.0, 0.445),   // 300m, 445mm drop
        (400.0, 0.871),   // 400m, 871mm drop
        (500.0, 1.493),   // 500m, 1493mm drop
        (600.0, 2.352),   // 600m, 2352mm drop
    ];
    
    println!("Observed Trajectory Data:");
    println!("  Distance (m) | Drop (m)");
    println!("  -------------|----------");
    for (dist, drop) in &observed_data {
        println!("  {:12.0} | {:8.3}", dist, drop);
    }
    println!();
    
    // Method 1: Simple iterative estimation
    println!("Method 1: Iterative BC Estimation");
    let estimated_bc = estimate_bc_iterative(&observed_data, velocity, mass, diameter);
    println!("  Estimated BC: {:.4}", estimated_bc);
    
    // Verify the estimation
    println!("\n  Verification:");
    println!("  Distance | Observed | Calculated | Error");
    println!("  ---------|----------|------------|-------");
    
    let mut total_error = 0.0;
    for (dist, observed_drop) in &observed_data {
        let calculated_drop = calculate_drop_at_distance(
            velocity, mass, diameter, estimated_bc, *dist
        );
        let error = (calculated_drop - observed_drop).abs();
        total_error += error;
        
        println!("  {:8.0} | {:8.3} | {:10.3} | {:6.3}",
            dist, observed_drop, calculated_drop, error);
    }
    
    let avg_error = total_error / observed_data.len() as f64;
    println!("  Average error: {:.3} m", avg_error);
    println!();
    
    // Method 2: Least squares estimation
    println!("Method 2: Least Squares BC Estimation");
    let ls_bc = estimate_bc_least_squares(&observed_data, velocity, mass, diameter);
    println!("  Estimated BC: {:.4}", ls_bc);
    
    // Compare different BC values
    println!("\nBC Sensitivity Analysis:");
    println!("  BC    | Avg Error (m) | Max Error (m)");
    println!("  ------|---------------|---------------");
    
    let bc_values = vec![0.35, 0.40, 0.45, 0.50, 0.55, 0.60];
    for bc in bc_values {
        let (avg_err, max_err) = calculate_bc_error(&observed_data, velocity, mass, diameter, bc);
        let marker = if (bc - estimated_bc).abs() < 0.01 { " <--" } else { "" };
        println!("  {:.2} | {:13.4} | {:13.4}{}",
            bc, avg_err, max_err, marker);
    }
    
    // Show how BC affects trajectory
    println!("\nTrajectory Comparison at Different BCs:");
    let distances = vec![200.0, 400.0, 600.0, 800.0, 1000.0];
    
    println!("  Distance | BC=0.40 | BC={:.2} | BC=0.60", estimated_bc);
    println!("  ---------|---------|---------|--------");
    
    for dist in distances {
        let drop_low = calculate_drop_at_distance(velocity, mass, diameter, 0.40, dist);
        let drop_est = calculate_drop_at_distance(velocity, mass, diameter, estimated_bc, dist);
        let drop_high = calculate_drop_at_distance(velocity, mass, diameter, 0.60, dist);
        
        println!("  {:8.0} | {:7.2} | {:7.2} | {:7.2}",
            dist, drop_low, drop_est, drop_high);
    }
}

fn estimate_bc_iterative(
    observed_data: &[(f64, f64)],
    velocity: f64,
    mass: f64,
    diameter: f64,
) -> f64 {
    let mut best_bc = 0.5;
    let mut best_error = f64::MAX;
    
    // Try different BC values
    for bc_test in (100..1000).step_by(5) {
        let bc = bc_test as f64 / 1000.0;
        
        let mut total_error = 0.0;
        for (dist, observed_drop) in observed_data {
            let calculated_drop = calculate_drop_at_distance(
                velocity, mass, diameter, bc, *dist
            );
            let error = (calculated_drop - observed_drop).powi(2);
            total_error += error;
        }
        
        if total_error < best_error {
            best_error = total_error;
            best_bc = bc;
        }
    }
    
    // Refine with smaller steps
    let start = ((best_bc - 0.05) * 10000.0) as i32;
    let end = ((best_bc + 0.05) * 10000.0) as i32;
    
    for bc_test in start..end {
        let bc = bc_test as f64 / 10000.0;
        
        let mut total_error = 0.0;
        for (dist, observed_drop) in observed_data {
            let calculated_drop = calculate_drop_at_distance(
                velocity, mass, diameter, bc, *dist
            );
            let error = (calculated_drop - observed_drop).powi(2);
            total_error += error;
        }
        
        if total_error < best_error {
            best_error = total_error;
            best_bc = bc;
        }
    }
    
    best_bc
}

fn estimate_bc_least_squares(
    observed_data: &[(f64, f64)],
    velocity: f64,
    mass: f64,
    diameter: f64,
) -> f64 {
    // Use gradient descent to minimize squared error
    let mut bc = 0.5;
    let learning_rate = 0.0001;
    let iterations = 1000;
    
    for _ in 0..iterations {
        let mut gradient = 0.0;
        
        for (dist, observed_drop) in observed_data {
            let calculated_drop = calculate_drop_at_distance(
                velocity, mass, diameter, bc, *dist
            );
            let error = calculated_drop - observed_drop;
            
            // Numerical gradient
            let bc_delta = 0.001;
            let drop_plus = calculate_drop_at_distance(
                velocity, mass, diameter, bc + bc_delta, *dist
            );
            let derivative = (drop_plus - calculated_drop) / bc_delta;
            
            gradient += 2.0 * error * derivative;
        }
        
        bc -= learning_rate * gradient;
        bc = bc.max(0.1).min(1.5); // Constrain to reasonable range
    }
    
    bc
}

fn calculate_drop_at_distance(
    velocity: f64,
    mass: f64,
    diameter: f64,
    bc: f64,
    target_distance: f64,
) -> f64 {
    let angle_rad: f64 = 0.0; // Assume zero angle (flat fire)
    let dt = 0.001;
    
    let mut x = 0.0;
    let mut y = 0.0;
    let mut vx = velocity * angle_rad.cos();
    let mut vy = velocity * angle_rad.sin();
    
    let g = 9.80665;
    let air_density = 1.225;
    let area = PI * (diameter / 2.0).powi(2);
    let cd = 0.47 / bc;
    
    // Store the height at target distance
    let mut drop_at_target = 0.0;
    let mut prev_x = 0.0;
    let mut prev_y = 0.0;
    
    while x < target_distance * 1.1 && y > -100.0 {
        let v = (vx * vx + vy * vy).sqrt();
        let drag = 0.5 * air_density * cd * area * v;
        
        let ax = -drag * vx / mass;
        let ay = -drag * vy / mass - g;
        
        // Check if we crossed the target distance
        if prev_x <= target_distance && x > target_distance {
            // Linear interpolation
            let t = (target_distance - prev_x) / (x - prev_x);
            let y_at_target = prev_y + t * (y - prev_y);
            drop_at_target = -y_at_target; // Drop is negative of height
            break;
        }
        
        prev_x = x;
        prev_y = y;
        
        vx += ax * dt;
        vy += ay * dt;
        x += vx * dt;
        y += vy * dt;
    }
    
    drop_at_target
}

fn calculate_bc_error(
    observed_data: &[(f64, f64)],
    velocity: f64,
    mass: f64,
    diameter: f64,
    bc: f64,
) -> (f64, f64) {
    let mut total_error = 0.0;
    let mut max_error: f64 = 0.0;
    
    for (dist, observed_drop) in observed_data {
        let calculated_drop = calculate_drop_at_distance(
            velocity, mass, diameter, bc, *dist
        );
        let error = (calculated_drop - observed_drop).abs();
        total_error += error;
        max_error = max_error.max(error);
    }
    
    let avg_error = total_error / observed_data.len() as f64;
    (avg_error, max_error)
}