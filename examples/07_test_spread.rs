use ballistics_engine::{
    BallisticInputs, MonteCarloParams, WindConditions,
    run_monte_carlo_with_wind, DragModel,
};

fn main() {
    // Set up base inputs
    let mut base_inputs = BallisticInputs::default();
    base_inputs.muzzle_velocity = 800.0;
    base_inputs.launch_angle = 0.785398;  // 45 degrees
    base_inputs.ballistic_coefficient = 0.5;
    base_inputs.mass = 0.01;
    base_inputs.diameter = 0.008;
    base_inputs.drag_model = DragModel::G1;
    base_inputs.azimuth_angle = 0.0;  // Shooting straight ahead
    
    // Set up wind
    let base_wind = WindConditions {
        speed: 0.0,
        direction: 0.0,
    };
    
    // Set up Monte Carlo parameters with both elevation and azimuth variation
    let mc_params = MonteCarloParams {
        num_simulations: 50,
        velocity_std_dev: 2.0,
        angle_std_dev: 0.002,  // ~0.11 degrees
        bc_std_dev: 0.005,
        wind_speed_std_dev: 0.5,
        target_distance: None,
        base_wind_speed: 0.0,
        base_wind_direction: 0.0,
        azimuth_std_dev: 0.002,  // ~0.11 degrees horizontal spread
    };
    
    // Run simulation
    match run_monte_carlo_with_wind(base_inputs, base_wind, mc_params) {
        Ok(results) => {
            println!("Simulation completed with {} shots", results.ranges.len());
            
            // Extract x and z coordinates from impact positions
            let mut x_coords = Vec::new();
            let mut z_coords = Vec::new();
            
            for pos in &results.impact_positions {
                x_coords.push(pos[0]);
                z_coords.push(pos[2]);
            }
            
            // Calculate statistics
            let x_min = x_coords.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
            let x_max = x_coords.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
            let z_min = z_coords.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
            let z_max = z_coords.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
            
            println!("\nImpact spread:");
            println!("X (range) spread: {:.2} to {:.2} m (spread: {:.2} m)", x_min, x_max, x_max - x_min);
            println!("Z (lateral) spread: {:.2} to {:.2} m (spread: {:.2} m)", z_min, z_max, z_max - z_min);
            
            // Simple ASCII visualization
            println!("\nImpact pattern (top-down view):");
            
            let grid_width = 40;
            let grid_height = 20;
            
            let x_range = x_max - x_min;
            let z_range = z_max - z_min;
            
            if x_range > 0.001 && z_range > 0.001 {
                let mut grid = vec![vec![' '; grid_width]; grid_height];
                
                for (x, z) in x_coords.iter().zip(z_coords.iter()) {
                    let grid_x = ((x - x_min) / x_range * (grid_width as f64 - 1.0)) as usize;
                    let grid_z = ((z - z_min) / z_range * (grid_height as f64 - 1.0)) as usize;
                    
                    if grid_x < grid_width && grid_z < grid_height {
                        grid[grid_z][grid_x] = '*';
                    }
                }
                
                // Print grid with border
                println!("{}", format!("{}{}{}",  "+", "-".repeat(grid_width), "+"));
                for row in grid {
                    print!("|");
                    for ch in row {
                        print!("{}", ch);
                    }
                    println!("|");
                }
                println!("{}", format!("{}{}{}", "+", "-".repeat(grid_width), "+"));
            } else if z_range < 0.001 {
                println!("WARNING: No horizontal spread detected!");
                println!("All shots landed in a vertical line.");
                println!("This suggests azimuth variation is not working correctly.");
            }
            
            // Show some individual impact positions
            println!("\nFirst 5 impact positions (x, y, z):");
            for (i, pos) in results.impact_positions.iter().take(5).enumerate() {
                println!("  Shot {}: ({:.2}, {:.2}, {:.2}) m", i+1, pos[0], pos[1], pos[2]);
            }
        },
        Err(e) => {
            eprintln!("Monte Carlo simulation failed: {}", e);
        }
    }
}