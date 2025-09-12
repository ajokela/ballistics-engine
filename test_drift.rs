use ballistics_engine::*;

fn main() {
    // Test parameters
    let muzzle_velocity = 800.0; // m/s
    let bc = 0.400;
    let mass_grains = 120.0;
    let diameter_inches = 0.284;
    let zero_range = 100.0; // meters
    let wind_speed = 10.0; // m/s
    let wind_direction = 89.0; // degrees
    
    let mut inputs = BallisticInputs::default();
    inputs.sight_height = 0.05; // 50mm above bore
    inputs.muzzle_height = 1.0; // meters above ground
    inputs.target_height = 1.0; // meters above ground
    inputs.bullet_mass = mass_grains;
    inputs.bullet_diameter = diameter_inches;
    inputs.ballistic_coefficient = bc;
    inputs.muzzle_velocity = muzzle_velocity;
    inputs.zero_range = Some(zero_range);
    inputs.enable_advanced_effects = true;
    inputs.twist_rate = 10.0;
    inputs.is_twist_right = true;
    inputs.latitude = 45.29;
    
    // Wind vector (89 degrees = almost pure crosswind)
    // x = speed * sin(89°), z = speed * cos(89°)
    let wind_rad = wind_direction.to_radians();
    inputs.wind_speed = wind_speed;
    inputs.wind_direction = wind_rad;
    
    // Calculate trajectory
    let mut solver = TrajectoryCalculator::new(inputs);
    solver.set_max_range(200.0);
    
    let trajectory = solver.calculate_trajectory();
    
    println!("Wind Vector:");
    println!("  Direction: {} degrees", wind_direction);
    println!("  Speed: {} m/s", wind_speed);
    println!("  X component: {:.3} m/s", wind_speed * wind_rad.sin());
    println!("  Z component: {:.3} m/s", wind_speed * wind_rad.cos());
    println!();
    
    println!("Range (m) | Drift (cm) | Height (m) | Time (s)");
    println!("----------|------------|------------|----------");
    
    // Sample points
    let ranges = vec![0.0, 50.0, 100.0, 142.0, 150.0, 200.0];
    
    for target_range in ranges {
        if let Some(point) = trajectory.points.iter().find(|p| p.position.z >= target_range) {
            println!("{:8.1} | {:10.1} | {:10.3} | {:8.3}",
                point.position.z,
                point.position.x * 100.0, // Convert to cm
                point.position.y,
                point.time
            );
        }
    }
    
    println!("\nMax range: {:.1} m", trajectory.max_range);
    println!("Ground strike at: {:.1} m", 
        trajectory.points.iter()
            .find(|p| p.position.y <= 0.0)
            .map(|p| p.position.z)
            .unwrap_or(trajectory.max_range)
    );
}