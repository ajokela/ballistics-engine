/// Wind Effects Example
/// 
/// This example demonstrates how wind conditions affect ballistic trajectories,
/// including crosswind, headwind, and tailwind scenarios.

use std::f64::consts::PI;

fn main() {
    println!("=== Wind Effects on Trajectory ===\n");
    
    // Projectile parameters
    let velocity = 850.0;      // m/s
    let angle_deg = 30.0;      // degrees
    let mass = 0.0194;         // kg (300 grains)
    let diameter = 0.00762;    // meters (7.62mm / .308)
    let bc = 0.45;
    
    println!("Projectile Parameters:");
    println!("  Velocity: {} m/s", velocity);
    println!("  Angle: {}°", angle_deg);
    println!("  Mass: {} kg ({:.0} grains)", mass, mass * 15432.358);
    println!("  Caliber: {:.2} mm", diameter * 1000.0);
    println!("  BC: {}", bc);
    println!();
    
    // Test different wind conditions
    let wind_scenarios = vec![
        ("No Wind", 0.0, 0.0),
        ("10 m/s Headwind", 10.0, 0.0),
        ("10 m/s Tailwind", 10.0, 180.0),
        ("10 m/s Crosswind (Right)", 10.0, 90.0),
        ("10 m/s Crosswind (Left)", 10.0, 270.0),
        ("15 m/s Quartering (45°)", 15.0, 45.0),
    ];
    
    println!("Wind Scenario Analysis:");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!(" Scenario              | Range (m) | Drift (m) | Impact V (m/s) | TOF (s)");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    
    let mut baseline_range = 0.0;
    
    for (name, wind_speed, wind_dir) in wind_scenarios {
        let result = calculate_with_wind(
            velocity, angle_deg, mass, diameter, bc,
            wind_speed, wind_dir
        );
        
        if name == "No Wind" {
            baseline_range = result.range;
        }
        
        println!(" {:20} | {:9.2} | {:9.2} | {:14.2} | {:7.3}",
            name, result.range, result.lateral_drift, 
            result.impact_velocity, result.time_of_flight);
    }
    
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!();
    
    // Detailed crosswind analysis
    println!("Crosswind Drift Analysis (10 m/s crosswind):");
    println!("  Distance (m) | Drift (m) | Drift (MOA)");
    println!("  -------------|-----------|------------");
    
    let distances = vec![100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0];
    
    for target_distance in distances {
        let drift = calculate_drift_at_distance(
            velocity, angle_deg, mass, diameter, bc,
            10.0, 90.0, target_distance
        );
        
        // Convert drift to MOA (Minutes of Angle)
        // 1 MOA ≈ 1.047 inches at 100 yards ≈ 0.291 mrad
        let moa = (drift / target_distance) * 3437.75; // Convert to MOA
        
        println!("  {:12.0} | {:9.3} | {:10.2}", target_distance, drift, moa);
    }
    
    println!();
    println!("Wind Correction Tips:");
    println!("  • Headwind: Increases trajectory arc, reduces range");
    println!("  • Tailwind: Flattens trajectory, increases range");
    println!("  • Crosswind: Causes lateral drift proportional to time of flight");
    println!("  • Wind drift increases non-linearly with distance");
}

struct TrajectoryResult {
    range: f64,
    lateral_drift: f64,
    impact_velocity: f64,
    time_of_flight: f64,
}

fn calculate_with_wind(
    velocity: f64,
    angle_deg: f64,
    mass: f64,
    diameter: f64,
    bc: f64,
    wind_speed: f64,
    wind_direction_deg: f64,
) -> TrajectoryResult {
    let angle_rad = angle_deg * PI / 180.0;
    let wind_dir_rad = wind_direction_deg * PI / 180.0;
    
    // Wind components (assuming projectile travels in +X direction)
    let wind_x = -wind_speed * wind_dir_rad.cos(); // Headwind is negative
    let wind_z = wind_speed * wind_dir_rad.sin();   // Crosswind
    
    let dt = 0.001; // Fine time step for accuracy
    
    // Initial conditions
    let mut t = 0.0;
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0; // Lateral position
    let mut vx = velocity * angle_rad.cos();
    let mut vy = velocity * angle_rad.sin();
    let mut vz = 0.0; // Lateral velocity
    
    // Physics constants
    let g = 9.80665;
    let air_density = 1.225;
    let area = PI * (diameter / 2.0).powi(2);
    let cd = 0.47 / bc;
    
    // Simulate trajectory
    while y >= 0.0 && t < 100.0 {
        // Velocity relative to air (including wind)
        let vx_rel = vx - wind_x;
        let vy_rel = vy; // No vertical wind
        let vz_rel = vz - wind_z;
        
        let v_rel = (vx_rel * vx_rel + vy_rel * vy_rel + vz_rel * vz_rel).sqrt();
        
        // Drag force components
        let drag_factor = 0.5 * air_density * cd * area * v_rel / mass;
        let ax = -drag_factor * vx_rel;
        let ay = -drag_factor * vy_rel - g;
        let az = -drag_factor * vz_rel;
        
        // Update state
        vx += ax * dt;
        vy += ay * dt;
        vz += az * dt;
        
        x += vx * dt;
        y += vy * dt;
        z += vz * dt;
        
        t += dt;
    }
    
    let impact_velocity = (vx * vx + vy * vy + vz * vz).sqrt();
    
    TrajectoryResult {
        range: x,
        lateral_drift: z,
        impact_velocity,
        time_of_flight: t,
    }
}

fn calculate_drift_at_distance(
    velocity: f64,
    angle_deg: f64,
    mass: f64,
    diameter: f64,
    bc: f64,
    wind_speed: f64,
    wind_direction_deg: f64,
    target_distance: f64,
) -> f64 {
    let angle_rad = angle_deg * PI / 180.0;
    let wind_dir_rad = wind_direction_deg * PI / 180.0;
    
    let wind_x = -wind_speed * wind_dir_rad.cos();
    let wind_z = wind_speed * wind_dir_rad.sin();
    
    let dt = 0.001;
    
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    let mut vx = velocity * angle_rad.cos();
    let mut vy = velocity * angle_rad.sin();
    let mut vz = 0.0;
    
    let g = 9.80665;
    let air_density = 1.225;
    let area = PI * (diameter / 2.0).powi(2);
    let cd = 0.47 / bc;
    
    // Simulate until target distance or ground impact
    while x < target_distance && y >= 0.0 {
        let vx_rel = vx - wind_x;
        let vy_rel = vy;
        let vz_rel = vz - wind_z;
        
        let v_rel = (vx_rel * vx_rel + vy_rel * vy_rel + vz_rel * vz_rel).sqrt();
        
        let drag_factor = 0.5 * air_density * cd * area * v_rel / mass;
        let ax = -drag_factor * vx_rel;
        let ay = -drag_factor * vy_rel - g;
        let az = -drag_factor * vz_rel;
        
        vx += ax * dt;
        vy += ay * dt;
        vz += az * dt;
        
        x += vx * dt;
        y += vy * dt;
        z += vz * dt;
    }
    
    z.abs()
}