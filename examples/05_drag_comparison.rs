/// Drag Model Comparison Example
/// 
/// This example compares different projectile types and their trajectories,
/// demonstrating the differences between various ballistic coefficients and masses.

use std::f64::consts::PI;

fn main() {
    println!("=== Projectile Comparison ===\n");
    
    // Define different projectile types
    let projectiles = vec![
        Projectile {
            name: ".22 LR (40gr)",
            velocity: 375.0,      // m/s
            mass: 0.00259,        // kg (40 grains)
            diameter: 0.00564,    // meters (.22 caliber)
            bc_g1: 0.138,
        },
        Projectile {
            name: "5.56 NATO (62gr)",
            velocity: 948.0,      // m/s
            mass: 0.00402,        // kg (62 grains)
            diameter: 0.00564,    // meters (.224 caliber)
            bc_g1: 0.307,
        },
        Projectile {
            name: ".308 Win (168gr)",
            velocity: 850.0,      // m/s
            mass: 0.01089,        // kg (168 grains)
            diameter: 0.00762,    // meters (.308 caliber)
            bc_g1: 0.462,
        },
        Projectile {
            name: ".338 Lapua (250gr)",
            velocity: 915.0,      // m/s
            mass: 0.01620,        // kg (250 grains)
            diameter: 0.00861,    // meters (.338 caliber)
            bc_g1: 0.675,
        },
        Projectile {
            name: ".50 BMG (750gr)",
            velocity: 860.0,      // m/s
            mass: 0.04860,        // kg (750 grains)
            diameter: 0.01270,    // meters (.50 caliber)
            bc_g1: 1.050,
        },
    ];
    
    // Compare at specific distances
    let check_distances = vec![100.0, 300.0, 500.0, 800.0, 1000.0];
    
    println!("Drop Comparison (meters):");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    print!(" Projectile          |");
    for dist in &check_distances {
        print!(" {:>6.0}m |", dist);
    }
    println!();
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    
    for proj in &projectiles {
        print!(" {:19} |", proj.name);
        for dist in &check_distances {
            let drop = calculate_drop(proj, *dist);
            if drop < 999.0 {
                print!(" {:>7.2} |", drop);
            } else {
                print!("    ---  |");
            }
        }
        println!();
    }
    
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!();
    
    // Time of flight comparison
    println!("Time of Flight (seconds):");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    print!(" Projectile          |");
    for dist in &check_distances {
        print!(" {:>6.0}m |", dist);
    }
    println!();
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    
    for proj in &projectiles {
        print!(" {:19} |", proj.name);
        for dist in &check_distances {
            let tof = calculate_time_to_distance(proj, *dist);
            if tof < 99.0 {
                print!(" {:>7.3} |", tof);
            } else {
                print!("    ---  |");
            }
        }
        println!();
    }
    
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!();
    
    // Remaining velocity comparison
    println!("Remaining Velocity (m/s):");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    print!(" Projectile          |");
    for dist in &check_distances {
        print!(" {:>6.0}m |", dist);
    }
    println!();
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    
    for proj in &projectiles {
        print!(" {:19} |", proj.name);
        for dist in &check_distances {
            let vel = calculate_velocity_at_distance(proj, *dist);
            if vel > 0.0 {
                print!(" {:>7.0} |", vel);
            } else {
                print!("    ---  |");
            }
        }
        println!();
    }
    
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!();
    
    // Energy comparison
    println!("Kinetic Energy (Joules):");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    print!(" Projectile          | Muzzle  |");
    for dist in &[100.0, 300.0, 500.0] {
        print!(" {:>6.0}m |", dist);
    }
    println!();
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    
    for proj in &projectiles {
        let muzzle_energy = 0.5 * proj.mass * proj.velocity * proj.velocity;
        print!(" {:19} | {:>7.0} |", proj.name, muzzle_energy);
        
        for dist in &[100.0, 300.0, 500.0] {
            let vel = calculate_velocity_at_distance(proj, *dist);
            if vel > 0.0 {
                let energy = 0.5 * proj.mass * vel * vel;
                print!(" {:>7.0} |", energy);
            } else {
                print!("    ---  |");
            }
        }
        println!();
    }
    
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!();
    
    // Maximum range comparison
    println!("Performance Summary:");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!(" Projectile          | Max Range | Supersonic | 50% Energy | Zero@100m");
    println!("                     |    (m)    | Range (m)  | Range (m)  | Drop@200m");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    
    for proj in &projectiles {
        let max_range = calculate_max_range(proj);
        let supersonic_range = calculate_supersonic_range(proj);
        let half_energy_range = calculate_half_energy_range(proj);
        let zero_drop = calculate_zero_drop_at_200(proj);
        
        println!(" {:19} | {:>9.0} | {:>10.0} | {:>10.0} | {:>9.2}",
            proj.name, max_range, supersonic_range, half_energy_range, zero_drop);
    }
    
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
}

struct Projectile {
    name: &'static str,
    velocity: f64,  // m/s
    mass: f64,      // kg
    diameter: f64,  // m
    bc_g1: f64,     // G1 ballistic coefficient
}

fn calculate_drop(proj: &Projectile, target_distance: f64) -> f64 {
    let dt = 0.001;
    let mut x = 0.0;
    let mut y = 0.0;
    let mut vx = proj.velocity;
    let mut vy = 0.0;
    
    let g = 9.80665;
    let air_density = 1.225;
    let area = PI * (proj.diameter / 2.0).powi(2);
    let cd = 0.47 / proj.bc_g1;
    
    while x < target_distance && y > -100.0 {
        let v = (vx * vx + vy * vy).sqrt();
        let drag = 0.5 * air_density * cd * area * v;
        
        let ax = -drag * vx / proj.mass;
        let ay = -drag * vy / proj.mass - g;
        
        vx += ax * dt;
        vy += ay * dt;
        x += vx * dt;
        y += vy * dt;
    }
    
    -y // Return drop (positive value)
}

fn calculate_time_to_distance(proj: &Projectile, target_distance: f64) -> f64 {
    let dt = 0.001;
    let mut t = 0.0;
    let mut x = 0.0;
    let mut y = 0.0;
    let mut vx = proj.velocity;
    let mut vy = 0.0;
    
    let g = 9.80665;
    let air_density = 1.225;
    let area = PI * (proj.diameter / 2.0).powi(2);
    let cd = 0.47 / proj.bc_g1;
    
    while x < target_distance && y > -100.0 && t < 100.0 {
        let v = (vx * vx + vy * vy).sqrt();
        let drag = 0.5 * air_density * cd * area * v;
        
        let ax = -drag * vx / proj.mass;
        let ay = -drag * vy / proj.mass - g;
        
        vx += ax * dt;
        vy += ay * dt;
        x += vx * dt;
        y += vy * dt;
        t += dt;
    }
    
    if x >= target_distance {
        t
    } else {
        999.0 // Didn't reach target
    }
}

fn calculate_velocity_at_distance(proj: &Projectile, target_distance: f64) -> f64 {
    let dt = 0.001;
    let mut x = 0.0;
    let mut y = 0.0;
    let mut vx = proj.velocity;
    let mut vy = 0.0;
    
    let g = 9.80665;
    let air_density = 1.225;
    let area = PI * (proj.diameter / 2.0).powi(2);
    let cd = 0.47 / proj.bc_g1;
    
    while x < target_distance && y > -100.0 {
        let v = (vx * vx + vy * vy).sqrt();
        let drag = 0.5 * air_density * cd * area * v;
        
        let ax = -drag * vx / proj.mass;
        let ay = -drag * vy / proj.mass - g;
        
        vx += ax * dt;
        vy += ay * dt;
        x += vx * dt;
        y += vy * dt;
    }
    
    if x >= target_distance {
        (vx * vx + vy * vy).sqrt()
    } else {
        0.0
    }
}

fn calculate_max_range(proj: &Projectile) -> f64 {
    // Calculate at 45 degree angle for maximum range
    let angle = 45.0 * PI / 180.0;
    let dt = 0.01;
    let mut x = 0.0;
    let mut y = 0.0;
    let mut vx = proj.velocity * angle.cos();
    let mut vy = proj.velocity * angle.sin();
    
    let g = 9.80665;
    let air_density = 1.225;
    let area = PI * (proj.diameter / 2.0).powi(2);
    let cd = 0.47 / proj.bc_g1;
    
    while y >= 0.0 && x < 50000.0 {
        let v = (vx * vx + vy * vy).sqrt();
        let drag = 0.5 * air_density * cd * area * v;
        
        let ax = -drag * vx / proj.mass;
        let ay = -drag * vy / proj.mass - g;
        
        vx += ax * dt;
        vy += ay * dt;
        x += vx * dt;
        y += vy * dt;
    }
    
    x
}

fn calculate_supersonic_range(proj: &Projectile) -> f64 {
    let speed_of_sound = 343.0; // m/s at sea level
    let dt = 0.001;
    let mut x = 0.0;
    let mut vx = proj.velocity;
    let mut vy = 0.0;
    
    let g = 9.80665;
    let air_density = 1.225;
    let area = PI * (proj.diameter / 2.0).powi(2);
    let cd = 0.47 / proj.bc_g1;
    
    while vx > speed_of_sound && x < 10000.0 {
        let v = (vx * vx + vy * vy).sqrt();
        let drag = 0.5 * air_density * cd * area * v;
        
        let ax = -drag * vx / proj.mass;
        let ay = -drag * vy / proj.mass - g;
        
        vx += ax * dt;
        vy += ay * dt;
        x += vx * dt;
    }
    
    x
}

fn calculate_half_energy_range(proj: &Projectile) -> f64 {
    let initial_energy = 0.5 * proj.mass * proj.velocity * proj.velocity;
    let half_energy = initial_energy / 2.0;
    let dt = 0.001;
    let mut x = 0.0;
    let mut vx = proj.velocity;
    let mut vy = 0.0;
    
    let g = 9.80665;
    let air_density = 1.225;
    let area = PI * (proj.diameter / 2.0).powi(2);
    let cd = 0.47 / proj.bc_g1;
    
    loop {
        let v = (vx * vx + vy * vy).sqrt();
        let current_energy = 0.5 * proj.mass * v * v;
        
        if current_energy <= half_energy || x > 10000.0 {
            break;
        }
        
        let drag = 0.5 * air_density * cd * area * v;
        let ax = -drag * vx / proj.mass;
        let ay = -drag * vy / proj.mass - g;
        
        vx += ax * dt;
        vy += ay * dt;
        x += vx * dt;
    }
    
    x
}

fn calculate_zero_drop_at_200(proj: &Projectile) -> f64 {
    // Calculate drop at 200m when zeroed at 100m
    // First find the angle needed for 100m zero
    let mut best_angle = 0.0;
    let mut best_error = f64::MAX;
    
    for angle_test in -100..100 {
        let angle = angle_test as f64 * 0.0001; // Small angle in radians
        let drop = calculate_drop_with_angle(proj, 100.0, angle);
        
        if drop.abs() < best_error {
            best_error = drop.abs();
            best_angle = angle;
        }
    }
    
    // Now calculate drop at 200m with this angle
    calculate_drop_with_angle(proj, 200.0, best_angle)
}

fn calculate_drop_with_angle(proj: &Projectile, target_distance: f64, angle: f64) -> f64 {
    let dt = 0.001;
    let mut x: f64 = 0.0;
    let mut y: f64 = 0.0;
    let mut vx = proj.velocity * angle.cos();
    let mut vy = proj.velocity * angle.sin();
    
    let g = 9.80665;
    let air_density = 1.225;
    let area = PI * (proj.diameter / 2.0).powi(2);
    let cd = 0.47 / proj.bc_g1;
    
    while x < target_distance && y.abs() < 100.0 {
        let v = (vx * vx + vy * vy).sqrt();
        let drag = 0.5 * air_density * cd * area * v;
        
        let ax = -drag * vx / proj.mass;
        let ay = -drag * vy / proj.mass - g;
        
        vx += ax * dt;
        vy += ay * dt;
        x += vx * dt;
        y += vy * dt;
    }
    
    -y
}