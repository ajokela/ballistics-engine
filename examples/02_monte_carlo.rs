/// Monte Carlo Simulation Example
/// 
/// This example demonstrates how to run Monte Carlo simulations to analyze
/// the statistical distribution of trajectories with parameter uncertainties.

use std::f64::consts::PI;

fn main() {
    println!("=== Monte Carlo Simulation Example ===\n");
    
    // Base parameters
    let base_velocity = 850.0;     // m/s
    let base_angle = 35.0;         // degrees
    let base_bc = 0.45;            // ballistic coefficient
    let mass = 0.0324;             // kg (500 grains)
    let diameter = 0.0127;         // meters (0.50 caliber)
    
    // Uncertainties (standard deviations)
    let velocity_std = 5.0;        // m/s
    let angle_std = 0.5;           // degrees
    let bc_std = 0.02;             // BC variation
    
    let num_simulations = 1000;
    
    println!("Base Parameters:");
    println!("  Velocity: {} ± {} m/s", base_velocity, velocity_std);
    println!("  Angle: {} ± {}°", base_angle, angle_std);
    println!("  BC: {} ± {}", base_bc, bc_std);
    println!("  Simulations: {}", num_simulations);
    println!();
    
    // Run Monte Carlo simulation
    let results = run_monte_carlo(
        base_velocity, base_angle, base_bc, mass, diameter,
        velocity_std, angle_std, bc_std, num_simulations
    );
    
    // Display results
    println!("Monte Carlo Results ({} successful runs):", results.successful_runs);
    println!();
    println!("Range Statistics:");
    println!("  Mean: {:.2} m", results.mean_range);
    println!("  Std Dev: {:.2} m", results.std_range);
    println!("  Min: {:.2} m", results.min_range);
    println!("  Max: {:.2} m", results.max_range);
    println!("  95% Confidence: {:.2} - {:.2} m", 
        results.mean_range - 2.0 * results.std_range,
        results.mean_range + 2.0 * results.std_range);
    println!();
    
    println!("Impact Velocity:");
    println!("  Mean: {:.2} m/s", results.mean_impact_velocity);
    println!("  Std Dev: {:.2} m/s", results.std_impact_velocity);
    println!();
    
    println!("Maximum Height:");
    println!("  Mean: {:.2} m", results.mean_max_height);
    println!("  Std Dev: {:.2} m", results.std_max_height);
    println!();
    
    println!("Dispersion Analysis:");
    println!("  CEP (50%): {:.2} m", results.cep);
    println!("  R95 (95%): {:.2} m", results.r95);
    
    // Show distribution histogram
    println!("\nRange Distribution:");
    show_histogram(&results.ranges, 10);
}

struct MonteCarloResults {
    successful_runs: usize,
    ranges: Vec<f64>,
    mean_range: f64,
    std_range: f64,
    min_range: f64,
    max_range: f64,
    mean_impact_velocity: f64,
    std_impact_velocity: f64,
    mean_max_height: f64,
    std_max_height: f64,
    cep: f64,
    r95: f64,
}

fn run_monte_carlo(
    base_velocity: f64,
    base_angle: f64,
    base_bc: f64,
    mass: f64,
    diameter: f64,
    velocity_std: f64,
    angle_std: f64,
    bc_std: f64,
    num_simulations: usize,
) -> MonteCarloResults {
    let mut rng = SimpleRng::new();
    let mut ranges = Vec::new();
    let mut impact_velocities = Vec::new();
    let mut max_heights = Vec::new();
    
    for _ in 0..num_simulations {
        // Generate random parameters with normal distribution
        let velocity = rng.normal(base_velocity, velocity_std).max(0.0);
        let angle = rng.normal(base_angle, angle_std);
        let bc = rng.normal(base_bc, bc_std).max(0.01);
        
        // Run trajectory simulation
        let result = simulate_trajectory(velocity, angle, bc, mass, diameter);
        
        ranges.push(result.range);
        impact_velocities.push(result.impact_velocity);
        max_heights.push(result.max_height);
    }
    
    // Calculate statistics
    let n = ranges.len() as f64;
    
    let mean_range = ranges.iter().sum::<f64>() / n;
    let std_range = (ranges.iter()
        .map(|r| (r - mean_range).powi(2))
        .sum::<f64>() / n).sqrt();
    
    let mean_velocity = impact_velocities.iter().sum::<f64>() / n;
    let std_velocity = (impact_velocities.iter()
        .map(|v| (v - mean_velocity).powi(2))
        .sum::<f64>() / n).sqrt();
    
    let mean_height = max_heights.iter().sum::<f64>() / n;
    let std_height = (max_heights.iter()
        .map(|h| (h - mean_height).powi(2))
        .sum::<f64>() / n).sqrt();
    
    let min_range = ranges.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let max_range = ranges.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    
    // CEP and R95 calculations
    let cep = std_range * 1.1774;  // 50% circular error
    let r95 = std_range * 2.4477;  // 95% circular error
    
    MonteCarloResults {
        successful_runs: ranges.len(),
        ranges: ranges.clone(),
        mean_range,
        std_range,
        min_range,
        max_range,
        mean_impact_velocity: mean_velocity,
        std_impact_velocity: std_velocity,
        mean_max_height: mean_height,
        std_max_height: std_height,
        cep,
        r95,
    }
}

struct TrajectoryResult {
    range: f64,
    max_height: f64,
    impact_velocity: f64,
}

fn simulate_trajectory(
    velocity: f64,
    angle_deg: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
) -> TrajectoryResult {
    let angle_rad = angle_deg * PI / 180.0;
    let dt = 0.01;
    
    let mut x = 0.0;
    let mut y = 0.0;
    let mut vx = velocity * angle_rad.cos();
    let mut vy = velocity * angle_rad.sin();
    let mut max_height = 0.0;
    
    let g = 9.80665;
    let air_density = 1.225;
    let area = PI * (diameter / 2.0).powi(2);
    let cd = 0.47 / bc;
    
    while y >= 0.0 {
        if y > max_height {
            max_height = y;
        }
        
        let v = (vx * vx + vy * vy).sqrt();
        let drag = 0.5 * air_density * cd * area * v;
        let ax = -drag * vx / mass;
        let ay = -drag * vy / mass - g;
        
        vx += ax * dt;
        vy += ay * dt;
        x += vx * dt;
        y += vy * dt;
    }
    
    let impact_velocity = (vx * vx + vy * vy).sqrt();
    
    TrajectoryResult {
        range: x,
        max_height,
        impact_velocity,
    }
}

// Simple random number generator
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
        self.seed = (self.seed.wrapping_mul(1664525).wrapping_add(1013904223)) & 0xFFFFFFFFu64;
        (self.seed as f64) / (0xFFFFFFFFu64 as f64)
    }
    
    fn normal(&mut self, mean: f64, std: f64) -> f64 {
        let u1 = self.next_f64();
        let u2 = self.next_f64();
        let z0 = (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();
        mean + std * z0
    }
}

fn show_histogram(data: &[f64], bins: usize) {
    let min = data.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let max = data.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let range = max - min;
    let bin_width = range / bins as f64;
    
    let mut histogram = vec![0; bins];
    
    for &value in data {
        let bin = ((value - min) / bin_width).floor() as usize;
        let bin = bin.min(bins - 1);
        histogram[bin] += 1;
    }
    
    let max_count = *histogram.iter().max().unwrap();
    let scale = 40.0 / max_count as f64;
    
    for i in 0..bins {
        let range_start = min + i as f64 * bin_width;
        let range_end = range_start + bin_width;
        let bar_length = (histogram[i] as f64 * scale) as usize;
        let bar = "█".repeat(bar_length);
        
        println!("  {:6.0}-{:6.0} m: {} ({})", 
            range_start, range_end, bar, histogram[i]);
    }
}