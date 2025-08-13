// Test program to demonstrate the units conversion system

#[derive(Debug, Clone, Copy)]
enum UnitSystem {
    Metric,
    Imperial,
}

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
}

fn main() {
    println!("Testing Unit Conversion System\n");
    
    // Test Imperial to Metric
    println!("Imperial Input -> Metric Processing:");
    let velocity_fps = 2800.0;
    let mass_grains = 180.0;
    let distance_yards = 1000.0;
    
    let velocity_mps = UnitConverter::velocity_to_metric(velocity_fps, UnitSystem::Imperial);
    let mass_kg = UnitConverter::mass_to_metric(mass_grains, UnitSystem::Imperial);
    let distance_m = UnitConverter::distance_to_metric(distance_yards, UnitSystem::Imperial);
    
    println!("  {} fps -> {:.2} m/s", velocity_fps, velocity_mps);
    println!("  {} grains -> {:.6} kg", mass_grains, mass_kg);
    println!("  {} yards -> {:.2} meters", distance_yards, distance_m);
    
    println!("\nMetric Processing -> Imperial Output:");
    let output_velocity = UnitConverter::velocity_from_metric(velocity_mps, UnitSystem::Imperial);
    let output_distance = UnitConverter::distance_from_metric(distance_m, UnitSystem::Imperial);
    
    println!("  {:.2} m/s -> {:.2} fps", velocity_mps, output_velocity);
    println!("  {:.2} meters -> {:.2} yards", distance_m, output_distance);
    
    println!("\nMetric In/Out (no conversion):");
    let metric_velocity = 850.0;
    let metric_mass = 0.01;
    let metric_distance = 900.0;
    
    let processed_velocity = UnitConverter::velocity_to_metric(metric_velocity, UnitSystem::Metric);
    let output_velocity_metric = UnitConverter::velocity_from_metric(processed_velocity, UnitSystem::Metric);
    
    println!("  {} m/s -> {} m/s (unchanged)", metric_velocity, output_velocity_metric);
}