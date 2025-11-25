use magnus::{Error, RHash, function};
use ballistics_engine::{
    DragModel, BallisticInputs, WindConditions, AtmosphericConditions, TrajectorySolver,
};

// Unit conversion constants
const GRAINS_TO_KG: f64 = 0.00006479891;
const FPS_TO_MPS: f64 = 0.3048;
const YARDS_TO_METERS: f64 = 0.9144;
const INCHES_TO_METERS: f64 = 0.0254;
const MPH_TO_MPS: f64 = 0.44704;
const DEGREES_TO_RADIANS: f64 = std::f64::consts::PI / 180.0;

/// Calculate trajectory from Ruby hash input
fn solve_trajectory(ruby: &magnus::Ruby, inputs_hash: RHash) -> Result<RHash, Error> {
    // Extract required values from Ruby hash
    let bc: f64 = inputs_hash.fetch("bc")?;
    let bullet_weight_grains: f64 = inputs_hash.fetch("bullet_weight_grains")?;
    let muzzle_velocity_fps: f64 = inputs_hash.fetch("muzzle_velocity_fps")?;
    let bullet_diameter_inches: f64 = inputs_hash.fetch("bullet_diameter_inches")?;
    let bullet_length_inches: f64 = inputs_hash.fetch("bullet_length_inches")?;
    let sight_height_inches: f64 = inputs_hash.fetch("sight_height_inches")?;
    let zero_distance_yards: f64 = inputs_hash.fetch("zero_distance_yards")?;

    // Optional values with defaults
    let shooting_angle_degrees: f64 = inputs_hash.lookup2("shooting_angle_degrees", 0.0)?;
    let twist_rate_inches: f64 = inputs_hash.lookup2("twist_rate_inches", 10.0)?;
    let is_right_twist: bool = inputs_hash.lookup2("is_right_twist", true)?;

    // Drag model (default to G7)
    let drag_model_str: String = inputs_hash.lookup2("drag_model", "G7")?;
    let drag_model = match drag_model_str.to_uppercase().as_str() {
        "G1" => DragModel::G1,
        "G7" => DragModel::G7,
        "G8" => DragModel::G8,
        _ => return Err(Error::new(ruby.exception_arg_error(), "Invalid drag_model, must be G1, G7, or G8")),
    };

    // Create ballistic inputs using defaults and override specific fields
    let ballistic_inputs = BallisticInputs {
        bc_type: drag_model,
        bc_value: bc,
        bullet_diameter: bullet_diameter_inches * INCHES_TO_METERS,
        bullet_mass: bullet_weight_grains * GRAINS_TO_KG,
        bullet_length: bullet_length_inches * INCHES_TO_METERS,
        muzzle_velocity: muzzle_velocity_fps * FPS_TO_MPS,
        sight_height: sight_height_inches * INCHES_TO_METERS,
        target_distance: zero_distance_yards * YARDS_TO_METERS,
        shooting_angle: shooting_angle_degrees * DEGREES_TO_RADIANS,
        twist_rate: twist_rate_inches,  // Already in inches per the struct definition
        is_twist_right: is_right_twist,
        caliber_inches: bullet_diameter_inches,
        weight_grains: bullet_weight_grains,
        ..Default::default()
    };

    // Optional wind conditions (default to no wind)
    let wind = if let Some(wind_hash) = inputs_hash.lookup::<_, Option<RHash>>("wind")? {
        let speed_mph: f64 = wind_hash.lookup2("speed_mph", 0.0)?;
        let direction_deg: f64 = wind_hash.lookup2("direction_degrees", 0.0)?;

        WindConditions {
            speed: speed_mph * MPH_TO_MPS,
            direction: direction_deg * DEGREES_TO_RADIANS,
        }
    } else {
        WindConditions {
            speed: 0.0,
            direction: 0.0,
        }
    };

    // Optional atmospheric conditions (default to standard conditions)
    let atmosphere = if let Some(atm_hash) = inputs_hash.lookup::<_, Option<RHash>>("atmosphere")? {
        let temp_f: f64 = atm_hash.lookup2("temperature_f", 59.0)?;
        let pressure_inhg: f64 = atm_hash.lookup2("pressure_inhg", 29.92)?;
        let humidity: f64 = atm_hash.lookup2("humidity_percent", 50.0)?;
        let altitude_ft: f64 = atm_hash.lookup2("altitude_feet", 0.0)?;

        // Convert to Celsius and other SI units
        let temp_c = (temp_f - 32.0) * 5.0 / 9.0;
        let pressure_pa = pressure_inhg * 3386.389;
        let altitude_m = altitude_ft * 0.3048;

        AtmosphericConditions {
            temperature: temp_c,
            pressure: pressure_pa,
            humidity,
            altitude: altitude_m,
        }
    } else {
        // Standard ICAO atmosphere
        AtmosphericConditions {
            temperature: 15.0,  // 15°C (59°F)
            pressure: 101325.0, // 1 atm in Pa
            humidity: 50.0,
            altitude: 0.0,
        }
    };

    // Solve trajectory - handle Result properly
    let solver = TrajectorySolver::new(ballistic_inputs, wind, atmosphere);
    let result = solver.solve()
        .map_err(|e| Error::new(ruby.exception_runtime_error(), e.to_string()))?;

    // Convert result to Ruby hash using Ruby context methods
    let result_hash = ruby.hash_new();

    result_hash.aset("max_range_yards", result.max_range / YARDS_TO_METERS)?;
    result_hash.aset("max_height_yards", result.max_height / YARDS_TO_METERS)?;
    result_hash.aset("time_of_flight", result.time_of_flight)?;
    result_hash.aset("impact_velocity_fps", result.impact_velocity / FPS_TO_MPS)?;
    result_hash.aset("impact_energy_ftlbs", result.impact_energy * 0.737562)?; // J to ft-lbs

    // Convert trajectory points to array of hashes - use correct field name "points"
    let points = ruby.ary_new();
    for point in result.points {
        let point_hash = ruby.hash_new();
        point_hash.aset("time", point.time)?;
        point_hash.aset("x", point.position.x / YARDS_TO_METERS)?;
        point_hash.aset("y", point.position.y / YARDS_TO_METERS)?;
        point_hash.aset("z", point.position.z / YARDS_TO_METERS)?;

        // Use the velocity_magnitude field directly
        point_hash.aset("velocity_fps", point.velocity_magnitude / FPS_TO_MPS)?;

        // Use the kinetic_energy field directly
        point_hash.aset("energy_ftlbs", point.kinetic_energy * 0.737562)?;

        points.push(point_hash)?;
    }

    result_hash.aset("points", points)?;

    Ok(result_hash)
}

#[magnus::init]
fn init(ruby: &magnus::Ruby) -> Result<(), Error> {
    let module = ruby.define_module("BallisticsEngine")?;
    module.define_module_function("solve", function!(solve_trajectory, 1))?;
    Ok(())
}
