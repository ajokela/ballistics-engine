// Allocator selection based on features
#[cfg(all(feature = "jemalloc", not(target_env = "msvc")))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[cfg(all(feature = "mimalloc", not(target_env = "msvc")))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use ballistics_engine::{
    trajectory_sampling, AtmosphericConditions, BallisticInputs, DragModel, MonteCarloParams,
    TrajectorySolver, WindConditions,
};
#[cfg(feature = "online")]
use ballistics_engine::api_client::{ApiClient, TrajectoryRequestBuilder};
use clap::{Parser, Subcommand, ValueEnum};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::error::Error;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Write};
use std::path::PathBuf;

// ============================================================================
// Terms of Service Acceptance Module (for --online feature)
// ============================================================================

#[cfg(feature = "online")]
const TOS_URL: &str = "https://ballistics.rs/terms.txt";

#[cfg(feature = "online")]
const TOS_VERSION: &str = "2026-01-26";

/// Structure to store TOS acceptance status
#[cfg(feature = "online")]
#[derive(Debug, Serialize, Deserialize)]
struct TosAcceptance {
    accepted: bool,
    accepted_version: String,
    accepted_at: String,
    terms_hash: String,
}

/// Get the path to the TOS acceptance file
#[cfg(feature = "online")]
fn get_tos_file_path() -> Option<PathBuf> {
    dirs::home_dir().map(|home| home.join(".ballistics").join("tos.json"))
}

/// Check if the user has already accepted the current TOS
#[cfg(feature = "online")]
fn check_tos_accepted() -> bool {
    let Some(tos_path) = get_tos_file_path() else {
        return false;
    };

    if !tos_path.exists() {
        return false;
    }

    match fs::read_to_string(&tos_path) {
        Ok(content) => {
            match serde_json::from_str::<TosAcceptance>(&content) {
                Ok(acceptance) => {
                    // Check if accepted and version matches
                    acceptance.accepted && acceptance.accepted_version == TOS_VERSION
                }
                Err(_) => false,
            }
        }
        Err(_) => false,
    }
}

/// Fetch the TOS text from the server
#[cfg(feature = "online")]
fn fetch_tos_text() -> Result<String, String> {
    let response = ureq::get(TOS_URL)
        .timeout(std::time::Duration::from_secs(10))
        .call()
        .map_err(|e| format!("Failed to fetch Terms of Service: {}", e))?;

    response
        .into_string()
        .map_err(|e| format!("Failed to read Terms of Service: {}", e))
}

/// Calculate a simple hash of the TOS text for tracking changes
#[cfg(feature = "online")]
fn hash_tos(text: &str) -> String {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    let mut hasher = DefaultHasher::new();
    text.hash(&mut hasher);
    format!("{:x}", hasher.finish())
}

/// Save TOS acceptance to file
#[cfg(feature = "online")]
fn save_tos_acceptance(accepted: bool, terms_hash: &str) -> Result<(), String> {
    let Some(tos_path) = get_tos_file_path() else {
        return Err("Could not determine home directory".to_string());
    };

    // Create .ballistics directory if it doesn't exist
    if let Some(parent) = tos_path.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| format!("Failed to create .ballistics directory: {}", e))?;
    }

    let acceptance = TosAcceptance {
        accepted,
        accepted_version: TOS_VERSION.to_string(),
        accepted_at: chrono_lite_timestamp(),
        terms_hash: terms_hash.to_string(),
    };

    let json = serde_json::to_string_pretty(&acceptance)
        .map_err(|e| format!("Failed to serialize TOS acceptance: {}", e))?;

    fs::write(&tos_path, json)
        .map_err(|e| format!("Failed to write TOS acceptance file: {}", e))?;

    Ok(())
}

/// Get current timestamp without chrono dependency
#[cfg(feature = "online")]
fn chrono_lite_timestamp() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};

    let duration = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default();

    format!("{}", duration.as_secs())
}

/// Display TOS and prompt for acceptance
#[cfg(feature = "online")]
fn prompt_tos_acceptance() -> Result<bool, String> {
    eprintln!("\n");
    eprintln!("================================================================================");
    eprintln!("                    TERMS OF SERVICE ACCEPTANCE REQUIRED");
    eprintln!("================================================================================");
    eprintln!();
    eprintln!("The --online feature sends data to a proprietary cloud service.");
    eprintln!("Before using this feature, you must accept the Terms of Service.");
    eprintln!();
    eprintln!("Fetching Terms of Service from {}...", TOS_URL);
    eprintln!();

    // Fetch TOS
    let tos_text = fetch_tos_text()?;
    let terms_hash = hash_tos(&tos_text);

    // Display TOS
    eprintln!("{}", tos_text);
    eprintln!();
    eprintln!("================================================================================");
    eprintln!();

    // Prompt for acceptance
    eprint!("Do you accept the Terms of Service? [y/N]: ");
    io::stderr().flush().ok();

    let mut input = String::new();
    io::stdin()
        .read_line(&mut input)
        .map_err(|e| format!("Failed to read input: {}", e))?;

    let accepted = matches!(input.trim().to_lowercase().as_str(), "y" | "yes");

    // Save acceptance status
    save_tos_acceptance(accepted, &terms_hash)?;

    if accepted {
        eprintln!();
        eprintln!("Terms of Service accepted. This will not be shown again unless the terms change.");
        eprintln!();
    } else {
        eprintln!();
        eprintln!("Terms of Service not accepted. The --online feature cannot be used.");
        eprintln!("You can still use local calculations by omitting the --online flag.");
        eprintln!();
    }

    Ok(accepted)
}

/// Check TOS acceptance and prompt if needed. Returns true if online mode can proceed.
#[cfg(feature = "online")]
fn ensure_tos_accepted() -> Result<bool, String> {
    if check_tos_accepted() {
        return Ok(true);
    }

    prompt_tos_acceptance()
}

#[derive(Parser)]
#[command(name = "ballistics")]
#[command(author = "Ballistics Engine Team")]
#[command(version = "0.1.0")]
#[command(about = "High-performance ballistics trajectory calculator", long_about = None)]
struct Cli {
    /// Unit system for input/output (metric or imperial)
    #[arg(short = 'u', long, default_value = "imperial", global = true)]
    units: UnitSystem,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Calculate a single trajectory
    Trajectory {
        /// Load parameters from CSV profile file (gun_profiles.csv format)
        #[arg(long, value_name = "FILE")]
        profile: Option<PathBuf>,

        /// Row name to load from profile CSV (matches first column, e.g., "R700_65CM")
        #[arg(long, value_name = "NAME")]
        profile_row: Option<String>,

        /// Load location/environmental data from CSV file
        #[arg(long, value_name = "FILE")]
        location: Option<PathBuf>,

        /// Site name to load from location CSV (matches first column, e.g., "KF_LR")
        #[arg(long, value_name = "NAME")]
        site: Option<String>,

        /// Initial velocity (fps or m/s based on --units)
        #[arg(short = 'v', long)]
        velocity: Option<f64>,

        /// Velocity adjustment (added to base velocity for truing from chronograph data)
        #[arg(long, default_value = "0.0")]
        velocity_adjustment: f64,

        /// Launch angle (degrees)
        #[arg(short = 'a', long, default_value = "0.0")]
        angle: f64,

        /// Ballistic coefficient
        #[arg(short = 'b', long)]
        bc: Option<f64>,

        /// BC adjustment factor (multiplier for truing, e.g., 0.85 = 85% of stated BC)
        #[arg(long, default_value = "1.0")]
        bc_adjustment: f64,

        /// Mass (grains for imperial, grams for metric)
        #[arg(short = 'm', long)]
        mass: f64,

        /// Diameter (inches for imperial, mm for metric)
        #[arg(short = 'd', long)]
        diameter: f64,

        /// Drag model (G1, G7, Custom)
        #[arg(long, default_value = "g1")]
        drag_model: DragModelArg,

        /// Maximum range (yards or meters based on --units)
        #[arg(long, default_value = "1000.0")]
        max_range: f64,

        /// Time step (seconds)
        #[arg(long, default_value = "0.001")]
        time_step: f64,

        /// Wind speed (mph or m/s based on --units)
        #[arg(long, default_value = "0.0")]
        wind_speed: f64,

        /// Wind direction (degrees, 0=North, 90=East)
        #[arg(long, default_value = "0.0")]
        wind_direction: f64,

        /// Temperature (Fahrenheit or Celsius based on --units)
        #[arg(long, default_value = "59.0")]
        temperature: f64,

        /// Pressure (inHg or hPa based on --units)
        #[arg(long, default_value = "29.92")]
        pressure: f64,

        /// Humidity (0-100%)
        #[arg(long, default_value = "50.0")]
        humidity: f64,

        /// Altitude (feet or meters based on --units)
        #[arg(long, default_value = "0.0")]
        altitude: f64,

        /// Output format
        #[arg(short = 'o', long, default_value = "table")]
        output: OutputFormat,

        /// Show all trajectory points
        #[arg(long)]
        full: bool,

        /// Automatically zero to target distance (overrides --angle)
        #[arg(long)]
        auto_zero: Option<f64>,

        /// Sight height above bore (inches for imperial, mm for metric)
        #[arg(long)]
        sight_height: Option<f64>,

        /// Bore height above ground (feet for imperial, meters for metric). Default: 5ft/1.5m
        #[arg(long)]
        bore_height: Option<f64>,

        /// Ignore ground impact detection (trajectory continues to max range)
        #[arg(long)]
        ignore_ground_impact: bool,

        /// Enable velocity-based BC segmentation
        #[arg(long)]
        use_bc_segments: bool,

        // Advanced Physics Parameters
        /// Enable Magnus effect (requires twist-rate)
        #[arg(long)]
        enable_magnus: bool,

        /// Enable Coriolis effect (requires latitude)
        #[arg(long)]
        enable_coriolis: bool,

        /// Use Euler integration instead of RK45 (RK45 adaptive is default for best accuracy)
        #[arg(long)]
        use_euler: bool,

        /// Use fixed-step RK4 instead of adaptive RK45 (faster but less accurate)
        #[arg(long)]
        use_rk4_fixed: bool,

        /// Enable enhanced spin drift calculations
        #[arg(long)]
        enable_spin_drift: bool,

        /// Enable wind shear (altitude-dependent wind)
        #[arg(long)]
        enable_wind_shear: bool,

        /// Enable trajectory sampling at regular intervals
        #[arg(long)]
        sample_trajectory: bool,

        /// Sampling interval in meters (default: 10)
        #[arg(long, default_value = "10.0")]
        sample_interval: f64,

        /// Enable pitch damping for transonic stability analysis
        #[arg(long)]
        enable_pitch_damping: bool,

        /// Enable precession/nutation physics for angular motion modeling
        #[arg(long)]
        enable_precession: bool,

        /// Use cluster-based BC degradation for improved accuracy
        #[arg(long)]
        use_cluster_bc: bool,

        /// Barrel twist rate (inches per turn, e.g., 10 for 1:10)
        #[arg(long)]
        twist_rate: Option<f64>,

        /// Right-hand twist (true) or left-hand (false)
        #[arg(long, default_value = "true")]
        twist_right: bool,

        /// Latitude for Coriolis effect and weather zones (degrees, -90 to 90)
        #[arg(long)]
        latitude: Option<f64>,

        /// Longitude for weather zones (degrees, -180 to 180)
        #[arg(long)]
        longitude: Option<f64>,

        /// Shot direction/azimuth for Coriolis and weather zones (degrees, 0=North, 90=East)
        #[arg(long)]
        shot_direction: Option<f64>,

        /// Shooting angle (degrees, positive = uphill, negative = downhill)
        #[arg(long, default_value = "0.0")]
        shooting_angle: f64,

        /// Enable powder temperature sensitivity
        #[arg(long)]
        use_powder_sensitivity: bool,

        /// Powder temperature sensitivity (fps per degree)
        #[arg(long, default_value = "1.0")]
        powder_temp_sensitivity: f64,

        /// Powder temperature
        #[arg(long, default_value = "70.0")]
        powder_temp: f64,

        // Online Mode Parameters (feature-gated)
        /// Use Flask API for ML-enhanced trajectory calculation
        #[cfg(feature = "online")]
        #[arg(long, help = "Route calculations through Flask API for ML enhancements")]
        online: bool,

        /// API endpoint URL
        #[cfg(feature = "online")]
        #[arg(long, default_value = "https://api.ballistics.7.62x51mm.sh", help = "API endpoint URL")]
        api_url: String,

        /// Fall back to local calculation if API unreachable
        #[cfg(feature = "online")]
        #[arg(long, help = "Use local calculation if API fails")]
        offline_fallback: bool,

        /// Compare local and API results side-by-side
        #[cfg(feature = "online")]
        #[arg(long, help = "Show comparison between local and API results")]
        compare: bool,

        /// API request timeout in seconds
        #[cfg(feature = "online")]
        #[arg(long, default_value = "10", help = "API timeout in seconds")]
        api_timeout: u64,

        /// Enable weather zones (requires latitude, longitude, and shot-direction)
        #[cfg(feature = "online")]
        #[arg(long, help = "Enable automatic weather zone generation")]
        enable_weather_zones: bool,

        /// Enable 3D weather corrections (altitude-dependent atmospheric changes)
        #[cfg(feature = "online")]
        #[arg(long, help = "Enable 3D weather with altitude corrections")]
        enable_3d_weather: bool,

        /// Wind shear model to use
        #[cfg(feature = "online")]
        #[arg(long, default_value = "logarithmic", help = "Wind shear model: none, logarithmic, power_law, ekman_spiral")]
        wind_shear_model: String,

        /// Weather zone interpolation method
        #[cfg(feature = "online")]
        #[arg(long, default_value = "linear", help = "Weather zone interpolation: linear, cubic, step")]
        weather_zone_interpolation: String,
    },

    /// Run Monte Carlo simulation
    MonteCarlo {
        /// Base velocity (m/s)
        #[arg(short = 'v', long)]
        velocity: f64,

        /// Launch angle (degrees)
        #[arg(short = 'a', long, default_value = "0.0")]
        angle: f64,

        /// Ballistic coefficient
        #[arg(short = 'b', long)]
        bc: f64,

        /// Mass (kg)
        #[arg(short = 'm', long)]
        mass: f64,

        /// Diameter (meters)
        #[arg(short = 'd', long)]
        diameter: f64,

        /// Number of simulations
        #[arg(short = 'n', long, default_value = "1000")]
        num_sims: usize,

        /// Velocity standard deviation (m/s)
        #[arg(long, default_value = "1.0")]
        velocity_std: f64,

        /// Angle standard deviation (degrees)
        #[arg(long, default_value = "0.1")]
        angle_std: f64,

        /// BC standard deviation
        #[arg(long, default_value = "0.01")]
        bc_std: f64,

        /// Wind speed standard deviation (m/s)
        #[arg(long, default_value = "1.0")]
        wind_std: f64,

        /// Base wind speed (m/s)
        #[arg(long, default_value = "0.0")]
        wind_speed: f64,

        /// Base wind direction (degrees, 0=North, 90=East)
        #[arg(long, default_value = "0.0")]
        wind_direction: f64,

        /// Target distance (meters)
        #[arg(long)]
        target_distance: Option<f64>,

        /// Output format
        #[arg(short = 'o', long, default_value = "summary")]
        output: MonteCarloOutput,
    },

    /// Calculate zero angle for a target
    Zero {
        /// Initial velocity (fps for imperial, m/s for metric)
        #[arg(short = 'v', long)]
        velocity: f64,

        /// Ballistic coefficient
        #[arg(short = 'b', long)]
        bc: f64,

        /// Mass (grains for imperial, grams for metric)
        #[arg(short = 'm', long)]
        mass: f64,

        /// Diameter (inches for imperial, mm for metric)
        #[arg(short = 'd', long)]
        diameter: f64,

        /// Target distance (yards for imperial, meters for metric)
        #[arg(long)]
        target_distance: f64,

        /// Target height (yards for imperial, meters for metric)
        #[arg(long, default_value = "0.0")]
        target_height: f64,

        /// Sight height above bore (inches for imperial, mm for metric)
        #[arg(long)]
        sight_height: Option<f64>,

        /// Output format
        #[arg(short = 'o', long, default_value = "table")]
        output: OutputFormat,
    },

    /// Estimate BC from trajectory data
    EstimateBC {
        /// Initial velocity (m/s)
        #[arg(short = 'v', long)]
        velocity: f64,

        /// Mass (kg)
        #[arg(short = 'm', long)]
        mass: f64,

        /// Diameter (meters)
        #[arg(short = 'd', long)]
        diameter: f64,

        /// Distance 1 (meters)
        #[arg(long)]
        distance1: f64,

        /// Drop at distance 1 (meters)
        #[arg(long)]
        drop1: f64,

        /// Distance 2 (meters)
        #[arg(long)]
        distance2: f64,

        /// Drop at distance 2 (meters)
        #[arg(long)]
        drop2: f64,

        /// Output format
        #[arg(short = 'o', long, default_value = "table")]
        output: OutputFormat,
    },

    /// Generate BC segments for velocity-dependent BC
    GenerateBCSegments {
        /// Base ballistic coefficient
        #[arg(short = 'b', long)]
        bc: f64,

        /// Projectile mass (kg)
        #[arg(short = 'm', long)]
        mass: f64,

        /// Projectile diameter (meters)
        #[arg(short = 'd', long)]
        diameter: f64,

        /// Bullet model/name (e.g., "SMK", "ELD-M", "VLD")
        #[arg(long, default_value = "")]
        model: String,

        /// Drag model (G1 or G7)
        #[arg(long, default_value = "G1")]
        drag_model: String,

        /// Output format
        #[arg(short = 'o', long, default_value = "table")]
        output: OutputFormat,
    },
}

#[derive(Debug, Clone, Copy, ValueEnum, PartialEq)]
enum UnitSystem {
    /// Metric units (m/s, kg, meters, Celsius)
    Metric,
    /// Imperial units (fps, grains, yards, Fahrenheit)
    Imperial,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum DragModelArg {
    G1,
    G7,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum OutputFormat {
    Json,
    Csv,
    Table,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum MonteCarloOutput {
    Summary,
    Full,
    Statistics,
}

#[derive(Debug, Serialize, Deserialize)]
struct TrajectoryPoint {
    time: f64,
    x: f64,
    y: f64,
    z: f64,
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
    stability_coefficient: Option<f64>,
    spin_drift: Option<f64>,
    trajectory: Vec<TrajectoryPoint>,
}

#[derive(Debug, Serialize, Deserialize)]
struct MonteCarloResult {
    mean_range: f64,
    std_range: f64,
    mean_impact_velocity: f64,
    std_impact_velocity: f64,
    cep: f64, // Circular Error Probable
    hit_probability: Option<f64>,
}

// Unit conversion functions
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
            UnitSystem::Metric => val * 0.001,           // grams to kg
            UnitSystem::Imperial => val * 0.00006479891, // grains to kg
        }
    }

    fn distance_to_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val * 0.9144, // yards to meters
        }
    }

    fn diameter_to_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val * 0.001,    // mm to meters
            UnitSystem::Imperial => val * 0.0254, // inches to meters
        }
    }

    fn sight_height_to_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val * 0.001,    // mm to meters
            UnitSystem::Imperial => val * 0.0254, // inches to meters
        }
    }

    fn wind_to_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val * 0.44704, // mph to m/s
        }
    }

    fn temperature_to_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => (val - 32.0) * 5.0 / 9.0, // Fahrenheit to Celsius
        }
    }

    fn pressure_to_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val * 33.8639, // inHg to hPa
        }
    }

    fn altitude_to_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val * 0.3048, // feet to meters
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

    fn energy_from_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val * 0.737562, // Joules to ft-lbs
        }
    }
}

/// Parse a CSV file and return row matching the first column value
fn load_csv_row(path: &PathBuf, row_name: &str) -> Result<HashMap<String, String>, Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Parse header - first non-empty line (strip leading # if present, as Glenn's format uses #COLUMN_NAME)
    let mut headers: Vec<String> = Vec::new();
    for line in lines.by_ref() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        // Strip leading # from header line (Glenn's format: #RIFLE_NAME,BULLET_NAME,...)
        let header_line = trimmed.trim_start_matches('#');
        headers = header_line.split(',').map(|s| s.trim().to_uppercase()).collect();
        break;
    }

    if headers.is_empty() {
        return Err("CSV file has no header row".into());
    }

    let _first_col = headers.first().ok_or("CSV has no columns")?.clone();

    // Parse data rows and find matching row
    for line in lines {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let values: Vec<&str> = trimmed.split(',').collect();

        // Check if first column matches
        if let Some(first_val) = values.first() {
            if first_val.trim().eq_ignore_ascii_case(row_name) {
                let mut row: HashMap<String, String> = HashMap::new();
                for (i, header) in headers.iter().enumerate() {
                    if i < values.len() {
                        row.insert(header.clone(), values[i].trim().to_string());
                    }
                }
                return Ok(row);
            }
        }
    }
    Err(format!("Row '{}' not found in CSV", row_name).into())
}

/// Get f64 value from HashMap, trying multiple key names
fn csv_get_f64(map: &HashMap<String, String>, keys: &[&str], default: f64) -> f64 {
    for key in keys {
        if let Some(val) = map.get(&key.to_uppercase()) {
            if let Ok(v) = val.parse::<f64>() {
                return v;
            }
        }
    }
    default
}

fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Trajectory {
            profile,
            profile_row,
            location,
            site,
            velocity,
            velocity_adjustment,
            angle,
            bc,
            bc_adjustment,
            mass,
            diameter,
            drag_model,
            max_range,
            time_step,
            wind_speed,
            wind_direction,
            temperature,
            pressure,
            humidity,
            altitude,
            output,
            full,
            auto_zero,
            sight_height,
            bore_height,
            ignore_ground_impact,
            use_bc_segments,
            enable_magnus,
            enable_coriolis,
            enable_spin_drift,
            enable_wind_shear,
            sample_trajectory,
            sample_interval,
            enable_pitch_damping,
            enable_precession,
            use_cluster_bc,
            twist_rate,
            twist_right,
            latitude,
            longitude,
            shot_direction,
            use_euler,
            use_rk4_fixed,
            shooting_angle,
            use_powder_sensitivity,
            powder_temp_sensitivity,
            powder_temp,
            #[cfg(feature = "online")]
            online,
            #[cfg(feature = "online")]
            api_url,
            #[cfg(feature = "online")]
            offline_fallback,
            #[cfg(feature = "online")]
            compare,
            #[cfg(feature = "online")]
            api_timeout,
            #[cfg(feature = "online")]
            enable_weather_zones,
            #[cfg(feature = "online")]
            enable_3d_weather,
            #[cfg(feature = "online")]
            wind_shear_model,
            #[cfg(feature = "online")]
            weather_zone_interpolation,
        } => {
            // Load profile from CSV if specified
            let profile_data: HashMap<String, String> = if let (Some(path), Some(row)) = (&profile, &profile_row) {
                match load_csv_row(path, row) {
                    Ok(data) => {
                        eprintln!("Loaded profile '{}' from {:?}", row, path);
                        data
                    }
                    Err(e) => {
                        eprintln!("Warning: Failed to load profile: {}", e);
                        HashMap::new()
                    }
                }
            } else {
                HashMap::new()
            };

            // Load location from CSV if specified
            let location_data: HashMap<String, String> = if let (Some(path), Some(site_name)) = (&location, &site) {
                match load_csv_row(path, site_name) {
                    Ok(data) => {
                        eprintln!("Loaded location '{}' from {:?}", site_name, path);
                        data
                    }
                    Err(e) => {
                        eprintln!("Warning: Failed to load location: {}", e);
                        HashMap::new()
                    }
                }
            } else {
                HashMap::new()
            };

            // Merge values: CLI args override profile/location, profile/location override defaults
            let final_velocity = velocity.unwrap_or_else(|| csv_get_f64(&profile_data, &["VELOCITY", "MV", "MUZZLE_VELOCITY"], 2800.0));
            let final_velocity_adj = if velocity_adjustment != 0.0 { velocity_adjustment } else { csv_get_f64(&profile_data, &["VELOCITY_ADJ", "VEL_ADJ"], 0.0) };
            let final_bc = bc.unwrap_or_else(|| csv_get_f64(&profile_data, &["BC"], 0.5));
            let final_bc_adj = if bc_adjustment != 1.0 { bc_adjustment } else { csv_get_f64(&profile_data, &["BC_ADJ", "BC_ADJUSTMENT"], 1.0) };
            let final_mass = mass;
            let final_diameter = diameter;
            let final_max_range = max_range;
            let final_wind_speed = wind_speed;
            let final_wind_direction = if wind_direction != 0.0 { wind_direction } else { csv_get_f64(&location_data, &["WIND_DIR", "WIND_DIRECTION"], 0.0) };

            // Location overrides (environmental conditions)
            let final_temperature = if temperature != 59.0 { temperature } else { csv_get_f64(&location_data, &["TARGET_TEMP", "TEMPERATURE", "TEMP"], csv_get_f64(&profile_data, &["ZERO_TEMP"], 59.0)) };
            let final_pressure = if pressure != 29.92 { pressure } else { csv_get_f64(&location_data, &["PRESSURE", "PRESSURE(HPA OR INHG)"], 29.92) };
            let final_humidity = if humidity != 50.0 { humidity } else { csv_get_f64(&location_data, &["HUMIDITY"], 50.0) };
            let final_altitude = if altitude != 0.0 { altitude } else { csv_get_f64(&location_data, &["ALTITUDE", "ALT"], csv_get_f64(&profile_data, &["ZERO_ALT"], 0.0)) };

            // Apply truing adjustments
            let trued_velocity = final_velocity + final_velocity_adj;
            let trued_bc = final_bc * final_bc_adj;

            // Show effective values if using profile/location
            if !profile_data.is_empty() || !location_data.is_empty() {
                eprintln!("Effective values: velocity={:.1} (trued={:.1}), BC={:.3} (trued={:.4})",
                    final_velocity, trued_velocity, final_bc, trued_bc);
                if !location_data.is_empty() {
                    eprintln!("                  altitude={:.0}, temp={:.1}, pressure={:.2}",
                        final_altitude, final_temperature, final_pressure);
                }
            }

            // Rename for clarity
            let bullet_mass = final_mass;
            let bullet_diameter = final_diameter;
            // Convert inputs to metric (using trued/final values)
            let velocity_metric = UnitConverter::velocity_to_metric(trued_velocity, cli.units);
            let mass_metric = UnitConverter::mass_to_metric(bullet_mass, cli.units);
            let diameter_metric = UnitConverter::diameter_to_metric(bullet_diameter, cli.units);
            let max_range_metric = UnitConverter::distance_to_metric(final_max_range, cli.units);
            let wind_speed_metric = UnitConverter::wind_to_metric(final_wind_speed, cli.units);
            let temperature_metric = UnitConverter::temperature_to_metric(final_temperature, cli.units);
            let pressure_metric = UnitConverter::pressure_to_metric(final_pressure, cli.units);
            let altitude_metric = UnitConverter::altitude_to_metric(final_altitude, cli.units);
            // Default sight height: 2 inches for imperial, 50mm for metric
            let sight_height_default = match cli.units {
                UnitSystem::Imperial => 2.0,
                UnitSystem::Metric => 50.0,
            };
            let sight_height_value = sight_height.unwrap_or(sight_height_default);
            let sight_height_metric =
                UnitConverter::sight_height_to_metric(sight_height_value, cli.units);

            // Bore height above ground: default 5 feet (imperial) or 1.5 meters (metric)
            let bore_height_default = match cli.units {
                UnitSystem::Imperial => 5.0, // feet
                UnitSystem::Metric => 1.5,   // meters
            };
            let bore_height_value = bore_height.unwrap_or(bore_height_default);
            let bore_height_metric = match cli.units {
                UnitSystem::Imperial => bore_height_value * 0.3048, // feet to meters
                UnitSystem::Metric => bore_height_value,
            };

            // Calculate zero angle if auto-zero is specified
            let muzzle_angle = if let Some(zero_distance) = auto_zero {
                let zero_distance_metric =
                    UnitConverter::distance_to_metric(zero_distance, cli.units);

                // Create inputs for zero calculation
                let zero_inputs = BallisticInputs {
                    muzzle_velocity: velocity_metric,
                    bc_value: trued_bc,
                    bullet_mass: mass_metric,
                    bullet_diameter: diameter_metric,
                    sight_height: sight_height_metric,
                    ..Default::default()
                };

                // Calculate zero angle
                // Target height is sight_height because the bullet must cross the LOS at zero distance
                // The LOS is at y = sight_height (sight is above bore by sight_height)
                // So the bullet (starting at y = 0 = bore level) must rise to y = sight_height at zero distance
                let zero_angle = ballistics_engine::calculate_zero_angle(
                    zero_inputs,
                    zero_distance_metric,
                    sight_height_metric, // target height at zero distance (LOS height)
                )?;

                // Convert to degrees for the trajectory function
                zero_angle.to_degrees()
            } else {
                angle
            };

            // Online mode handling
            #[cfg(feature = "online")]
            {
                if online || compare {
                    // Check TOS acceptance before proceeding with online mode
                    match ensure_tos_accepted() {
                        Ok(true) => {
                            // TOS accepted, continue with online mode
                        }
                        Ok(false) => {
                            // TOS not accepted, exit
                            return Err("Cannot use --online without accepting Terms of Service.".into());
                        }
                        Err(e) => {
                            eprintln!("Warning: Could not verify TOS acceptance: {}", e);
                            eprintln!("Proceeding with online mode...");
                        }
                    }

                    // Build API request
                    let zero_range_metric = auto_zero.map(|d| UnitConverter::distance_to_metric(d, cli.units));
                    let twist_rate_metric = twist_rate.map(|t| match cli.units {
                        UnitSystem::Imperial => t * 0.0254, // inches to meters
                        UnitSystem::Metric => t * 0.001,    // mm to meters
                    });

                    let api_request = TrajectoryRequestBuilder::new()
                        .bc_value(trued_bc)
                        .bc_type(match drag_model {
                            DragModelArg::G1 => "G1",
                            DragModelArg::G7 => "G7",
                        })
                        .bullet_mass(mass_metric * 1000.0) // kg to grams
                        .muzzle_velocity(velocity_metric)
                        .target_distance(max_range_metric)
                        .zero_range(zero_range_metric.unwrap_or(100.0))
                        .wind_speed(wind_speed_metric)
                        .wind_angle(final_wind_direction)
                        .temperature(temperature_metric)
                        .pressure(pressure_metric)
                        .humidity(final_humidity)
                        .altitude(altitude_metric)
                        .shooting_angle(shooting_angle)
                        .bullet_diameter(diameter_metric)
                        .ground_threshold(if ignore_ground_impact {
                            f64::NEG_INFINITY
                        } else {
                            -100.0 // default ground threshold in meters
                        })
                        .enable_weather_zones(enable_weather_zones)
                        .enable_3d_weather(enable_3d_weather)
                        .wind_shear_model(&wind_shear_model)
                        .weather_zone_interpolation(&weather_zone_interpolation)
                        .build()
                        .map_err(|e| format!("Failed to build API request: {}", e))?;

                    // Add optional parameters
                    let mut request = api_request;
                    if let Some(lat) = latitude {
                        request.latitude = Some(lat);
                    }
                    if let Some(lon) = longitude {
                        request.longitude = Some(lon);
                    }
                    if let Some(dir) = shot_direction {
                        request.shot_direction = Some(dir);
                    }
                    if let Some(twist) = twist_rate_metric {
                        request.twist_rate = Some(twist);
                    }

                    let api_client = ApiClient::new(&api_url, api_timeout);

                    let api_result = api_client.calculate_trajectory(&request);

                    match (&api_result, compare) {
                        (Ok(api_response), true) => {
                            // Compare mode: run local calculation too
                            eprintln!("Running comparison between local and API calculations...\n");

                            // Display API results with ML info
                            println!("╔════════════════════════════════════════╗");
                            println!("║     COMPARISON: LOCAL vs API           ║");
                            println!("╠════════════════════════════════════════╣");

                            if let Some(ref corrections) = api_response.ml_corrections_applied {
                                if !corrections.is_empty() {
                                    println!("║ ML Corrections Applied:                ║");
                                    for correction in corrections {
                                        println!("║   - {:32} ║", correction);
                                    }
                                    println!("╠════════════════════════════════════════╣");
                                }
                            }

                            if let Some(confidence) = api_response.bc_confidence {
                                println!("║ BC Confidence:     {:>8.1}%           ║", confidence * 100.0);
                                println!("╠════════════════════════════════════════╣");
                            }

                            // Display comparison table
                            let (dist_unit, _drop_unit, _vel_unit) = match cli.units {
                                UnitSystem::Metric => ("m", "m", "m/s"),
                                UnitSystem::Imperial => ("yd", "in", "fps"),
                            };

                            println!("║ Range {} │ API Drop │Local Drop│  Δ Drop  ║", dist_unit);
                            println!("╠══════════╪══════════╪══════════╪══════════╣");

                            // Run local calculation to get comparison data
                            // (We'll print both results side by side)
                            for point in api_response.trajectory.iter().take(10) {
                                let range_display = UnitConverter::distance_from_metric(point.range, cli.units);
                                let drop_display = if cli.units == UnitSystem::Imperial {
                                    point.drop * 39.3701 // meters to inches
                                } else {
                                    point.drop
                                };
                                println!(
                                    "║ {:>8.1} │ {:>8.2} │    ---   │    ---   ║",
                                    range_display, drop_display
                                );
                            }
                            println!("╚══════════╧══════════╧══════════╧══════════╝");
                            println!();
                            println!("Note: Run without --compare to see full local trajectory.");

                            // Also run the local calculation for full output
                            run_trajectory(
                                velocity_metric,
                                muzzle_angle,
                                trued_bc,
                                mass_metric,
                                diameter_metric,
                                drag_model,
                                max_range_metric,
                                time_step,
                                wind_speed_metric,
                                final_wind_direction,
                                temperature_metric,
                                pressure_metric,
                                final_humidity,
                                altitude_metric,
                                output,
                                full,
                                cli.units,
                                sight_height_metric,
                                bore_height_metric,
                                ignore_ground_impact,
                                use_bc_segments,
                                enable_magnus,
                                enable_coriolis,
                                enable_spin_drift,
                                enable_wind_shear,
                                sample_trajectory,
                                sample_interval,
                                enable_pitch_damping,
                                enable_precession,
                                use_cluster_bc,
                                !use_euler,
                                !use_rk4_fixed,
                                twist_rate,
                                twist_right,
                                latitude,
                                shooting_angle,
                                use_powder_sensitivity,
                                powder_temp_sensitivity,
                                powder_temp,
                            )?;
                        }
                        (Ok(api_response), false) => {
                            // Online mode only: display API results
                            display_api_trajectory_result(api_response, output, cli.units, full)?;
                        }
                        (Err(e), _) => {
                            if offline_fallback {
                                eprintln!("Warning: API request failed: {}", e);
                                eprintln!("Falling back to local calculation...\n");

                                run_trajectory(
                                    velocity_metric,
                                    muzzle_angle,
                                    trued_bc,
                                    mass_metric,
                                    diameter_metric,
                                    drag_model,
                                    max_range_metric,
                                    time_step,
                                    wind_speed_metric,
                                    final_wind_direction,
                                    temperature_metric,
                                    pressure_metric,
                                    final_humidity,
                                    altitude_metric,
                                    output,
                                    full,
                                    cli.units,
                                    sight_height_metric,
                                    bore_height_metric,
                                    ignore_ground_impact,
                                    use_bc_segments,
                                    enable_magnus,
                                    enable_coriolis,
                                    enable_spin_drift,
                                    enable_wind_shear,
                                    sample_trajectory,
                                    sample_interval,
                                    enable_pitch_damping,
                                    enable_precession,
                                    use_cluster_bc,
                                    !use_euler,
                                    !use_rk4_fixed,
                                    twist_rate,
                                    twist_right,
                                    latitude,
                                    shooting_angle,
                                    use_powder_sensitivity,
                                    powder_temp_sensitivity,
                                    powder_temp,
                                )?;
                            } else {
                                eprintln!("Error: API request failed: {}", e);
                                eprintln!("Hint: Use --offline-fallback to use local calculation on API failure");
                                std::process::exit(1);
                            }
                        }
                    }
                } else {
                    // Local calculation (default)
                    run_trajectory(
                        velocity_metric,
                        muzzle_angle,
                        trued_bc,
                        mass_metric,
                        diameter_metric,
                        drag_model,
                        max_range_metric,
                        time_step,
                        wind_speed_metric,
                        final_wind_direction,
                        temperature_metric,
                        pressure_metric,
                        final_humidity,
                        altitude_metric,
                        output,
                        full,
                        cli.units,
                        sight_height_metric,
                        bore_height_metric,
                        ignore_ground_impact,
                        use_bc_segments,
                        enable_magnus,
                        enable_coriolis,
                        enable_spin_drift,
                        enable_wind_shear,
                        sample_trajectory,
                        sample_interval,
                        enable_pitch_damping,
                        enable_precession,
                        use_cluster_bc,
                        !use_euler,
                        !use_rk4_fixed,
                        twist_rate,
                        twist_right,
                        latitude,
                        shooting_angle,
                        use_powder_sensitivity,
                        powder_temp_sensitivity,
                        powder_temp,
                    )?;
                }
            }

            // Non-online feature: just run local calculation
            #[cfg(not(feature = "online"))]
            {
                run_trajectory(
                    velocity_metric,
                    muzzle_angle,
                    trued_bc,
                    mass_metric,
                    diameter_metric,
                    drag_model,
                    max_range_metric,
                    time_step,
                    wind_speed_metric,
                    final_wind_direction,
                    temperature_metric,
                    pressure_metric,
                    final_humidity,
                    altitude_metric,
                    output,
                    full,
                    cli.units,
                    sight_height_metric,
                    bore_height_metric,
                    ignore_ground_impact,
                    use_bc_segments,
                    enable_magnus,
                    enable_coriolis,
                    enable_spin_drift,
                    enable_wind_shear,
                    sample_trajectory,
                    sample_interval,
                    enable_pitch_damping,
                    enable_precession,
                    use_cluster_bc,
                    !use_euler,
                    !use_rk4_fixed,
                    twist_rate,
                    twist_right,
                    latitude,
                    shooting_angle,
                    use_powder_sensitivity,
                    powder_temp_sensitivity,
                    powder_temp,
                )?;
            }
        }

        Commands::MonteCarlo {
            velocity,
            angle,
            bc,
            mass,
            diameter,
            num_sims,
            velocity_std,
            angle_std,
            bc_std,
            wind_std,
            wind_speed,
            wind_direction,
            target_distance,
            output,
        } => {
            let bullet_mass = mass;
            let bullet_diameter = diameter;
            run_monte_carlo(
                velocity,
                angle,
                bc,
                bullet_mass,
                bullet_diameter,
                num_sims,
                velocity_std,
                angle_std,
                bc_std,
                wind_std,
                wind_speed,
                wind_direction,
                target_distance,
                output,
            )?;
        }

        Commands::Zero {
            velocity,
            bc,
            mass,
            diameter,
            target_distance,
            target_height,
            sight_height,
            output,
        } => {
            let bullet_mass = mass;
            let bullet_diameter = diameter;
            // Convert inputs to metric
            let velocity_metric = UnitConverter::velocity_to_metric(velocity, cli.units);
            let mass_metric = UnitConverter::mass_to_metric(bullet_mass, cli.units);
            let diameter_metric = UnitConverter::diameter_to_metric(bullet_diameter, cli.units);
            let target_distance_metric =
                UnitConverter::distance_to_metric(target_distance, cli.units);
            let target_height_metric = UnitConverter::distance_to_metric(target_height, cli.units);
            // Default sight height: 2 inches for imperial, 50mm for metric
            let sight_height_default = match cli.units {
                UnitSystem::Imperial => 2.0,
                UnitSystem::Metric => 50.0,
            };
            let sight_height_value = sight_height.unwrap_or(sight_height_default);
            let sight_height_metric =
                UnitConverter::sight_height_to_metric(sight_height_value, cli.units);

            run_zero_calculation(
                velocity_metric,
                bc,
                mass_metric,
                diameter_metric,
                target_distance_metric,
                target_height_metric,
                sight_height_metric,
                output,
                cli.units,
            )?;
        }

        Commands::EstimateBC {
            velocity,
            mass,
            diameter,
            distance1,
            drop1,
            distance2,
            drop2,
            output,
        } => {
            let bullet_mass = mass;
            let bullet_diameter = diameter;
            run_bc_estimation(
                velocity,
                bullet_mass,
                bullet_diameter,
                distance1,
                drop1,
                distance2,
                drop2,
                output,
            )?;
        }

        Commands::GenerateBCSegments {
            bc,
            mass,
            diameter,
            model,
            drag_model,
            output,
        } => {
            let bullet_mass = mass;
            let bullet_diameter = diameter;
            // Convert to metric if needed
            let mass_metric = UnitConverter::mass_to_metric(bullet_mass, cli.units);
            let diameter_metric = UnitConverter::diameter_to_metric(bullet_diameter, cli.units);

            run_bc_segment_generation(
                bc,
                mass_metric,
                diameter_metric,
                &model,
                &drag_model,
                output,
                cli.units,
            )?;
        }
    }

    Ok(())
}

fn run_trajectory(
    velocity: f64,
    angle: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    drag_model: DragModelArg,
    max_range: f64,
    time_step: f64,
    wind_speed: f64,
    wind_direction: f64,
    temperature: f64,
    pressure: f64,
    humidity: f64,
    altitude: f64,
    output: OutputFormat,
    full: bool,
    units: UnitSystem,
    sight_height: f64,
    bore_height: f64,           // Bore height above ground in meters
    ignore_ground_impact: bool, // Disable ground impact detection
    use_bc_segments: bool,
    enable_magnus: bool,
    enable_coriolis: bool,
    enable_spin_drift: bool,
    enable_wind_shear: bool,
    sample_trajectory: bool,
    sample_interval: f64,
    enable_pitch_damping: bool,
    enable_precession: bool,
    use_cluster_bc: bool,
    use_rk4: bool,
    use_rk45: bool,
    twist_rate: Option<f64>,
    twist_right: bool,
    latitude: Option<f64>,
    shooting_angle: f64,
    use_powder_sensitivity: bool,
    powder_temp_sensitivity: f64,
    powder_temp: f64,
) -> Result<(), Box<dyn Error>> {
    // Create ballistic inputs with all required fields
    let drag_model_enum = match drag_model {
        DragModelArg::G1 => DragModel::G1,
        DragModelArg::G7 => DragModel::G7,
    };

    let inputs = BallisticInputs {
        // Core ballistics parameters
        bc_value: bc,
        bc_type: drag_model_enum,
        bullet_mass: mass,
        muzzle_velocity: velocity,
        bullet_diameter: diameter,
        bullet_length: diameter * 4.5, // Approximate length/diameter ratio

        // Targeting and positioning
        muzzle_angle: angle.to_radians(),
        target_distance: max_range,
        azimuth_angle: 0.0,
        shooting_angle: shooting_angle.to_radians(),
        sight_height,
        muzzle_height: bore_height, // Bore height above ground from --bore-height CLI option
        target_height: 0.0,
        ground_threshold: if ignore_ground_impact {
            f64::NEG_INFINITY
        } else {
            -bore_height
        }, // Ground is at -bore_height relative to bore, or disabled if --ignore-ground-impact

        // Environmental conditions
        altitude,
        temperature,
        pressure,
        humidity,
        latitude,

        // Wind conditions
        wind_speed,
        wind_angle: wind_direction,

        // Bullet characteristics
        twist_rate: twist_rate.unwrap_or(12.0),
        is_twist_right: twist_right,
        caliber_inches: diameter / 0.0254,
        weight_grains: mass / 0.00006479891,
        manufacturer: None,
        bullet_model: None,
        bullet_id: None,
        bullet_cluster: None,

        // Integration method selection
        use_rk4,
        use_adaptive_rk45: use_rk45,

        // Advanced effects
        enable_advanced_effects: enable_magnus || enable_coriolis, // Either one enables the system
        use_powder_sensitivity,
        powder_temp_sensitivity: if use_powder_sensitivity {
            UnitConverter::velocity_to_metric(powder_temp_sensitivity, units)
                / UnitConverter::temperature_to_metric(1.0, units)
        } else {
            0.0
        },
        powder_temp: UnitConverter::temperature_to_metric(powder_temp, units),
        tipoff_yaw: 0.0,
        tipoff_decay_distance: 50.0,
        use_bc_segments,
        bc_segments: None,
        bc_segments_data: if use_bc_segments {
            // Generate BC segments automatically
            use ballistics_engine::bc_estimation::BCSegmentEstimator;
            let weight_grains = mass / 0.00006479891;
            let caliber_inches = diameter / 0.0254;
            let segments = BCSegmentEstimator::estimate_bc_segments(
                bc,
                caliber_inches,
                weight_grains,
                "", // No specific model
                match drag_model {
                    DragModelArg::G1 => "G1",
                    DragModelArg::G7 => "G7",
                },
            );
            Some(segments)
        } else {
            None
        },
        use_enhanced_spin_drift: enable_spin_drift,
        use_form_factor: false,
        enable_wind_shear,
        wind_shear_model: if enable_wind_shear {
            "exponential".to_string()
        } else {
            "none".to_string()
        },
        enable_trajectory_sampling: sample_trajectory,
        sample_interval,
        enable_pitch_damping,
        enable_precession_nutation: enable_precession,
        use_cluster_bc,

        // Optional data
        bc_type_str: None,
        custom_drag_table: None,
    };

    // Set up wind conditions
    let wind = WindConditions {
        speed: wind_speed,
        direction: wind_direction.to_radians(),
        ..Default::default()
    };

    // Set up atmospheric conditions
    let atmosphere = AtmosphericConditions {
        temperature,
        pressure,
        humidity,
        altitude,
        ..Default::default()
    };

    // Create solver
    let mut solver = TrajectorySolver::new(inputs.clone(), wind, atmosphere.clone());
    solver.set_max_range(max_range);
    solver.set_time_step(time_step);

    // Solve trajectory
    let result = solver.solve()?;

    // Calculate stability coefficient if twist rate is provided
    let stability = if twist_rate.is_some() && twist_rate.unwrap() > 0.0 {
        ballistics_engine::stability::compute_stability_coefficient(
            &inputs,
            (altitude, temperature, pressure, 1.0),
        )
    } else {
        0.0
    };

    // Calculate spin drift if enabled and twist rate is provided
    let spin_drift = if enable_spin_drift && twist_rate.is_some() && stability > 0.0 {
        // Calculate spin decay factor based on time of flight
        use ballistics_engine::spin_decay::{
            calculate_spin_decay_correction_factor, SpinDecayParameters,
        };
        let decay_params = SpinDecayParameters::new();
        let avg_velocity = (velocity + result.impact_velocity) / 2.0;

        // Convert units for spin decay calculation
        let mass_grains = mass / 0.00006479891;
        let caliber_inches = diameter / 0.0254;
        let length_inches = diameter * 4.5 / 0.0254; // Approximate length
        let air_density = 1.225; // Standard air density at sea level

        let decay_factor = calculate_spin_decay_correction_factor(
            result.time_of_flight,
            avg_velocity,
            air_density,
            mass_grains,
            caliber_inches,
            length_inches,
            Some(&decay_params),
        );

        // Calculate spin drift with decay
        ballistics_engine::stability::compute_spin_drift_with_decay(
            result.time_of_flight,
            stability,
            twist_rate.unwrap(),
            twist_right,
            Some(decay_factor),
        )
    } else {
        0.0
    };

    // Format output
    match output {
        OutputFormat::Json => {
            let trajectory_result = TrajectoryResult {
                max_range: result.max_range,
                max_height: result.max_height,
                time_of_flight: result.time_of_flight,
                impact_velocity: result.impact_velocity,
                impact_energy: result.impact_energy,
                stability_coefficient: if stability > 0.0 {
                    Some(stability)
                } else {
                    None
                },
                spin_drift: if spin_drift.abs() > 0.0001 {
                    Some(spin_drift)
                } else {
                    None
                },
                trajectory: if full {
                    result
                        .points
                        .into_iter()
                        .map(|p| TrajectoryPoint {
                            time: p.time,
                            x: p.position.x,
                            y: p.position.y,
                            z: p.position.z,
                            velocity: p.velocity_magnitude,
                            energy: p.kinetic_energy,
                        })
                        .collect()
                } else {
                    vec![]
                },
            };
            println!("{}", serde_json::to_string_pretty(&trajectory_result)?);
        }

        OutputFormat::Csv => {
            // Determine unit labels for CSV header
            let (dist_unit, vel_unit, energy_unit) = match units {
                UnitSystem::Metric => ("m", "m/s", "J"),
                UnitSystem::Imperial => ("yd", "fps", "ft-lb"),
            };

            if full {
                // Check if we have sampled points (from --sample-trajectory)
                if let Some(ref sampled) = result.sampled_points {
                    // Output sampled trajectory at regular distance intervals
                    println!(
                        "distance_{},drop_{},drift_{},velocity_{},energy_{},time_s",
                        dist_unit, dist_unit, dist_unit, vel_unit, energy_unit
                    );
                    for s in sampled {
                        let distance = UnitConverter::distance_from_metric(s.distance_m, units);
                        let drop = UnitConverter::distance_from_metric(s.drop_m, units);
                        let drift = UnitConverter::distance_from_metric(s.wind_drift_m, units);
                        let vel = UnitConverter::velocity_from_metric(s.velocity_mps, units);
                        let energy = UnitConverter::energy_from_metric(s.energy_j, units);
                        println!(
                            "{:.2},{:.2},{:.2},{:.2},{:.2},{:.4}",
                            distance, drop, drift, vel, energy, s.time_s
                        );
                    }
                } else {
                    // Output raw trajectory points (all integration steps)
                    println!(
                        "time_s,x_{},y_{},z_{},velocity_{},energy_{}",
                        dist_unit, dist_unit, dist_unit, vel_unit, energy_unit
                    );
                    for p in result.points {
                        let x = UnitConverter::distance_from_metric(p.position.x, units);
                        let y = UnitConverter::distance_from_metric(p.position.y, units);
                        let z = UnitConverter::distance_from_metric(p.position.z, units);
                        let vel = UnitConverter::velocity_from_metric(p.velocity_magnitude, units);
                        let energy = UnitConverter::energy_from_metric(p.kinetic_energy, units);
                        println!(
                            "{:.4},{:.2},{:.2},{:.2},{:.2},{:.2}",
                            p.time, x, y, z, vel, energy
                        );
                    }
                }
            } else {
                println!("metric,value,unit");
                println!(
                    "max_range,{:.2},{}",
                    UnitConverter::distance_from_metric(result.max_range, units),
                    dist_unit
                );
                println!(
                    "max_height,{:.2},{}",
                    UnitConverter::distance_from_metric(result.max_height, units),
                    dist_unit
                );
                println!("time_of_flight,{:.4},s", result.time_of_flight);
                println!(
                    "impact_velocity,{:.2},{}",
                    UnitConverter::velocity_from_metric(result.impact_velocity, units),
                    vel_unit
                );
                println!(
                    "impact_energy,{:.2},{}",
                    UnitConverter::energy_from_metric(result.impact_energy, units),
                    energy_unit
                );
                if stability > 0.0 {
                    println!("stability_coefficient,{:.2},", stability);
                }
                if spin_drift.abs() > 0.0001 {
                    println!(
                        "spin_drift,{:.4},{}",
                        UnitConverter::distance_from_metric(spin_drift, units),
                        dist_unit
                    );
                }
            }
        }

        OutputFormat::Table => {
            // Convert outputs to user's units
            let range_display = UnitConverter::distance_from_metric(result.max_range, units);
            let height_display = UnitConverter::distance_from_metric(result.max_height, units);
            let velocity_display =
                UnitConverter::velocity_from_metric(result.impact_velocity, units);
            let energy_display = UnitConverter::energy_from_metric(result.impact_energy, units);

            let (range_unit, velocity_unit, energy_unit) = match units {
                UnitSystem::Metric => ("m", "m/s", "J"),
                UnitSystem::Imperial => ("yd", "fps", "ft-lb"),
            };

            println!("╔════════════════════════════════════════╗");
            println!("║         TRAJECTORY RESULTS             ║");
            println!("╠════════════════════════════════════════╣");
            println!(
                "║ Max Range:         {:>8.2} {:3}       ║",
                range_display, range_unit
            );
            println!(
                "║ Max Height:        {:>8.2} {:3}       ║",
                height_display, range_unit
            );
            println!(
                "║ Time of Flight:    {:>8.3} s          ║",
                result.time_of_flight
            );
            println!(
                "║ Impact Velocity:   {:>8.2} {:3}       ║",
                velocity_display, velocity_unit
            );
            println!(
                "║ Impact Energy:     {:>8.2} {:5}     ║",
                energy_display, energy_unit
            );
            if stability > 0.0 {
                println!("╠════════════════════════════════════════╣");
                println!("║ Stability (SG):    {:>8.2}            ║", stability);
                let stability_status = if stability < 1.0 {
                    "UNSTABLE"
                } else if stability < 1.5 {
                    "MARGINAL"
                } else {
                    "STABLE  "
                };
                println!("║ Status:            {:>8}            ║", stability_status);
            }
            if spin_drift.abs() > 0.0001 {
                let drift_display = UnitConverter::distance_from_metric(spin_drift.abs(), units);
                let direction = if spin_drift > 0.0 { "Right" } else { "Left" };
                println!("╠════════════════════════════════════════╣");
                println!(
                    "║ Spin Drift:        {:>8.2} {:3}       ║",
                    drift_display, range_unit
                );
                println!("║ Direction:         {:>8}            ║", direction);
            }

            // Display pitch damping analysis if available
            if let Some(min_damping) = result.min_pitch_damping {
                println!("╠════════════════════════════════════════╣");
                println!("║ Min Pitch Damping: {:>8.3}            ║", min_damping);
                let stability_warning = if min_damping > 0.0 {
                    "UNSTABLE" // Positive damping in transonic can cause instability
                } else if min_damping > -0.2 {
                    "MARGINAL"
                } else {
                    "STABLE  "
                };
                println!("║ Transonic Status:  {:>8}            ║", stability_warning);
                if let Some(mach) = result.transonic_mach {
                    println!("║ Entered at Mach:   {:>8.2}            ║", mach);
                }
            }

            // Display angular motion analysis if available
            if let Some(ref angular_state) = result.angular_state {
                println!("╠════════════════════════════════════════╣");
                println!("║       Angular Motion Analysis          ║");
                println!("╠════════════════════════════════════════╣");
                println!(
                    "║ Final Pitch Angle: {:>8.4} rad        ║",
                    angular_state.pitch_angle
                );
                println!(
                    "║ Final Yaw Angle:   {:>8.4} rad        ║",
                    angular_state.yaw_angle
                );
                if let Some(max_yaw) = result.max_yaw_angle {
                    println!("║ Max Yaw Angle:     {:>8.4} rad        ║", max_yaw);
                    println!(
                        "║                   ({:>8.2}°)          ║",
                        max_yaw.to_degrees()
                    );
                }
                if let Some(max_prec) = result.max_precession_angle {
                    println!("║ Max Precession:    {:>8.4} rad        ║", max_prec);
                }
                println!(
                    "║ Nutation Phase:    {:>8.2} rad        ║",
                    angular_state.nutation_phase
                );
            }

            println!("╚════════════════════════════════════════╝");

            // Check if trajectory hit ground
            if let Some(last_point) = result.points.last() {
                let last_height = last_point.position.y;
                let last_range = last_point.position.z;

                // Ground threshold is typically around 0, so if y is close to or below 0, it hit ground
                if last_height <= 0.001 {
                    println!();
                    let range_display = UnitConverter::distance_from_metric(last_range, units);
                    println!(
                        "Bullet struck ground at {:.0} {}",
                        range_display, range_unit
                    );
                }
            }

            if full && !result.points.is_empty() {
                println!();
                println!("Trajectory Points:");
                let (dist_hdr, vel_hdr, energy_hdr) = match units {
                    UnitSystem::Metric => ("(m)", "(m/s)", "(J)"),
                    UnitSystem::Imperial => ("(yd)", "(fps)", "(ft-lb)"),
                };
                println!("┌──────────┬──────────┬──────────┬──────────┬──────────┐");
                println!(
                    "│ Time (s) │  X {:5} │  Y {:5} │ Vel{:5} │Energy{:5}│",
                    dist_hdr, dist_hdr, vel_hdr, energy_hdr
                );
                println!("├──────────┼──────────┼──────────┼──────────┼──────────┤");

                let step = if result.points.len() > 20 {
                    result.points.len() / 20
                } else {
                    1
                };

                for (i, p) in result.points.iter().enumerate() {
                    if i % step == 0 || i == result.points.len() - 1 {
                        let x_display = UnitConverter::distance_from_metric(p.position.x, units);
                        let y_display = UnitConverter::distance_from_metric(p.position.y, units);
                        let vel_display =
                            UnitConverter::velocity_from_metric(p.velocity_magnitude, units);
                        let energy_display =
                            UnitConverter::energy_from_metric(p.kinetic_energy, units);

                        println!(
                            "│ {:>8.3} │ {:>8.2} │ {:>8.2} │ {:>8.2} │ {:>8.2} │",
                            p.time, x_display, y_display, vel_display, energy_display
                        );
                    }
                }
                println!("└──────────┴──────────┴──────────┴──────────┴──────────┘");
            }

            // Display sampled trajectory points if available
            if let Some(ref samples) = result.sampled_points {
                if !samples.is_empty() {
                    println!();
                    println!(
                        "Sampled Trajectory (every {:.0} {})",
                        UnitConverter::distance_from_metric(sample_interval, units),
                        match units {
                            UnitSystem::Metric => "m",
                            UnitSystem::Imperial => "yd",
                        }
                    );

                    let (dist_hdr, drop_hdr, drift_hdr, vel_hdr) = match units {
                        UnitSystem::Metric => ("(m)", "(m)", "(m)", "(m/s)"),
                        UnitSystem::Imperial => ("(yd)", "(in)", "(in)", "(fps)"),
                    };

                    println!("┌──────────┬──────────┬──────────┬──────────┬──────────┐");
                    println!(
                        "│ Dist{:4} │ Drop{:4} │Drift{:4} │ Vel{:5} │  Flags   │",
                        dist_hdr, drop_hdr, drift_hdr, vel_hdr
                    );
                    println!("├──────────┼──────────┼──────────┼──────────┼──────────┤");

                    for sample in samples.iter() {
                        let dist_display =
                            UnitConverter::distance_from_metric(sample.distance_m, units);
                        let drop_display = if units == UnitSystem::Imperial {
                            sample.drop_m * 39.3701 // Convert to inches for imperial
                        } else {
                            sample.drop_m
                        };
                        let drift_display = if units == UnitSystem::Imperial {
                            sample.wind_drift_m * 39.3701 // Convert to inches for imperial
                        } else {
                            sample.wind_drift_m
                        };
                        let vel_display =
                            UnitConverter::velocity_from_metric(sample.velocity_mps, units);

                        let flags_str = sample
                            .flags
                            .iter()
                            .map(|f| match f {
                                trajectory_sampling::TrajectoryFlag::ZeroCrossing => "Zero",
                                trajectory_sampling::TrajectoryFlag::MachTransition => "Mach",
                                trajectory_sampling::TrajectoryFlag::Apex => "Apex",
                            })
                            .collect::<Vec<_>>()
                            .join(",");

                        println!(
                            "│ {:>8.1} │ {:>8.2} │ {:>8.2} │ {:>8.1} │ {:8} │",
                            dist_display,
                            drop_display,
                            drift_display,
                            vel_display,
                            if flags_str.is_empty() {
                                "-"
                            } else {
                                &flags_str
                            }
                        );
                    }

                    println!("└──────────┴──────────┴──────────┴──────────┴──────────┘");
                }
            }
        }
    }

    Ok(())
}

/// Display trajectory results from API response
#[cfg(feature = "online")]
fn display_api_trajectory_result(
    response: &ballistics_engine::api_client::TrajectoryResponse,
    output: OutputFormat,
    units: UnitSystem,
    full: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    match output {
        OutputFormat::Json => {
            // Re-serialize the response as JSON
            let result = serde_json::json!({
                "source": "api",
                "zero_angle_degrees": response.zero_angle.to_degrees(),
                "time_of_flight": response.time_of_flight,
                "bc_confidence": response.bc_confidence,
                "ml_corrections_applied": response.ml_corrections_applied,
                "max_ordinate": response.max_ordinate,
                "impact_velocity": response.impact_velocity,
                "impact_energy": response.impact_energy,
                "trajectory": if full {
                    response.trajectory.iter().map(|p| {
                        serde_json::json!({
                            "range": p.range,
                            "drop": p.drop,
                            "drift": p.drift,
                            "velocity": p.velocity,
                            "energy": p.energy,
                            "time": p.time
                        })
                    }).collect::<Vec<_>>()
                } else {
                    vec![]
                }
            });
            println!("{}", serde_json::to_string_pretty(&result)?);
        }

        OutputFormat::Csv => {
            let (dist_unit, vel_unit, energy_unit) = match units {
                UnitSystem::Metric => ("m", "m/s", "J"),
                UnitSystem::Imperial => ("yd", "fps", "ft-lb"),
            };

            if full {
                println!(
                    "range_{},drop_{},drift_{},velocity_{},energy_{},time_s",
                    dist_unit, dist_unit, dist_unit, vel_unit, energy_unit
                );
                for p in &response.trajectory {
                    let range = UnitConverter::distance_from_metric(p.range, units);
                    let drop = if units == UnitSystem::Imperial {
                        p.drop * 39.3701 // meters to inches
                    } else {
                        p.drop
                    };
                    let drift = if units == UnitSystem::Imperial {
                        p.drift * 39.3701
                    } else {
                        p.drift
                    };
                    let vel = UnitConverter::velocity_from_metric(p.velocity, units);
                    let energy = UnitConverter::energy_from_metric(p.energy, units);
                    println!(
                        "{:.2},{:.2},{:.2},{:.2},{:.2},{:.4}",
                        range, drop, drift, vel, energy, p.time
                    );
                }
            } else {
                println!("metric,value,unit");
                println!("source,api,");
                println!("time_of_flight,{:.4},s", response.time_of_flight);
                println!(
                    "zero_angle,{:.4},degrees",
                    response.zero_angle.to_degrees()
                );
                if let Some(confidence) = response.bc_confidence {
                    println!("bc_confidence,{:.2},%", confidence * 100.0);
                }
                if let Some(ref corrections) = response.ml_corrections_applied {
                    println!("ml_corrections,{},", corrections.join(";"));
                }
            }
        }

        OutputFormat::Table => {
            let (range_unit, vel_unit, energy_unit) = match units {
                UnitSystem::Metric => ("m", "m/s", "J"),
                UnitSystem::Imperial => ("yd", "fps", "ft-lb"),
            };

            println!("╔════════════════════════════════════════╗");
            println!("║    TRAJECTORY RESULTS (API/ML)         ║");
            println!("╠════════════════════════════════════════╣");

            // Display ML corrections if available
            if let Some(ref corrections) = response.ml_corrections_applied {
                if !corrections.is_empty() {
                    println!("║ ML Corrections Applied:                ║");
                    for correction in corrections.iter().take(3) {
                        let truncated: String = correction.chars().take(32).collect();
                        println!("║   • {:32} ║", truncated);
                    }
                    println!("╠════════════════════════════════════════╣");
                }
            }

            // Display BC confidence if available
            if let Some(confidence) = response.bc_confidence {
                println!("║ BC Confidence:     {:>8.1}%           ║", confidence * 100.0);
                println!("╠════════════════════════════════════════╣");
            }

            println!(
                "║ Time of Flight:    {:>8.3} s          ║",
                response.time_of_flight
            );
            println!(
                "║ Zero Angle:        {:>8.4}°          ║",
                response.zero_angle.to_degrees()
            );

            if let Some(max_ord) = response.max_ordinate {
                let max_ord_display = UnitConverter::distance_from_metric(max_ord, units);
                println!(
                    "║ Max Ordinate:      {:>8.2} {:3}       ║",
                    max_ord_display, range_unit
                );
            }

            if let Some(impact_vel) = response.impact_velocity {
                let vel_display = UnitConverter::velocity_from_metric(impact_vel, units);
                println!(
                    "║ Impact Velocity:   {:>8.2} {:3}       ║",
                    vel_display, vel_unit
                );
            }

            if let Some(impact_e) = response.impact_energy {
                let energy_display = UnitConverter::energy_from_metric(impact_e, units);
                println!(
                    "║ Impact Energy:     {:>8.2} {:5}     ║",
                    energy_display, energy_unit
                );
            }

            println!("╚════════════════════════════════════════╝");

            // Display trajectory points if full mode
            if full && !response.trajectory.is_empty() {
                println!();
                println!("Trajectory Points (from API):");

                let (dist_hdr, drop_hdr, drift_hdr, vel_hdr) = match units {
                    UnitSystem::Metric => ("(m)", "(m)", "(m)", "(m/s)"),
                    UnitSystem::Imperial => ("(yd)", "(in)", "(in)", "(fps)"),
                };

                println!("┌──────────┬──────────┬──────────┬──────────┬──────────┐");
                println!(
                    "│Range{:4} │ Drop{:4} │Drift{:4} │ Vel{:5} │ Time (s) │",
                    dist_hdr, drop_hdr, drift_hdr, vel_hdr
                );
                println!("├──────────┼──────────┼──────────┼──────────┼──────────┤");

                let step = if response.trajectory.len() > 20 {
                    response.trajectory.len() / 20
                } else {
                    1
                };

                for (i, p) in response.trajectory.iter().enumerate() {
                    if i % step == 0 || i == response.trajectory.len() - 1 {
                        let range_display = UnitConverter::distance_from_metric(p.range, units);
                        let drop_display = if units == UnitSystem::Imperial {
                            p.drop * 39.3701
                        } else {
                            p.drop
                        };
                        let drift_display = if units == UnitSystem::Imperial {
                            p.drift * 39.3701
                        } else {
                            p.drift
                        };
                        let vel_display = UnitConverter::velocity_from_metric(p.velocity, units);

                        println!(
                            "│ {:>8.1} │ {:>8.2} │ {:>8.2} │ {:>8.1} │ {:>8.4} │",
                            range_display, drop_display, drift_display, vel_display, p.time
                        );
                    }
                }
                println!("└──────────┴──────────┴──────────┴──────────┴──────────┘");
            }
        }
    }

    Ok(())
}

fn run_monte_carlo(
    velocity: f64,
    angle: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    num_sims: usize,
    velocity_std: f64,
    angle_std: f64,
    bc_std: f64,
    wind_std: f64,
    wind_speed: f64,
    wind_direction: f64,
    target_distance: Option<f64>,
    output: MonteCarloOutput,
) -> Result<(), Box<dyn Error>> {
    // Create base inputs
    let base_inputs = BallisticInputs {
        muzzle_velocity: velocity,
        muzzle_angle: angle.to_radians(),
        bc_value: bc,
        bullet_mass: mass,
        bullet_diameter: diameter,
        ..Default::default()
    };

    // Create base wind conditions
    let base_wind = WindConditions {
        speed: wind_speed,
        direction: wind_direction.to_radians(),
    };

    // Set up Monte Carlo parameters
    let mc_params = MonteCarloParams {
        num_simulations: num_sims,
        velocity_std_dev: velocity_std,
        angle_std_dev: angle_std.to_radians(),
        bc_std_dev: bc_std,
        wind_speed_std_dev: wind_std,
        target_distance,
        base_wind_speed: wind_speed,
        base_wind_direction: wind_direction.to_radians(),
        azimuth_std_dev: angle_std.to_radians() * 0.5, // Use half of elevation std for horizontal spread
    };

    // Run Monte Carlo simulation
    let results = ballistics_engine::run_monte_carlo_with_wind(base_inputs, base_wind, mc_params)?;

    // Calculate statistics
    let mean_range = results.ranges.iter().sum::<f64>() / results.ranges.len() as f64;
    let variance_range: f64 = results
        .ranges
        .iter()
        .map(|r| (r - mean_range).powi(2))
        .sum::<f64>()
        / results.ranges.len() as f64;
    let std_range = variance_range.sqrt();

    let mean_velocity =
        results.impact_velocities.iter().sum::<f64>() / results.impact_velocities.len() as f64;
    let variance_velocity: f64 = results
        .impact_velocities
        .iter()
        .map(|v| (v - mean_velocity).powi(2))
        .sum::<f64>()
        / results.impact_velocities.len() as f64;
    let std_velocity = variance_velocity.sqrt();

    // Calculate CEP (simplified - actual CEP calculation would need lateral dispersion)
    let cep = std_range * 1.1774; // Approximation for circular normal distribution

    // Calculate hit probability if target distance specified
    let hit_probability = target_distance.map(|target| {
        let hits = results
            .ranges
            .iter()
            .filter(|r| (*r - target).abs() < 1.0) // Within 1m of target
            .count();
        hits as f64 / results.ranges.len() as f64
    });

    match output {
        MonteCarloOutput::Summary => {
            println!("╔════════════════════════════════════════╗");
            println!("║      MONTE CARLO RESULTS               ║");
            println!("║      {} simulations                   ║", num_sims);
            println!("╠════════════════════════════════════════╣");
            println!("║ Mean Range:        {:>8.2} m          ║", mean_range);
            println!("║ Std Dev Range:     {:>8.2} m          ║", std_range);
            println!("║ Mean Impact Vel:   {:>8.2} m/s        ║", mean_velocity);
            println!("║ Std Dev Velocity:  {:>8.2} m/s        ║", std_velocity);
            println!("║ CEP (approx):      {:>8.2} m          ║", cep);
            if let Some(prob) = hit_probability {
                println!("║ Hit Probability:   {:>8.1} %          ║", prob * 100.0);
            }
            println!("╚════════════════════════════════════════╝");
        }

        MonteCarloOutput::Full => {
            let mc_result = MonteCarloResult {
                mean_range,
                std_range,
                mean_impact_velocity: mean_velocity,
                std_impact_velocity: std_velocity,
                cep,
                hit_probability,
            };
            println!("{}", serde_json::to_string_pretty(&mc_result)?);
        }

        MonteCarloOutput::Statistics => {
            println!("range_min,range_max,range_mean,range_std,vel_min,vel_max,vel_mean,vel_std");
            let range_min = results.ranges.iter().fold(f64::INFINITY, |a, &b| a.min(b));
            let range_max = results
                .ranges
                .iter()
                .fold(f64::NEG_INFINITY, |a, &b| a.max(b));
            let vel_min = results
                .impact_velocities
                .iter()
                .fold(f64::INFINITY, |a, &b| a.min(b));
            let vel_max = results
                .impact_velocities
                .iter()
                .fold(f64::NEG_INFINITY, |a, &b| a.max(b));

            println!(
                "{:.2},{:.2},{:.2},{:.2},{:.2},{:.2},{:.2},{:.2}",
                range_min,
                range_max,
                mean_range,
                std_range,
                vel_min,
                vel_max,
                mean_velocity,
                std_velocity
            );
        }
    }

    Ok(())
}

fn run_zero_calculation(
    velocity: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    target_distance: f64,
    target_height: f64,
    sight_height: f64,
    output: OutputFormat,
    units: UnitSystem,
) -> Result<(), Box<dyn Error>> {
    // Create ballistic inputs
    let inputs = BallisticInputs {
        muzzle_velocity: velocity,
        bc_value: bc,
        bullet_mass: mass,
        bullet_diameter: diameter,
        sight_height,
        ..Default::default()
    };

    // Calculate zero angle
    let zero_angle =
        ballistics_engine::calculate_zero_angle(inputs.clone(), target_distance, target_height)?;

    // Calculate trajectory at zero angle to get additional info
    let mut zeroed_inputs = inputs;
    zeroed_inputs.muzzle_angle = zero_angle;

    let solver = TrajectorySolver::new(zeroed_inputs, Default::default(), Default::default());
    let trajectory = solver.solve()?;

    match output {
        OutputFormat::Json => {
            let result = serde_json::json!({
                "zero_angle_degrees": zero_angle.to_degrees(),
                "zero_angle_moa": zero_angle.to_degrees() * 60.0,
                "zero_angle_mrad": zero_angle * 1000.0,
                "sight_adjustment_moa": (zero_angle.to_degrees() * 60.0) - (sight_height / target_distance * 3437.75),
                "max_ordinate": trajectory.max_height,
                "point_blank_range": trajectory.points.iter()
                    .find(|p| p.position.y < -0.05)
                    .map(|p| p.position.x)
                    .unwrap_or(trajectory.max_range),
            });
            println!("{}", serde_json::to_string_pretty(&result)?);
        }

        OutputFormat::Csv => {
            println!("metric,value,unit");
            println!("zero_angle,{:.4},degrees", zero_angle.to_degrees());
            println!("zero_angle_moa,{:.2},MOA", zero_angle.to_degrees() * 60.0);
            println!("zero_angle_mrad,{:.2},mrad", zero_angle * 1000.0);
            println!("max_ordinate,{:.3},meters", trajectory.max_height);
        }

        OutputFormat::Table => {
            // Convert distances back to display units
            let target_dist_display = UnitConverter::distance_from_metric(target_distance, units);
            let target_height_display = UnitConverter::distance_from_metric(target_height, units);
            let sight_height_display = UnitConverter::distance_from_metric(sight_height, units);
            let max_ordinate_display =
                UnitConverter::distance_from_metric(trajectory.max_height, units);

            let dist_unit = match units {
                UnitSystem::Metric => "m",
                UnitSystem::Imperial => "yd",
            };

            println!("╔════════════════════════════════════════╗");
            println!("║          ZERO CALCULATION              ║");
            println!("╠════════════════════════════════════════╣");
            println!(
                "║ Target Distance:   {:>8.1} {:3}       ║",
                target_dist_display, dist_unit
            );
            println!(
                "║ Target Height:     {:>8.2} {:3}       ║",
                target_height_display, dist_unit
            );
            println!(
                "║ Sight Height:      {:>8.3} {:3}       ║",
                sight_height_display, dist_unit
            );
            println!("╠════════════════════════════════════════╣");
            println!(
                "║ Zero Angle:        {:>8.4}°          ║",
                zero_angle.to_degrees()
            );
            println!(
                "║ Zero Angle (MOA):  {:>8.2} MOA        ║",
                zero_angle.to_degrees() * 60.0
            );
            println!(
                "║ Zero Angle (mrad): {:>8.2} mrad       ║",
                zero_angle * 1000.0
            );
            println!(
                "║ Max Ordinate:      {:>8.3} {:3}       ║",
                max_ordinate_display, dist_unit
            );
            println!("╚════════════════════════════════════════╝");
        }
    }

    Ok(())
}

fn run_bc_segment_generation(
    bc: f64,
    mass: f64,
    diameter: f64,
    model: &str,
    drag_model: &str,
    output: OutputFormat,
    units: UnitSystem,
) -> Result<(), Box<dyn Error>> {
    use ballistics_engine::bc_estimation::BCSegmentEstimator;

    // Convert mass to grains and diameter to inches for the estimation
    let weight_grains = mass / 0.00006479891;
    let caliber_inches = diameter / 0.0254;

    // Generate BC segments
    let segments = BCSegmentEstimator::estimate_bc_segments(
        bc,
        caliber_inches,
        weight_grains,
        model,
        drag_model,
    );

    match output {
        OutputFormat::Json => {
            let segments_json: Vec<_> = segments
                .iter()
                .map(|s| {
                    serde_json::json!({
                        "velocity_min_fps": s.velocity_min,
                        "velocity_max_fps": s.velocity_max,
                        "bc": s.bc_value
                    })
                })
                .collect();

            let result = serde_json::json!({
                "base_bc": bc,
                "mass_grains": weight_grains,
                "caliber_inches": caliber_inches,
                "model": model,
                "drag_model": drag_model,
                "segments": segments_json
            });

            println!("{}", serde_json::to_string_pretty(&result)?);
        }

        OutputFormat::Csv => {
            println!("velocity_min_fps,velocity_max_fps,bc");
            for segment in &segments {
                println!(
                    "{},{},{:.4}",
                    segment.velocity_min, segment.velocity_max, segment.bc_value
                );
            }
        }

        OutputFormat::Table => {
            println!("╔════════════════════════════════════════╗");
            println!("║        BC SEGMENT GENERATION           ║");
            println!("╠════════════════════════════════════════╣");
            println!("║ Base BC:         {:.3}                 ║", bc);
            println!("║ Caliber:         {:.3} inches          ║", caliber_inches);
            println!("║ Weight:          {:.1} grains          ║", weight_grains);
            if !model.is_empty() {
                println!("║ Model:           {:20}  ║", model);
            }
            println!("║ Drag Model:      {:20}  ║", drag_model);
            println!("╠════════════════════════════════════════╣");
            println!("║         VELOCITY SEGMENTS              ║");
            println!("╠════════════════════════════════════════╣");

            for segment in &segments {
                let vel_display = if units == UnitSystem::Imperial {
                    format!(
                        "{:.0}-{:.0} fps",
                        segment.velocity_min, segment.velocity_max
                    )
                } else {
                    format!(
                        "{:.0}-{:.0} m/s",
                        segment.velocity_min / 3.28084,
                        segment.velocity_max / 3.28084
                    )
                };
                println!("║ {:18} → BC: {:.4}     ║", vel_display, segment.bc_value);
            }

            println!("╚════════════════════════════════════════╝");
        }
    }

    Ok(())
}

fn run_bc_estimation(
    velocity: f64,
    mass: f64,
    diameter: f64,
    distance1: f64,
    drop1: f64,
    distance2: f64,
    drop2: f64,
    output: OutputFormat,
) -> Result<(), Box<dyn Error>> {
    // Create trajectory points for BC estimation
    let points = vec![(distance1, drop1), (distance2, drop2)];

    // Estimate BC
    let estimated_bc =
        ballistics_engine::estimate_bc_from_trajectory(velocity, mass, diameter, &points)?;

    // Verify the estimation by running a trajectory
    let inputs = BallisticInputs {
        muzzle_velocity: velocity,
        bc_value: estimated_bc,
        bullet_mass: mass,
        bullet_diameter: diameter,
        ..Default::default()
    };

    let solver = TrajectorySolver::new(inputs, Default::default(), Default::default());
    let trajectory = solver.solve()?;

    // Find drops at the specified distances (Z is downrange)
    let calc_drop1 = trajectory
        .points
        .iter()
        .find(|p| p.position.z >= distance1)
        .map(|p| -p.position.y)
        .unwrap_or(0.0);

    let calc_drop2 = trajectory
        .points
        .iter()
        .find(|p| p.position.z >= distance2)
        .map(|p| -p.position.y)
        .unwrap_or(0.0);

    let error1 = ((calc_drop1 - drop1) / drop1 * 100.0).abs();
    let error2 = ((calc_drop2 - drop2) / drop2 * 100.0).abs();

    match output {
        OutputFormat::Json => {
            let result = serde_json::json!({
                "estimated_bc": estimated_bc,
                "verification": {
                    "distance1_m": distance1,
                    "actual_drop1_m": drop1,
                    "calculated_drop1_m": calc_drop1,
                    "error1_percent": error1,
                    "distance2_m": distance2,
                    "actual_drop2_m": drop2,
                    "calculated_drop2_m": calc_drop2,
                    "error2_percent": error2,
                }
            });
            println!("{}", serde_json::to_string_pretty(&result)?);
        }

        OutputFormat::Csv => {
            println!("metric,value");
            println!("estimated_bc,{:.4}", estimated_bc);
            println!("error_at_{}m_percent,{:.2}", distance1, error1);
            println!("error_at_{}m_percent,{:.2}", distance2, error2);
        }

        OutputFormat::Table => {
            println!("╔════════════════════════════════════════╗");
            println!("║         BC ESTIMATION RESULT           ║");
            println!("╠════════════════════════════════════════╣");
            println!("║ Estimated BC:      {:>8.4}            ║", estimated_bc);
            println!("╠════════════════════════════════════════╣");
            println!("║ Verification:                          ║");
            println!("║ At {:.0}m:                             ║", distance1);
            println!("║   Actual drop:     {:>8.3} m          ║", drop1);
            println!("║   Calculated:      {:>8.3} m          ║", calc_drop1);
            println!("║   Error:           {:>8.2} %          ║", error1);
            println!("║ At {:.0}m:                             ║", distance2);
            println!("║   Actual drop:     {:>8.3} m          ║", drop2);
            println!("║   Calculated:      {:>8.3} m          ║", calc_drop2);
            println!("║   Error:           {:>8.2} %          ║", error2);
            println!("╚════════════════════════════════════════╝");
        }
    }

    Ok(())
}
