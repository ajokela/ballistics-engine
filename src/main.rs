// Allocator selection based on features. jemalloc and mimalloc each define a #[global_allocator],
// and a crate may have at most one — so when both features are enabled (e.g. `--all-features`,
// which docs.rs uses by default) jemalloc yields to mimalloc to keep the build compiling.
#[cfg(all(
    feature = "jemalloc",
    not(feature = "mimalloc"),
    not(target_env = "msvc")
))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[cfg(all(feature = "mimalloc", not(target_env = "msvc")))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[cfg(feature = "pdf")]
mod pdf_dope_card;
#[cfg(feature = "pdf")]
use pdf_dope_card::{DopeCardConfig, DopeCardRow, FontSizePreset, calculate_density_altitude, yards_to_mil, calculate_lead_mil};

use ballistics_engine::{
    trajectory_sampling, AtmosphericConditions, BallisticInputs, BCSegmentData, DragModel,
    MonteCarloParams, TrajectorySolver, WindConditions,
};
#[cfg(feature = "online")]
use ballistics_engine::api_client::{ApiClient, TrajectoryRequestBuilder, TrueVelocityRequest};
#[cfg(feature = "online")]
use ballistics_engine::bc_table_download::Bc5dDownloader;
use ballistics_engine::bc_table::BcCorrectionTable;
use ballistics_engine::bc_table_5d::Bc5dTableManager;
use clap::{CommandFactory, Parser, Subcommand, ValueEnum};
use serde::{Deserialize, Serialize};
use strsim::levenshtein;
use std::collections::HashMap;
use std::error::Error;
use std::fs;
use std::io::{self, Write};
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

/// Stable FNV-1a hash (not affected by Rust version changes unlike DefaultHasher)
#[cfg(feature = "online")]
fn fnv1a_hash(data: &[u8]) -> u64 {
    let mut hash: u64 = 0xcbf29ce484222325;
    for &byte in data {
        hash ^= byte as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

/// Calculate a stable hash of the TOS text for tracking changes
#[cfg(feature = "online")]
fn hash_tos(text: &str) -> String {
    format!("{:x}", fnv1a_hash(text.as_bytes()))
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

// ============================================================================
// Custom Value Parsers for f64 Range Validation (MBA-732)
// ============================================================================

/// Create a value parser that validates an f64 is within [min, max].
fn f64_range(min: f64, max: f64) -> clap::builder::ValueParser {
    clap::builder::ValueParser::new(F64RangeParser { min, max })
}

#[derive(Clone, Debug)]
struct F64RangeParser {
    min: f64,
    max: f64,
}

impl clap::builder::TypedValueParser for F64RangeParser {
    type Value = f64;

    fn parse_ref(
        &self,
        _cmd: &clap::Command,
        arg: Option<&clap::Arg>,
        value: &std::ffi::OsStr,
    ) -> Result<Self::Value, clap::Error> {
        let s = value.to_str().ok_or_else(|| {
            let mut err = clap::Error::new(clap::error::ErrorKind::InvalidUtf8);
            if let Some(a) = arg {
                err.insert(
                    clap::error::ContextKind::InvalidArg,
                    clap::error::ContextValue::String(a.get_id().to_string()),
                );
            }
            err
        })?;
        let inner: f64 = s.parse().map_err(|_| {
            let mut err = clap::Error::new(clap::error::ErrorKind::ValueValidation);
            if let Some(a) = arg {
                err.insert(
                    clap::error::ContextKind::InvalidArg,
                    clap::error::ContextValue::String(a.get_id().to_string()),
                );
            }
            err.insert(
                clap::error::ContextKind::InvalidValue,
                clap::error::ContextValue::String(format!("'{}' is not a valid number", s)),
            );
            err
        })?;
        if !inner.is_finite() || inner < self.min || inner > self.max {
            let mut err = clap::Error::new(clap::error::ErrorKind::ValueValidation);
            if let Some(a) = arg {
                err.insert(
                    clap::error::ContextKind::InvalidArg,
                    clap::error::ContextValue::String(a.get_id().to_string()),
                );
            }
            err.insert(
                clap::error::ContextKind::InvalidValue,
                clap::error::ContextValue::String(format!(
                    "{} is not in range {:.6}..={:.6}",
                    inner, self.min, self.max
                )),
            );
            return Err(err);
        }
        Ok(inner)
    }
}

#[derive(Parser)]
#[command(name = "ballistics")]
#[command(author = "Ballistics Engine Team")]
#[command(version = env!("CARGO_PKG_VERSION"))]
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

        /// Load a saved profile by name (from ~/.ballistics/profiles/)
        #[arg(long, value_name = "NAME")]
        saved_profile: Option<String>,

        /// Load location/environmental data from CSV file
        #[arg(long, value_name = "FILE")]
        location: Option<PathBuf>,

        /// Site name to load from location CSV (matches first column, e.g., "KF_LR")
        #[arg(long, value_name = "NAME")]
        site: Option<String>,

        /// Initial velocity (fps or m/s based on --units)
        #[arg(short = 'v', long, value_parser = f64_range(0.0, 6000.0))]
        velocity: Option<f64>,

        /// Velocity adjustment (added to base velocity for truing from chronograph data)
        #[arg(long, default_value = "0.0", allow_hyphen_values = true)]
        velocity_adjustment: f64,

        /// Launch angle (degrees)
        #[arg(short = 'a', long, default_value = "0.0")]
        angle: f64,

        /// Ballistic coefficient
        #[arg(short = 'b', long, value_parser = f64_range(0.001, 2.0))]
        bc: Option<f64>,

        /// BC adjustment factor (multiplier for truing, e.g., 0.85 = 85% of stated BC)
        #[arg(long, default_value = "1.0")]
        bc_adjustment: f64,

        /// Mass (grains for imperial, grams for metric)
        #[arg(short = 'm', long, value_parser = f64_range(0.1, 2000.0))]
        mass: Option<f64>,

        /// Diameter (inches for imperial, mm for metric)
        #[arg(short = 'd', long, value_parser = f64_range(0.01, 60.0))]
        diameter: Option<f64>,

        /// Drag model (G1, G7, Custom)
        #[arg(long, default_value = "g1")]
        drag_model: DragModelArg,

        /// Maximum range (yards or meters based on --units)
        #[arg(long, default_value = "1000.0")]
        max_range: f64,

        /// Integration time step in seconds — RK4/Euler only (the adaptive RK45 default ignores it and steps adaptively)
        #[arg(long, default_value = "0.001", value_parser = f64_range(0.00001, 0.1))]
        time_step: f64,

        /// Wind speed (mph or m/s based on --units)
        #[arg(long, default_value = "0.0")]
        wind_speed: f64,

        /// Wind direction (degrees; wind-FROM: 0=headwind, 90=from right, 180=tailwind, 270=from left)
        #[arg(long, default_value = "0.0")]
        wind_direction: f64,

        /// Downrange wind segment "SPEED:ANGLE:UNTIL_DISTANCE" (repeatable). SPEED/UNTIL_DISTANCE
        /// follow --units (mph & yd imperial, m/s & m metric); ANGLE is degrees, same convention
        /// as --wind-direction. Each segment applies from the previous boundary out to
        /// UNTIL_DISTANCE; wind is zero beyond the last segment. Overrides --wind-speed/-direction.
        /// Not compatible with --enable-wind-shear.
        #[arg(long = "wind-segment", value_name = "SPEED:ANGLE:DIST", action = clap::ArgAction::Append)]
        wind_segment: Vec<String>,

        /// Temperature (Fahrenheit or Celsius based on --units)
        #[arg(long, default_value = "59.0", value_parser = f64_range(-100.0, 200.0))]
        temperature: f64,

        /// Pressure (inHg or hPa based on --units)
        #[arg(long, default_value = "29.92", value_parser = f64_range(15.0, 1200.0))]
        pressure: f64,

        /// Humidity (0-100%)
        #[arg(long, default_value = "50.0", value_parser = f64_range(0.0, 100.0))]
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

        /// Apply aerodynamic jump as a muzzle launch-angle perturbation (EXPERIMENTAL — heuristic model, MBA-959)
        #[arg(long)]
        enable_aerodynamic_jump: bool,

        /// Enable wind shear (altitude-dependent wind)
        #[arg(long)]
        enable_wind_shear: bool,

        /// Enable trajectory sampling at regular intervals
        #[arg(long)]
        sample_trajectory: bool,

        /// Sampling interval in meters (default: 10)
        /// Trajectory sampling interval in meters (used with --sample-trajectory; always metric, not unit-system dependent)
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

        /// Path to BC correction table file for offline ML-enhanced corrections
        #[arg(long, value_name = "FILE", help = "Use precomputed BC correction table instead of online API")]
        bc_table: Option<PathBuf>,

        /// Directory containing caliber-specific BC5D tables (e.g., bc5d_308.bin)
        #[arg(long, value_name = "DIR", help = "Use caliber-specific 5D BC correction tables")]
        bc_table_dir: Option<PathBuf>,

        /// Enable auto-download of BC5D correction tables when needed
        #[cfg(feature = "online")]
        #[arg(long, help = "Auto-download BC5D tables when needed")]
        bc_table_auto: bool,

        /// Base URL for BC5D table downloads
        #[cfg(feature = "online")]
        #[arg(long, default_value = "https://ballistics.tools/downloads/bc5d", help = "Base URL for BC5D table downloads")]
        bc_table_url: String,

        /// Force re-download of BC5D tables even if cached
        #[cfg(feature = "online")]
        #[arg(long, help = "Force re-download even if cached")]
        bc_table_refresh: bool,

        /// Bullet length in inches (for BC table lookup; estimated from diameter if not provided)
        #[arg(long)]
        bullet_length: Option<f64>,

        /// Barrel twist rate (inches per turn, e.g., 10 for 1:10)
        #[arg(long)]
        twist_rate: Option<f64>,

        /// Twist hand: `--twist-right` or `--twist-right true` = right-hand (default);
        /// `--twist-right false` = left-hand
        #[arg(
            long,
            num_args = 0..=1,
            default_value_t = true,
            default_missing_value = "true",
            action = clap::ArgAction::Set
        )]
        twist_right: bool,

        /// Latitude for Coriolis effect and weather zones (degrees, -90 to 90)
        #[arg(long, value_parser = f64_range(-90.0, 90.0), allow_hyphen_values = true)]
        latitude: Option<f64>,

        /// Longitude for weather zones (degrees, -180 to 180)
        #[arg(long, value_parser = f64_range(-180.0, 180.0), allow_hyphen_values = true)]
        longitude: Option<f64>,

        /// Shot direction/azimuth for Coriolis and weather zones (degrees, 0=North, 90=East)
        #[arg(long)]
        shot_direction: Option<f64>,

        /// Shooting angle (degrees, positive = uphill, negative = downhill)
        #[arg(long, default_value = "0.0", allow_hyphen_values = true)]
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

        // PDF Dope Card Parameters
        /// Target speed in mph (for lead calculation in PDF output)
        #[arg(long, default_value = "0.0")]
        target_speed: f64,

        /// Powder type/name (for PDF metadata)
        #[arg(long)]
        powder: Option<String>,

        /// Bullet name (for PDF metadata)
        #[arg(long)]
        bullet_name: Option<String>,

        /// Location name for PDF header (overrides --site for display)
        #[arg(long)]
        location_name: Option<String>,

        /// Output file path (required for PDF format)
        #[arg(long, value_name = "FILE")]
        output_file: Option<PathBuf>,

        /// Font scale factor for PDF output (1.0 = default, 1.5 = 50% larger)
        #[arg(long, default_value = "1.0")]
        font_scale: f32,

        /// Font size preset (small, medium, large) — overridden by --font-scale if both given
        #[arg(long, value_name = "SIZE")]
        font_preset: Option<String>,

        /// Use bold font for data cells in PDF output
        #[arg(long)]
        bold_data: bool,
    },

    /// Run Monte Carlo simulation
    MonteCarlo {
        /// Base velocity (m/s)
        #[arg(short = 'v', long, value_parser = f64_range(0.0, 6000.0))]
        velocity: f64,

        /// Launch angle (degrees)
        #[arg(short = 'a', long, default_value = "0.0")]
        angle: f64,

        /// Ballistic coefficient
        #[arg(short = 'b', long, value_parser = f64_range(0.001, 2.0))]
        bc: f64,

        /// Mass (kg)
        #[arg(short = 'm', long, value_parser = f64_range(0.1, 2000.0))]
        mass: f64,

        /// Diameter (meters)
        #[arg(short = 'd', long, value_parser = f64_range(0.01, 60.0))]
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
        #[arg(short = 'v', long, value_parser = f64_range(0.0, 6000.0))]
        velocity: f64,

        /// Ballistic coefficient
        #[arg(short = 'b', long, value_parser = f64_range(0.001, 2.0))]
        bc: f64,

        /// Mass (grains for imperial, grams for metric)
        #[arg(short = 'm', long, value_parser = f64_range(0.1, 2000.0))]
        mass: f64,

        /// Diameter (inches for imperial, mm for metric)
        #[arg(short = 'd', long, value_parser = f64_range(0.01, 60.0))]
        diameter: f64,

        /// Target distance (yards for imperial, meters for metric)
        #[arg(long)]
        target_distance: f64,

        /// Target height (yards for imperial, meters for metric)
        #[arg(long, default_value = "0.0", allow_hyphen_values = true)]
        target_height: f64,

        /// Sight height above bore (inches for imperial, mm for metric)
        #[arg(long)]
        sight_height: Option<f64>,

        /// Temperature (Fahrenheit or Celsius based on --units)
        #[arg(long, default_value = "59.0", value_parser = f64_range(-100.0, 200.0))]
        temperature: f64,

        /// Pressure (inHg or hPa based on --units)
        #[arg(long, default_value = "29.92", value_parser = f64_range(15.0, 1200.0))]
        pressure: f64,

        /// Humidity (0-100%)
        #[arg(long, default_value = "50.0", value_parser = f64_range(0.0, 100.0))]
        humidity: f64,

        /// Altitude (feet or meters based on --units)
        #[arg(long, default_value = "0.0")]
        altitude: f64,

        /// Output format
        #[arg(short = 'o', long, default_value = "table")]
        output: OutputFormat,
    },

    /// Estimate BC from trajectory data
    EstimateBC {
        /// Initial velocity (m/s)
        #[arg(short = 'v', long, value_parser = f64_range(0.0, 6000.0))]
        velocity: f64,

        /// Mass (kg)
        #[arg(short = 'm', long, value_parser = f64_range(0.1, 2000.0))]
        mass: f64,

        /// Diameter (meters)
        #[arg(short = 'd', long, value_parser = f64_range(0.01, 60.0))]
        diameter: f64,

        /// Distance 1 (meters)
        #[arg(long)]
        distance1: f64,

        /// Drop at distance 1 (meters)
        #[arg(long, allow_hyphen_values = true)]
        drop1: f64,

        /// Distance 2 (meters)
        #[arg(long)]
        distance2: f64,

        /// Drop at distance 2 (meters)
        #[arg(long, allow_hyphen_values = true)]
        drop2: f64,

        /// Output format
        #[arg(short = 'o', long, default_value = "table")]
        output: OutputFormat,
    },

    /// Generate BC segments for velocity-dependent BC
    GenerateBCSegments {
        /// Base ballistic coefficient
        #[arg(short = 'b', long, value_parser = f64_range(0.001, 2.0))]
        bc: f64,

        /// Projectile mass (kg)
        #[arg(short = 'm', long, value_parser = f64_range(0.1, 2000.0))]
        mass: f64,

        /// Projectile diameter (meters)
        #[arg(short = 'd', long, value_parser = f64_range(0.01, 60.0))]
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

    /// Calculate effective muzzle velocity from observed drop
    TrueVelocity {
        /// Measured drop in MILs at the target range
        #[arg(long)]
        measured_drop: f64,

        /// Range at which drop was measured (yards for imperial, meters for metric)
        #[arg(long)]
        range: f64,

        /// Ballistic coefficient
        #[arg(short = 'b', long, value_parser = f64_range(0.001, 2.0))]
        bc: f64,

        /// Drag model (G1, G7)
        #[arg(long, default_value = "g1")]
        drag_model: DragModelArg,

        /// Bullet weight (grains for imperial, grams for metric)
        #[arg(short = 'm', long, value_parser = f64_range(0.1, 2000.0))]
        mass: f64,

        /// Bullet diameter/caliber (inches for imperial, mm for metric)
        #[arg(short = 'd', long, value_parser = f64_range(0.01, 60.0))]
        diameter: f64,

        /// Chronograph velocity for comparison (fps for imperial, m/s for metric)
        #[arg(long)]
        chrono_velocity: Option<f64>,

        /// Zero distance (yards for imperial, meters for metric)
        #[arg(long, default_value = "100.0")]
        zero_distance: f64,

        /// Sight height above bore (inches for imperial, mm for metric)
        #[arg(long, default_value = "2.0")]
        sight_height: f64,

        /// Temperature (Fahrenheit for imperial, Celsius for metric)
        #[arg(long, default_value = "59.0", value_parser = f64_range(-100.0, 200.0))]
        temperature: f64,

        /// Pressure (inHg for imperial, hPa for metric)
        #[arg(long, default_value = "29.92", value_parser = f64_range(15.0, 1200.0))]
        pressure: f64,

        /// Humidity (0-100%)
        #[arg(long, default_value = "50.0", value_parser = f64_range(0.0, 100.0))]
        humidity: f64,

        /// Altitude (feet for imperial, meters for metric)
        #[arg(long, default_value = "0.0")]
        altitude: f64,

        /// Unit system
        #[arg(long, default_value = "imperial")]
        units: UnitSystem,

        /// Directory containing caliber-specific BC5D tables (e.g., bc5d_308.bin)
        #[arg(long, value_name = "DIR", help = "Use caliber-specific 5D BC correction tables")]
        bc_table_dir: Option<PathBuf>,

        /// Enable auto-download of BC5D correction tables when needed
        #[cfg(feature = "online")]
        #[arg(long, help = "Auto-download BC5D tables when needed")]
        bc_table_auto: bool,

        /// Base URL for BC5D table downloads
        #[cfg(feature = "online")]
        #[arg(long, default_value = "https://ballistics.tools/downloads/bc5d", help = "Base URL for BC5D table downloads")]
        bc_table_url: String,

        /// Force offline mode (use local calculation, skip API)
        #[arg(long, help = "Force offline mode using local calculation")]
        offline: bool,

        /// Fall back to local calculation if API fails
        #[cfg(feature = "online")]
        #[arg(long, help = "Use local calculation if API fails")]
        offline_fallback: bool,

        /// API endpoint URL
        #[cfg(feature = "online")]
        #[arg(long, default_value = "https://api.ballistics.7.62x51mm.sh")]
        api_url: String,

        /// API request timeout in seconds
        #[cfg(feature = "online")]
        #[arg(long, default_value = "10")]
        api_timeout: u64,

        /// Bullet length in inches (for BC table lookup; estimated from diameter if not provided)
        #[arg(long)]
        bullet_length: Option<f64>,

        /// Output format
        #[arg(short = 'o', long, default_value = "table")]
        output: OutputFormat,
    },

    /// Calculate Maximum Point-Blank Range (MPBR)
    Mpbr {
        /// Load parameters from saved profile
        #[arg(long, value_name = "NAME")]
        profile: Option<String>,

        /// Initial velocity (fps or m/s based on --units)
        #[arg(short = 'v', long, value_parser = f64_range(0.0, 6000.0))]
        velocity: Option<f64>,

        /// Ballistic coefficient
        #[arg(short = 'b', long, value_parser = f64_range(0.001, 2.0))]
        bc: Option<f64>,

        /// Mass (grains for imperial, grams for metric)
        #[arg(short = 'm', long, value_parser = f64_range(0.1, 2000.0))]
        mass: Option<f64>,

        /// Diameter (inches for imperial, mm for metric)
        #[arg(short = 'd', long, value_parser = f64_range(0.01, 60.0))]
        diameter: Option<f64>,

        /// Drag model (G1, G7)
        #[arg(long, default_value = "g1")]
        drag_model: DragModelArg,

        /// Vital zone diameter (inches for imperial, cm for metric)
        #[arg(long, default_value = "8.0")]
        vital_zone: f64,

        /// Sight height above bore (inches for imperial, mm for metric)
        #[arg(long)]
        sight_height: Option<f64>,

        /// Temperature (Fahrenheit or Celsius based on --units)
        #[arg(long, default_value = "59.0", value_parser = f64_range(-100.0, 200.0))]
        temperature: f64,

        /// Pressure (inHg or hPa based on --units)
        #[arg(long, default_value = "29.92", value_parser = f64_range(15.0, 1200.0))]
        pressure: f64,

        /// Humidity (0-100%)
        #[arg(long, default_value = "50.0", value_parser = f64_range(0.0, 100.0))]
        humidity: f64,

        /// Altitude (feet or meters based on --units)
        #[arg(long, default_value = "0.0")]
        altitude: f64,

        /// Output format
        #[arg(short = 'o', long, default_value = "table")]
        output: OutputFormat,
    },

    /// Generate come-up (elevation adjustment) table
    ComeUps {
        /// Load parameters from saved profile
        #[arg(long, value_name = "NAME")]
        profile: Option<String>,

        /// Initial velocity (fps or m/s based on --units)
        #[arg(short = 'v', long, value_parser = f64_range(0.0, 6000.0))]
        velocity: Option<f64>,

        /// Ballistic coefficient
        #[arg(short = 'b', long, value_parser = f64_range(0.001, 2.0))]
        bc: Option<f64>,

        /// Mass (grains for imperial, grams for metric)
        #[arg(short = 'm', long, value_parser = f64_range(0.1, 2000.0))]
        mass: Option<f64>,

        /// Diameter (inches for imperial, mm for metric)
        #[arg(short = 'd', long, value_parser = f64_range(0.01, 60.0))]
        diameter: Option<f64>,

        /// Drag model (G1, G7)
        #[arg(long, default_value = "g1")]
        drag_model: DragModelArg,

        /// Zero distance (yards for imperial, meters for metric)
        #[arg(long)]
        zero_distance: f64,

        /// Start range (yards or meters)
        #[arg(long, default_value = "100.0")]
        start: f64,

        /// End range (yards or meters)
        #[arg(long, default_value = "1200.0")]
        end: f64,

        /// Range step (yards or meters)
        #[arg(long, default_value = "50.0")]
        step: f64,

        /// Adjustment unit (mil or moa)
        #[arg(long, default_value = "mil")]
        adjustment_unit: AdjustmentUnit,

        /// Sight height above bore (inches for imperial, mm for metric)
        #[arg(long)]
        sight_height: Option<f64>,

        /// Temperature (Fahrenheit or Celsius based on --units)
        #[arg(long, default_value = "59.0", value_parser = f64_range(-100.0, 200.0))]
        temperature: f64,

        /// Pressure (inHg or hPa based on --units)
        #[arg(long, default_value = "29.92", value_parser = f64_range(15.0, 1200.0))]
        pressure: f64,

        /// Humidity (0-100%)
        #[arg(long, default_value = "50.0", value_parser = f64_range(0.0, 100.0))]
        humidity: f64,

        /// Altitude (feet or meters based on --units)
        #[arg(long, default_value = "0.0")]
        altitude: f64,

        /// Wind speed (mph or m/s based on --units)
        #[arg(long, default_value = "0.0")]
        wind_speed: f64,

        /// Wind direction (degrees, 0=headwind, 90=from right)
        #[arg(long, default_value = "0.0")]
        wind_direction: f64,

        /// Output format
        #[arg(short = 'o', long, default_value = "table")]
        output: OutputFormat,
    },

    /// Generate wind drift card (wind deflection at multiple speeds)
    WindCard {
        /// Load parameters from saved profile
        #[arg(long, value_name = "NAME")]
        profile: Option<String>,

        /// Initial velocity (fps or m/s based on --units)
        #[arg(short = 'v', long, value_parser = f64_range(0.0, 6000.0))]
        velocity: Option<f64>,

        /// Ballistic coefficient
        #[arg(short = 'b', long, value_parser = f64_range(0.001, 2.0))]
        bc: Option<f64>,

        /// Mass (grains for imperial, grams for metric)
        #[arg(short = 'm', long, value_parser = f64_range(0.1, 2000.0))]
        mass: Option<f64>,

        /// Diameter (inches for imperial, mm for metric)
        #[arg(short = 'd', long, value_parser = f64_range(0.01, 60.0))]
        diameter: Option<f64>,

        /// Drag model (G1, G7)
        #[arg(long, default_value = "g1")]
        drag_model: DragModelArg,

        /// Zero distance (yards for imperial, meters for metric)
        #[arg(long)]
        zero_distance: f64,

        /// Comma-separated wind speeds to calculate (mph or m/s)
        #[arg(long, default_value = "5,10,15,20")]
        wind_speeds: String,

        /// Start range (yards or meters)
        #[arg(long, default_value = "100.0")]
        start: f64,

        /// End range (yards or meters)
        #[arg(long, default_value = "1000.0")]
        end: f64,

        /// Range step (yards or meters)
        #[arg(long, default_value = "100.0")]
        step: f64,

        /// Adjustment unit (mil or moa)
        #[arg(long, default_value = "mil")]
        adjustment_unit: AdjustmentUnit,

        /// Sight height above bore (inches for imperial, mm for metric)
        #[arg(long)]
        sight_height: Option<f64>,

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
    },

    /// Analyze gyroscopic stability (Miller stability factor)
    Stability {
        /// Load parameters from saved profile
        #[arg(long, value_name = "NAME")]
        profile: Option<String>,

        /// Bullet mass (grains for imperial, grams for metric)
        #[arg(short = 'm', long)]
        mass: Option<f64>,

        /// Bullet diameter (inches for imperial, mm for metric)
        #[arg(short = 'd', long)]
        diameter: Option<f64>,

        /// Bullet length (inches for imperial, mm for metric)
        #[arg(short = 'l', long)]
        length: f64,

        /// Barrel twist rate (inches/turn for imperial, mm/turn for metric)
        #[arg(short = 't', long)]
        twist_rate: f64,

        /// Muzzle velocity (fps for imperial, m/s for metric; default: 2700 fps / 823 m/s)
        #[arg(short = 'v', long)]
        velocity: Option<f64>,

        /// Temperature (Fahrenheit or Celsius based on --units)
        #[arg(long, default_value = "59.0")]
        temperature: f64,

        /// Pressure (inHg or hPa based on --units)
        #[arg(long, default_value = "29.92")]
        pressure: f64,

        /// Altitude (feet or meters based on --units)
        #[arg(long, default_value = "0.0")]
        altitude: f64,

        /// Output format
        #[arg(short = 'o', long, default_value = "table")]
        output: OutputFormat,
    },

    /// Generate comprehensive range table (drop, wind, velocity, energy, ToF)
    RangeTable {
        /// Load parameters from saved profile
        #[arg(long, value_name = "NAME")]
        profile: Option<String>,

        /// Initial velocity (fps or m/s based on --units)
        #[arg(short = 'v', long)]
        velocity: Option<f64>,

        /// Ballistic coefficient
        #[arg(short = 'b', long)]
        bc: Option<f64>,

        /// Mass (grains for imperial, grams for metric)
        #[arg(short = 'm', long)]
        mass: Option<f64>,

        /// Diameter (inches for imperial, mm for metric)
        #[arg(short = 'd', long)]
        diameter: Option<f64>,

        /// Drag model (G1, G7)
        #[arg(long, default_value = "g1")]
        drag_model: DragModelArg,

        /// Zero distance (yards for imperial, meters for metric)
        #[arg(long)]
        zero_distance: f64,

        /// Start range (yards or meters)
        #[arg(long, default_value = "100.0")]
        start: f64,

        /// End range (yards or meters)
        #[arg(long, default_value = "1200.0")]
        end: f64,

        /// Range step (yards or meters)
        #[arg(long, default_value = "50.0")]
        step: f64,

        /// Wind speed (mph or m/s based on --units)
        #[arg(long, default_value = "10.0")]
        wind_speed: f64,

        /// Wind direction (degrees, 0=headwind, 90=from right)
        #[arg(long, default_value = "90.0")]
        wind_direction: f64,

        /// Adjustment unit (mil or moa)
        #[arg(long, default_value = "mil")]
        adjustment_unit: AdjustmentUnit,

        /// Sight height above bore (inches for imperial, mm for metric)
        #[arg(long)]
        sight_height: Option<f64>,

        /// Temperature (Fahrenheit or Celsius based on --units)
        #[arg(long, default_value = "59.0", value_parser = f64_range(-100.0, 200.0))]
        temperature: f64,

        /// Pressure (inHg or hPa based on --units)
        #[arg(long, default_value = "29.92", value_parser = f64_range(15.0, 1200.0))]
        pressure: f64,

        /// Humidity (0-100%)
        #[arg(long, default_value = "50.0", value_parser = f64_range(0.0, 100.0))]
        humidity: f64,

        /// Altitude (feet or meters based on --units)
        #[arg(long, default_value = "0.0")]
        altitude: f64,

        /// Output format
        #[arg(short = 'o', long, default_value = "table")]
        output: OutputFormat,
    },

    /// Manage saved ballistic profiles
    Profile {
        #[command(subcommand)]
        action: ProfileAction,
    },

    /// Generate shell completions
    Completions {
        /// Shell to generate completions for
        #[arg(value_enum)]
        shell: clap_complete::Shell,
    },
}

#[derive(Subcommand)]
enum ProfileAction {
    /// Save a new profile
    Save {
        /// Profile name (used for recall)
        name: String,

        /// Initial velocity (fps or m/s based on --units)
        #[arg(short = 'v', long, value_parser = f64_range(0.0, 6000.0))]
        velocity: f64,

        /// Ballistic coefficient
        #[arg(short = 'b', long, value_parser = f64_range(0.001, 2.0))]
        bc: f64,

        /// Mass (grains for imperial, grams for metric)
        #[arg(short = 'm', long, value_parser = f64_range(0.1, 2000.0))]
        mass: f64,

        /// Diameter (inches for imperial, mm for metric)
        #[arg(short = 'd', long, value_parser = f64_range(0.01, 60.0))]
        diameter: f64,

        /// Drag model (G1, G7)
        #[arg(long, default_value = "g1")]
        drag_model: DragModelArg,

        /// Barrel twist rate (inches per turn)
        #[arg(long)]
        twist_rate: Option<f64>,

        /// Sight height above bore (inches for imperial, mm for metric)
        #[arg(long)]
        sight_height: Option<f64>,

        /// Default zero distance (yards for imperial, meters for metric)
        #[arg(long)]
        zero_distance: Option<f64>,

        /// Default temperature
        #[arg(long, default_value = "59.0", value_parser = f64_range(-100.0, 200.0))]
        temperature: f64,

        /// Default pressure
        #[arg(long, default_value = "29.92", value_parser = f64_range(15.0, 1200.0))]
        pressure: f64,

        /// Default humidity
        #[arg(long, default_value = "50.0", value_parser = f64_range(0.0, 100.0))]
        humidity: f64,

        /// Default altitude
        #[arg(long, default_value = "0.0")]
        altitude: f64,

        /// Bullet name/description
        #[arg(long)]
        bullet_name: Option<String>,

        /// Wind speed (mph or m/s based on --units)
        #[arg(long)]
        wind_speed: Option<f64>,

        /// Wind direction (degrees, 0=headwind, 90=from right)
        #[arg(long)]
        wind_direction: Option<f64>,

        /// Shooting angle (degrees, positive = uphill, negative = downhill)
        #[arg(long, allow_hyphen_values = true)]
        shooting_angle: Option<f64>,

        /// Auto-zero / zero distance for trajectory
        #[arg(long)]
        auto_zero: Option<f64>,

        /// Right-hand twist (default: true)
        #[arg(long)]
        twist_right: Option<bool>,

        /// Enable velocity-based BC segmentation
        #[arg(long)]
        use_bc_segments: Option<bool>,

        /// Bullet length in inches (for BC table lookup)
        #[arg(long)]
        bullet_length: Option<f64>,
    },

    /// List all saved profiles
    List,

    /// Show details of a saved profile
    Show {
        /// Profile name
        name: String,
    },

    /// Delete a saved profile
    Delete {
        /// Profile name
        name: String,
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
    Pdf,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum MonteCarloOutput {
    Summary,
    Full,
    Statistics,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum AdjustmentUnit {
    /// Milliradians (1 MIL = 3.6 inches at 100 yards)
    Mil,
    /// Minutes of Angle (1 MOA = 1.047 inches at 100 yards)
    Moa,
}

/// Saved ballistic profile for quick recall
#[derive(Debug, Clone, Serialize, Deserialize)]
struct ProfileData {
    name: String,
    velocity: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    drag_model: String,
    #[serde(default)]
    twist_rate: Option<f64>,
    #[serde(default)]
    sight_height: Option<f64>,
    #[serde(default)]
    zero_distance: Option<f64>,
    #[serde(default = "default_unit_system")]
    units: String,
    #[serde(default = "default_temperature")]
    temperature: f64,
    #[serde(default = "default_pressure")]
    pressure: f64,
    #[serde(default = "default_humidity")]
    humidity: f64,
    #[serde(default)]
    altitude: f64,
    #[serde(default)]
    bullet_name: Option<String>,
    #[serde(default)]
    created: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    wind_speed: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    wind_direction: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    shooting_angle: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    auto_zero: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    twist_right: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    use_bc_segments: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    bullet_length: Option<f64>,
}

fn default_unit_system() -> String { "imperial".to_string() }
fn default_temperature() -> f64 { 59.0 }
fn default_pressure() -> f64 { 29.92 }
fn default_humidity() -> f64 { 50.0 }

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

/// Configuration struct for run_trajectory(), replacing the 31+ parameter signature.
/// All distance/velocity/mass/temperature values should be in metric units
/// (converted before construction), except where noted.
struct TrajectoryConfig {
    // Core ballistic params (metric)
    velocity: f64,
    angle: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    drag_model: DragModelArg,
    max_range: f64,
    time_step: f64,

    // Environment (metric)
    temperature: f64,
    pressure: f64,
    humidity: f64,
    altitude: f64,

    // Wind (metric)
    wind_speed: f64,
    wind_direction: f64,
    // Downrange-segmented wind (engine units: speed km/h, angle deg, distance m).
    // When non-empty, overrides the scalar wind above.
    wind_segments: Vec<ballistics_engine::wind::WindSegment>,

    // Output
    output: OutputFormat,
    full: bool,
    units: UnitSystem,

    // Sight / bore geometry (metric)
    sight_height: f64,
    bore_height: f64,
    ignore_ground_impact: bool,

    // BC options
    use_bc_segments: bool,
    use_cluster_bc: bool,
    bc_table_segments: Option<Vec<BCSegmentData>>,

    // Advanced physics toggles
    enable_magnus: bool,
    enable_coriolis: bool,
    enable_spin_drift: bool,
    enable_aerodynamic_jump: bool,
    enable_wind_shear: bool,
    enable_pitch_damping: bool,
    enable_precession: bool,

    // Sampling
    sample_trajectory: bool,
    sample_interval: f64,

    // Integration method
    use_rk4: bool,
    use_rk45: bool,

    // Bullet / rifling
    twist_rate: Option<f64>,
    twist_right: bool,
    latitude: Option<f64>,
    shooting_angle: f64,

    // Powder sensitivity
    use_powder_sensitivity: bool,
    powder_temp_sensitivity: f64,
    powder_temp: f64,

    // PDF metadata
    pdf_metadata: Option<PdfMetadata>,
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

// ============================================================================
// Fuzzy Column Name Matching (MBA-614)
// ============================================================================

/// Known column names for gun profiles and locations
const KNOWN_COLUMNS: &[&str] = &[
    // Gun profile columns
    "RIFLE_NAME", "BULLET_NAME", "CALIBER", "BULLET_LENGTH", "VELOCITY",
    "MUZZLE_VELOCITY", "MV", "ZERO_TEMP", "ZERO_ALT", "ZERO_RANGE",
    "ZERO_DISTANCE", "ZERO", "BC", "BC_TYPE", "BC_ADJ", "BC_ADJUSTMENT",
    "BULLET_WEIGHT", "WIND_SPEED", "TARGET_SPEED", "WIND_ANGLE",
    "SIGHT_HEIGHT", "BARREL_TWIST", "TWIST_RATE", "TWIST", "VELOCITY_ADJ",
    "VEL_ADJ", "RANGE_MIN", "RANGE_MAX", "RANGE_INCR", "GUN_MPH", "COA",
    "POWDER_NAME", "POWDER_CHARGE", "PLASTIC_TIP_LENGTH", "JBM_BULLET_ID",
    "SPIN_DRIFT", "CHART_TYPE", "REGENERATE", "V_OFFSET_MIL", "H_OFFSET_MIL",
    "LATITUDE", "LONGITUDE", "DRAG_MODEL",
    // Location columns
    "LOCATION_NAME", "ALTITUDE", "ALT", "PRESSURE", "PRESSURE(HPA OR INHG)",
    "TARGET_TEMP", "TEMPERATURE", "TEMP", "HUMIDITY", "DA", "DENSITY_ALTITUDE",
    "WIND_DIR", "WIND_DIRECTION",
];

/// Find the closest matching known column name using Levenshtein distance
/// Returns Some((known_name, distance)) if a match within max_distance is found
fn find_closest_column(unknown: &str, max_distance: usize) -> Option<(&'static str, usize)> {
    let unknown_upper = unknown.to_uppercase();

    // First check for exact match
    for &known in KNOWN_COLUMNS {
        if known == unknown_upper {
            return Some((known, 0));
        }
    }

    // Find closest match within max_distance
    KNOWN_COLUMNS
        .iter()
        .filter_map(|&known| {
            let dist = levenshtein(&unknown_upper, known);
            if dist <= max_distance && dist > 0 {
                Some((known, dist))
            } else {
                None
            }
        })
        .min_by_key(|(_, dist)| *dist)
}

/// Normalize a column header, suggesting corrections for typos
/// Returns (normalized_name, was_corrected, original_name)
fn normalize_column_header(header: &str) -> (String, bool, String) {
    let upper = header.trim().to_uppercase();

    // Check if it's a known column or close match
    if let Some((known, dist)) = find_closest_column(&upper, 2) {
        if dist == 0 {
            // Exact match
            (known.to_string(), false, header.to_string())
        } else {
            // Fuzzy match - return corrected name
            (known.to_string(), true, header.to_string())
        }
    } else {
        // Unknown column, keep as-is
        (upper, false, header.to_string())
    }
}

/// Parse a CSV file and return row matching the first column value
/// Parse a CSV file and return row matching the first column value.
/// Uses the `csv` crate for proper quoting/escaping support (MBA-730).
fn load_csv_row(path: &PathBuf, row_name: &str) -> Result<HashMap<String, String>, Box<dyn Error>> {
    use csv::ReaderBuilder;

    let mut reader = ReaderBuilder::new()
        .has_headers(true)
        .comment(Some(b'#'))
        .flexible(true)
        .from_path(path)?;

    // Normalize headers and track fuzzy corrections
    let mut corrections: Vec<(String, String)> = Vec::new();
    let headers: Vec<String> = reader.headers()?
        .iter()
        .map(|h| {
            let (normalized, was_corrected, original) = normalize_column_header(h);
            if was_corrected {
                corrections.push((original, normalized.clone()));
            }
            normalized
        })
        .collect();

    // Print warnings for any corrected column names
    for (original, corrected) in &corrections {
        eprintln!("Warning: Column '{}' not recognized - did you mean '{}'? Using {}.",
                  original.trim(), corrected, corrected);
    }

    if headers.is_empty() {
        return Err("CSV file has no header row".into());
    }

    // Parse data rows and find matching row
    for record in reader.records() {
        let record = record?;
        if let Some(name) = record.get(0) {
            if name.eq_ignore_ascii_case(row_name) {
                let mut map = HashMap::new();
                for (i, header) in headers.iter().enumerate() {
                    if let Some(val) = record.get(i) {
                        if !val.is_empty() {
                            map.insert(header.clone(), val.to_string());
                        }
                    }
                }
                return Ok(map);
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

// ============================================================================
// Profile Management Helpers
// ============================================================================

/// Get the profiles directory (~/.ballistics/profiles/)
fn get_profiles_dir() -> Result<PathBuf, Box<dyn Error>> {
    let home = dirs::home_dir().ok_or("Could not determine home directory")?;
    let dir = home.join(".ballistics").join("profiles");
    fs::create_dir_all(&dir)?;
    Ok(dir)
}

/// Save a profile to disk
fn save_profile(profile: &ProfileData) -> Result<PathBuf, Box<dyn Error>> {
    let dir = get_profiles_dir()?;
    let path = dir.join(format!("{}.json", profile.name));
    let json = serde_json::to_string_pretty(profile)?;
    fs::write(&path, json)?;
    Ok(path)
}

/// Load a profile by name
fn load_profile(name: &str) -> Result<ProfileData, Box<dyn Error>> {
    let dir = get_profiles_dir()?;
    let path = dir.join(format!("{}.json", name));
    if !path.exists() {
        return Err(format!("Profile '{}' not found. Use 'ballistics profile list' to see available profiles.", name).into());
    }
    let content = fs::read_to_string(&path)?;
    let profile: ProfileData = serde_json::from_str(&content)?;
    Ok(profile)
}

/// List all saved profiles
fn list_profiles() -> Result<Vec<ProfileData>, Box<dyn Error>> {
    let dir = get_profiles_dir()?;
    let mut profiles = Vec::new();
    for entry in fs::read_dir(&dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().and_then(|e| e.to_str()) == Some("json") {
            if let Ok(content) = fs::read_to_string(&path) {
                if let Ok(profile) = serde_json::from_str::<ProfileData>(&content) {
                    profiles.push(profile);
                }
            }
        }
    }
    profiles.sort_by(|a, b| a.name.cmp(&b.name));
    Ok(profiles)
}

/// Delete a profile by name
fn delete_profile(name: &str) -> Result<(), Box<dyn Error>> {
    let dir = get_profiles_dir()?;
    let path = dir.join(format!("{}.json", name));
    if !path.exists() {
        return Err(format!("Profile '{}' not found.", name).into());
    }
    fs::remove_file(&path)?;
    Ok(())
}

/// Convert drop to adjustment unit (MIL or MOA)
fn drop_to_adjustment(drop_yd: f64, range_yd: f64, unit: AdjustmentUnit) -> f64 {
    if range_yd < 1.0 {
        return 0.0;
    }
    match unit {
        AdjustmentUnit::Mil => (drop_yd / range_yd) * 1000.0,
        AdjustmentUnit::Moa => (drop_yd / range_yd) * 3438.0,
    }
}

/// Get a timestamp string without chrono
fn timestamp_string() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    format!("{}", secs)
}

/// Parse the drag model string from a profile
fn parse_drag_model_arg(s: &str) -> DragModelArg {
    match s.to_uppercase().as_str() {
        "G7" => DragModelArg::G7,
        _ => DragModelArg::G1,
    }
}

/// Parse a `--wind-segment` value `"SPEED:ANGLE:UNTIL_DISTANCE"` into an engine
/// `WindSegment` `(speed_kmh, angle_deg, until_distance_m)`. SPEED and UNTIL_DISTANCE
/// are interpreted in the CLI display units (mph & yards imperial, m/s & meters
/// metric); ANGLE is degrees in the wind-FROM convention (same as `--wind-direction`).
fn parse_wind_segment(
    s: &str,
    units: UnitSystem,
) -> Result<ballistics_engine::wind::WindSegment, String> {
    ballistics_engine::wind::parse_wind_segment_str(s, matches!(units, UnitSystem::Imperial))
}

fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Trajectory {
            profile,
            profile_row,
            saved_profile,
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
            wind_segment,
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
            enable_aerodynamic_jump,
            enable_wind_shear,
            sample_trajectory,
            sample_interval,
            enable_pitch_damping,
            enable_precession,
            use_cluster_bc,
            bc_table,
            bc_table_dir,
            #[cfg(feature = "online")]
            bc_table_auto,
            #[cfg(feature = "online")]
            bc_table_url,
            #[cfg(feature = "online")]
            bc_table_refresh,
            bullet_length,
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
            // PDF dope card parameters
            target_speed,
            powder,
            bullet_name,
            location_name,
            output_file,
            font_scale,
            font_preset,
            bold_data,
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

            // Load saved JSON profile if specified
            let saved_profile_data: Option<ProfileData> = saved_profile.as_ref().map(|name| {
                load_profile(name).unwrap_or_else(|e| { eprintln!("Error loading saved profile: {}", e); std::process::exit(1); })
            });
            if let Some(ref sp) = saved_profile_data { eprintln!("Loaded saved profile '{}'", sp.name); }

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
            let final_velocity = velocity
                .or_else(|| { let v = csv_get_f64(&profile_data, &["VELOCITY", "MV", "MUZZLE_VELOCITY"], 0.0); if v > 0.0 { Some(v) } else { None } })
                .or_else(|| saved_profile_data.as_ref().map(|p| p.velocity))
                .unwrap_or(2800.0);
            let final_velocity_adj = if velocity_adjustment != 0.0 { velocity_adjustment } else { csv_get_f64(&profile_data, &["VELOCITY_ADJ", "VEL_ADJ"], 0.0) };
            let final_bc = bc
                .or_else(|| { let v = csv_get_f64(&profile_data, &["BC"], 0.0); if v > 0.0 { Some(v) } else { None } })
                .or_else(|| saved_profile_data.as_ref().map(|p| p.bc))
                .unwrap_or(0.5);
            let final_bc_adj = if bc_adjustment != 1.0 { bc_adjustment } else { csv_get_f64(&profile_data, &["BC_ADJ", "BC_ADJUSTMENT"], 1.0) };
            let final_mass = mass
                .or_else(|| saved_profile_data.as_ref().map(|p| p.mass))
                .unwrap_or_else(|| { eprintln!("Error: --mass is required (or use --saved-profile)"); std::process::exit(1); });
            let final_diameter = diameter
                .or_else(|| saved_profile_data.as_ref().map(|p| p.diameter))
                .unwrap_or_else(|| { eprintln!("Error: --diameter is required (or use --saved-profile)"); std::process::exit(1); });
            let final_max_range = max_range;
            let final_wind_speed = if wind_speed != 0.0 { wind_speed } else { saved_profile_data.as_ref().and_then(|p| p.wind_speed).unwrap_or(0.0) };
            let final_wind_direction = if wind_direction != 0.0 { wind_direction } else { csv_get_f64(&location_data, &["WIND_DIR", "WIND_DIRECTION"], 0.0) };

            // Location overrides (environmental conditions)
            let final_temperature = if temperature != 59.0 { temperature } else { csv_get_f64(&location_data, &["TARGET_TEMP", "TEMPERATURE", "TEMP"], csv_get_f64(&profile_data, &["ZERO_TEMP"], 59.0)) };
            let final_pressure = if pressure != 29.92 { pressure } else { csv_get_f64(&location_data, &["PRESSURE", "PRESSURE(HPA OR INHG)"], 29.92) };
            let final_humidity = if humidity != 50.0 { humidity } else { csv_get_f64(&location_data, &["HUMIDITY"], 50.0) };
            let final_altitude = if altitude != 0.0 { altitude } else { csv_get_f64(&location_data, &["ALTITUDE", "ALT"], csv_get_f64(&profile_data, &["ZERO_ALT"], 0.0)) };

            // Get zero range: CLI --auto-zero overrides profile ZERO_RANGE
            let final_auto_zero: Option<f64> = auto_zero.or_else(|| {
                let zero_from_csv = csv_get_f64(&profile_data, &["ZERO_RANGE", "ZERO_DISTANCE", "ZERO"], 0.0);
                if zero_from_csv > 0.0 { Some(zero_from_csv) } else { saved_profile_data.as_ref().and_then(|p| p.auto_zero.or(p.zero_distance)) }
            });

            // Resolve additional params from saved profile (if not explicitly set via CLI)
            let drag_model = if saved_profile_data.is_some() && velocity.is_none() && bc.is_none() {
                saved_profile_data.as_ref().map(|p| parse_drag_model_arg(&p.drag_model)).unwrap_or(drag_model)
            } else { drag_model };
            let use_bc_segments = use_bc_segments || saved_profile_data.as_ref().and_then(|p| p.use_bc_segments).unwrap_or(false);
            let twist_right = saved_profile_data.as_ref().and_then(|p| p.twist_right).unwrap_or(twist_right);
            let shooting_angle = if shooting_angle != 0.0 { shooting_angle } else { saved_profile_data.as_ref().and_then(|p| p.shooting_angle).unwrap_or(0.0) };
            let bullet_length = bullet_length.or_else(|| saved_profile_data.as_ref().and_then(|p| p.bullet_length));
            let sight_height = sight_height.or_else(|| saved_profile_data.as_ref().and_then(|p| p.sight_height));
            let twist_rate = twist_rate.or_else(|| saved_profile_data.as_ref().and_then(|p| p.twist_rate));

            // Apply truing adjustments
            let trued_velocity = final_velocity + final_velocity_adj;
            let mut trued_bc = final_bc * final_bc_adj;

            // Apply BC correction from table if provided
            // Generate velocity-dependent BC segments from the table for accurate trajectory
            let mut bc_table_segments: Option<Vec<BCSegmentData>> = None;
            let bc_table_correction: Option<f64> = if let Some(table_path) = &bc_table {
                match BcCorrectionTable::load(table_path) {
                    Ok(table) => {
                        // Get bullet length: CLI arg > CSV profile > estimate from diameter
                        let bullet_length_in = bullet_length
                            .or_else(|| {
                                let len = csv_get_f64(&profile_data, &["BULLET_LENGTH", "LENGTH"], 0.0);
                                if len > 0.0 { Some(len) } else { None }
                            })
                            .unwrap_or_else(|| {
                                // Estimate: typical rifle bullets are ~3.5 calibers long
                                final_diameter * 3.5
                            });

                        // Get bullet mass in grains
                        let mass_grains = match cli.units {
                            UnitSystem::Imperial => final_mass, // already in grains
                            UnitSystem::Metric => final_mass * 15.4324, // grams to grains
                        };

                        // BC type string
                        let bc_type_str = match drag_model {
                            DragModelArg::G1 => "G1",
                            DragModelArg::G7 => "G7",
                        };

                        // Velocity breakpoints for BC segments (denser in transonic region)
                        // These match the table's velocity bins for optimal interpolation
                        let velocity_breakpoints: Vec<f64> = vec![
                            trued_velocity, // Muzzle velocity
                            2700.0, 2500.0, 2300.0, 2100.0, 2000.0, 1900.0, 1800.0,
                            1700.0, 1600.0, 1500.0, 1400.0, 1350.0, 1300.0, 1250.0,
                            1200.0, 1150.0, 1100.0, 1050.0, 1000.0, 950.0, 900.0,
                            850.0, 800.0, 700.0, 600.0, 500.0,
                        ];

                        // Filter breakpoints to be below muzzle velocity and sort descending
                        let mut valid_velocities: Vec<f64> = velocity_breakpoints
                            .into_iter()
                            .filter(|&v| v <= trued_velocity && v >= 500.0)
                            .collect();
                        valid_velocities.sort_by(|a, b| b.partial_cmp(a).unwrap());
                        valid_velocities.dedup();

                        // Generate BC segments from table
                        let mut segments: Vec<BCSegmentData> = Vec::new();
                        for i in 0..valid_velocities.len().saturating_sub(1) {
                            let vel_max = valid_velocities[i];
                            let vel_min = valid_velocities[i + 1];
                            let vel_mid = (vel_max + vel_min) / 2.0;

                            // Lookup correction at midpoint velocity
                            let correction = table.lookup(
                                final_bc,
                                bc_type_str,
                                mass_grains,
                                bullet_length_in,
                                vel_mid,
                            );

                            // Calculate corrected BC for this segment
                            let segment_bc = final_bc * correction;

                            segments.push(BCSegmentData {
                                velocity_min: vel_min,
                                velocity_max: vel_max,
                                bc_value: segment_bc,
                            });
                        }

                        // Get correction at muzzle velocity for display
                        let muzzle_correction = table.lookup(
                            final_bc,
                            bc_type_str,
                            mass_grains,
                            bullet_length_in,
                            trued_velocity,
                        );

                        eprintln!("BC Table: Loaded {} (v{}, {})",
                            table_path.display(), table.version(), table.dimensions_str());
                        eprintln!("BC Table: Generated {} velocity-dependent BC segments", segments.len());
                        eprintln!("BC Table: Muzzle correction={:.4} for BC={:.3} {} {}gr {:.3}\" @ {:.0}fps",
                            muzzle_correction, final_bc, bc_type_str, mass_grains, bullet_length_in, trued_velocity);

                        // Show BC range if segments were generated
                        if !segments.is_empty() {
                            let min_bc = segments.iter().map(|s| s.bc_value).fold(f64::INFINITY, f64::min);
                            let max_bc = segments.iter().map(|s| s.bc_value).fold(f64::NEG_INFINITY, f64::max);
                            eprintln!("BC Table: BC range {:.5} - {:.5} across velocity envelope", min_bc, max_bc);
                            bc_table_segments = Some(segments);
                        }

                        // Still apply correction to trued_bc for display purposes
                        trued_bc *= muzzle_correction;
                        Some(muzzle_correction)
                    }
                    Err(e) => {
                        eprintln!("Warning: Failed to load BC table: {}", e);
                        eprintln!("         Continuing without BC correction table.");
                        None
                    }
                }
            } else {
                None
            };

            // Apply BC correction from 5D caliber-specific table if provided (and bc_table wasn't used)
            // Determine effective table directory: explicit --bc-table-dir, or auto-download cache
            #[cfg(feature = "online")]
            let effective_bc_table_dir: Option<PathBuf> = if bc_table_dir.is_some() {
                bc_table_dir.clone()
            } else if bc_table_auto {
                // Auto-download mode: ensure table is available, use cache directory
                match Bc5dDownloader::new(&bc_table_url, bc_table_refresh) {
                    Ok(mut downloader) => {
                        match downloader.ensure_table(final_diameter) {
                            Ok(table_path) => {
                                // Return the parent directory (cache dir) as table directory
                                table_path.parent().map(|p| p.to_path_buf())
                            }
                            Err(e) => {
                                eprintln!("Warning: {}", e);
                                eprintln!("         Continuing without BC5D correction table.");
                                None
                            }
                        }
                    }
                    Err(e) => {
                        eprintln!("Warning: BC5D table download failed: {}", e);
                        eprintln!("         Continuing without BC5D correction table.");
                        None
                    }
                }
            } else {
                None
            };

            #[cfg(not(feature = "online"))]
            let effective_bc_table_dir: Option<PathBuf> = bc_table_dir.clone();

            // Apply BC correction from 5D caliber-specific table if provided (MBA-744: uses shared helper)
            let bc_table_5d_correction: Option<f64> = if bc_table_correction.is_none() {
                if let Some(table_dir) = &effective_bc_table_dir {
                    let mut manager = Bc5dTableManager::new(table_dir);

                    let mass_grains = match cli.units {
                        UnitSystem::Imperial => final_mass,
                        UnitSystem::Metric => final_mass * 15.4324,
                    };

                    let bc_type_str = match drag_model {
                        DragModelArg::G1 => "G1",
                        DragModelArg::G7 => "G7",
                    };

                    // Print table info if available
                    if let Ok(table) = manager.get_table(final_diameter) {
                        eprintln!("BC5D Table: Loaded caliber {:.3} (v{}, API {}, {})",
                            table.caliber(), table.version(), table.api_version(), table.dimensions_str());
                    }

                    let bullet_length_in = bullet_length
                        .or_else(|| {
                            let csv_len = csv_get_f64(&profile_data, &["BULLET_LENGTH", "LENGTH"], 0.0);
                            if csv_len > 0.0 { Some(csv_len) } else { None }
                        })
                        .unwrap_or(final_diameter * 3.5);

                    if let Some(segments) = generate_bc5d_segments(
                        &mut manager,
                        final_diameter,
                        final_bc,
                        bc_type_str,
                        mass_grains,
                        Some(trued_velocity),
                        Some(bullet_length_in),
                    ) {
                        bc_table_segments = Some(segments);

                        let muzzle_correction = manager.lookup(
                            final_diameter, mass_grains, final_bc, trued_velocity, trued_velocity, bc_type_str
                        ).unwrap_or(1.0);
                        eprintln!("BC5D Table: Muzzle correction={:.4} for BC={:.3} {} {}gr @ {:.0}fps",
                            muzzle_correction, final_bc, bc_type_str, mass_grains, trued_velocity);
                        trued_bc *= muzzle_correction;
                        Some(muzzle_correction)
                    } else {
                        None
                    }
                } else {
                    None
                }
            } else {
                None // bc_table was already used
            };

            // Combine corrections for display
            let combined_bc_correction = bc_table_correction.or(bc_table_5d_correction);

            // Show effective values if using profile/location or BC table
            if !profile_data.is_empty() || !location_data.is_empty() || saved_profile_data.is_some() || combined_bc_correction.is_some() {
                let bc_info = if let Some(corr) = combined_bc_correction {
                    format!("BC={:.3} (table-corrected={:.4}, factor={:.4})", final_bc, trued_bc, corr)
                } else {
                    format!("BC={:.3} (trued={:.4})", final_bc, trued_bc)
                };
                eprintln!("Effective values: velocity={:.1} (trued={:.1}), {}",
                    final_velocity, trued_velocity, bc_info);
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

            // Construct PDF metadata if PDF output is requested
            let pdf_metadata = if matches!(output, OutputFormat::Pdf) {
                // Get rifle name from profile_row or use a default
                let rifle_name = profile_row.clone().unwrap_or_else(|| "Rifle".to_string());
                // Get location name: CLI > site > default
                let loc_name = location_name.clone()
                    .or_else(|| site.clone())
                    .unwrap_or_else(|| "Field".to_string());
                // Get powder and bullet from CLI, or from profile if available
                let powder_name = powder.clone()
                    .or_else(|| profile_data.get("POWDER_NAME").cloned())
                    .unwrap_or_default();
                let bullet_display = bullet_name.clone()
                    .or_else(|| profile_data.get("BULLET_NAME").cloned())
                    .unwrap_or_default();

                // Resolve effective font scale: --font-scale overrides --font-preset
                #[cfg(feature = "pdf")]
                let effective_font_scale = if font_scale != 1.0 {
                    font_scale
                } else if let Some(ref preset) = font_preset {
                    FontSizePreset::from_str(preset).map(|p| p.scale()).unwrap_or_else(|| {
                        eprintln!("Warning: Unknown font preset '{}', using medium", preset);
                        1.0
                    })
                } else {
                    1.0
                };
                #[cfg(not(feature = "pdf"))]
                let effective_font_scale = font_scale;

                Some(PdfMetadata {
                    rifle_name,
                    location_name: loc_name,
                    powder: powder_name,
                    bullet_name: bullet_display,
                    target_speed_mph: target_speed,
                    output_file: output_file.clone(),
                    velocity_fps: trued_velocity,
                    temperature_f: final_temperature,
                    pressure_inhg: final_pressure,
                    altitude_ft: final_altitude,
                    wind_speed_mph: final_wind_speed,
                    weight_gr: bullet_mass,
                    font_scale: effective_font_scale,
                    bold_data,
                })
            } else {
                None
            };

            // Calculate zero angle if auto-zero is specified (from CLI or profile)
            let muzzle_angle = if let Some(zero_distance) = final_auto_zero {
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


            // Parse downrange wind segments (display units -> engine units).
            let wind_segments: Vec<ballistics_engine::wind::WindSegment> = wind_segment
                .iter()
                .map(|s| parse_wind_segment(s, cli.units))
                .collect::<Result<Vec<_>, String>>()?;

            // Construct TrajectoryConfig once, used by all code paths
            let traj_config = TrajectoryConfig {
                velocity: velocity_metric,
                angle: muzzle_angle,
                bc: trued_bc,
                mass: mass_metric,
                diameter: diameter_metric,
                drag_model,
                max_range: max_range_metric,
                time_step,
                temperature: temperature_metric,
                pressure: pressure_metric,
                humidity: final_humidity,
                altitude: altitude_metric,
                wind_speed: wind_speed_metric,
                wind_direction: final_wind_direction,
                wind_segments,
                output,
                full,
                units: cli.units,
                sight_height: sight_height_metric,
                bore_height: bore_height_metric,
                ignore_ground_impact,
                use_bc_segments,
                use_cluster_bc,
                bc_table_segments: bc_table_segments.clone(),
                enable_magnus,
                enable_coriolis,
                enable_spin_drift,
                enable_aerodynamic_jump,
                enable_wind_shear,
                enable_pitch_damping,
                enable_precession,
                sample_trajectory,
                sample_interval,
                use_rk4: !use_euler,
                use_rk45: !use_rk4_fixed,
                twist_rate,
                twist_right,
                latitude,
                shooting_angle,
                use_powder_sensitivity,
                powder_temp_sensitivity,
                powder_temp,
                pdf_metadata: pdf_metadata.clone(),
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
                    let zero_range_metric = final_auto_zero.map(|d| UnitConverter::distance_to_metric(d, cli.units));

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
                        .sample_interval(sample_interval)
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
                    // twist_rate is inches-per-turn for ALL unit systems (documented at the
                    // --twist-rate flag, and used as-is by the local solver / TrajectoryConfig and
                    // the compare-mode local inputs), and api_client sends it verbatim — so forward
                    // the RAW value. The original code sent meters (~33x too small); converting
                    // metric by /25.4 here would instead make the API disagree with local under
                    // --compare. No conversion: all paths use the documented inches/turn.
                    if let Some(twist) = twist_rate {
                        request.twist_rate = Some(twist);
                    }

                    let api_client = ApiClient::new(&api_url, api_timeout);

                    let api_result = api_client.calculate_trajectory(&request);

                    match (&api_result, compare) {
                        (Ok(api_response), true) => {
                            // Compare mode: run local calculation first, then display side-by-side (MBA-720)
                            eprintln!("Running comparison between local and API calculations...\n");

                            // Run local trajectory to get comparison data BEFORE displaying the table
                            let local_result = {
                                let drag_model_enum = match drag_model {
                                    DragModelArg::G1 => DragModel::G1,
                                    DragModelArg::G7 => DragModel::G7,
                                };

                                let local_inputs = BallisticInputs {
                                    bc_value: trued_bc,
                                    bc_type: drag_model_enum,
                                    bullet_mass: mass_metric,
                                    muzzle_velocity: velocity_metric,
                                    bullet_diameter: diameter_metric,
                                    bullet_length: diameter_metric * 4.5,
                                    muzzle_angle: muzzle_angle.to_radians(),
                                    target_distance: max_range_metric,
                                    azimuth_angle: 0.0,
                                    shooting_angle: shooting_angle.to_radians(),
                                    sight_height: sight_height_metric,
                                    muzzle_height: bore_height_metric,
                                    target_height: 0.0,
                                    ground_threshold: if ignore_ground_impact { f64::NEG_INFINITY } else { -bore_height_metric },
                                    altitude: altitude_metric,
                                    temperature: temperature_metric,
                                    pressure: pressure_metric,
                                    humidity: final_humidity,
                                    latitude,
                                    wind_speed: wind_speed_metric,
                                    wind_angle: final_wind_direction,
                                    twist_rate: twist_rate.unwrap_or(12.0),
                                    is_twist_right: twist_right,
                                    caliber_inches: diameter_metric / 0.0254,
                                    weight_grains: mass_metric / 0.00006479891,
                                    manufacturer: None,
                                    bullet_model: None,
                                    bullet_id: None,
                                    bullet_cluster: None,
                                    use_rk4: !use_euler,
                                    use_adaptive_rk45: !use_rk4_fixed,
                                    enable_advanced_effects: enable_magnus || enable_coriolis,
                                    enable_magnus,
                                    enable_coriolis,
                                    use_powder_sensitivity,
                                    powder_temp_sensitivity: if use_powder_sensitivity {
                                        UnitConverter::velocity_to_metric(powder_temp_sensitivity, cli.units)
                                            / UnitConverter::temperature_to_metric(1.0, cli.units)
                                    } else { 0.0 },
                                    powder_temp: UnitConverter::temperature_to_metric(powder_temp, cli.units),
                                    tipoff_yaw: 0.0,
                                    tipoff_decay_distance: 50.0,
                                    use_bc_segments: use_bc_segments || bc_table_segments.is_some(),
                                    bc_segments: None,
                                    bc_segments_data: bc_table_segments.clone().or_else(|| {
                                        if use_bc_segments {
                                            use ballistics_engine::bc_estimation::BCSegmentEstimator;
                                            Some(BCSegmentEstimator::estimate_bc_segments(
                                                trued_bc,
                                                diameter_metric / 0.0254,
                                                mass_metric / 0.00006479891,
                                                "",
                                                match drag_model { DragModelArg::G1 => "G1", DragModelArg::G7 => "G7" },
                                            ))
                                        } else { None }
                                    }),
                                    use_enhanced_spin_drift: enable_spin_drift,
                                    enable_aerodynamic_jump,
                                    use_form_factor: false,
                                    enable_wind_shear,
                                    wind_shear_model: if enable_wind_shear { "exponential".to_string() } else { "none".to_string() },
                                    enable_trajectory_sampling: sample_trajectory,
                                    sample_interval,
                                    enable_pitch_damping,
                                    enable_precession_nutation: enable_precession,
                                    use_cluster_bc,
                                    bc_type_str: None,
                                    custom_drag_table: None,
                                };

                                let local_wind = WindConditions {
                                    speed: wind_speed_metric,
                                    direction: final_wind_direction.to_radians(),
                                    ..Default::default()
                                };
                                let local_atmo = AtmosphericConditions {
                                    temperature: temperature_metric,
                                    pressure: pressure_metric,
                                    humidity: final_humidity,
                                    altitude: altitude_metric,
                                    ..Default::default()
                                };

                                let mut local_solver = TrajectorySolver::new(local_inputs, local_wind, local_atmo);
                                local_solver.set_max_range(max_range_metric);
                                local_solver.set_time_step(time_step);
                                local_solver.solve().ok()
                            };

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

                            let (dist_unit, _drop_unit, _vel_unit) = match cli.units {
                                UnitSystem::Metric => ("m", "m", "m/s"),
                                UnitSystem::Imperial => ("yd", "in", "fps"),
                            };

                            println!("║ Range {} │ API Drop │Local Drop│  Δ Drop  ║", dist_unit);
                            println!("╠══════════╪══════════╪══════════╪══════════╣");

                            for point in api_response.trajectory.iter().take(10) {
                                let range_display = UnitConverter::distance_from_metric(point.range, cli.units);
                                let api_drop_display = if cli.units == UnitSystem::Imperial {
                                    point.drop * 39.3701
                                } else {
                                    point.drop
                                };

                                if let Some(ref local_res) = local_result {
                                    if let Some(pos) = local_res.position_at_range(point.range) {
                                        let local_drop_display = if cli.units == UnitSystem::Imperial {
                                            pos.y * 39.3701
                                        } else {
                                            pos.y
                                        };
                                        let delta = api_drop_display - local_drop_display;
                                        println!(
                                            "║ {:>8.1} │ {:>8.2} │ {:>8.2} │ {:>+8.2} ║",
                                            range_display, api_drop_display, local_drop_display, delta
                                        );
                                    } else {
                                        println!(
                                            "║ {:>8.1} │ {:>8.2} │    ---   │    ---   ║",
                                            range_display, api_drop_display
                                        );
                                    }
                                } else {
                                    println!(
                                        "║ {:>8.1} │ {:>8.2} │  error   │    ---   ║",
                                        range_display, api_drop_display
                                    );
                                }
                            }
                            println!("╚══════════╧══════════╧══════════╧══════════╝");
                            println!();
                        }
                        (Ok(api_response), false) => {
                            // Online mode only: display API results
                            display_api_trajectory_result(api_response, output, cli.units, full)?;
                        }
                        (Err(e), _) => {
                            if offline_fallback {
                                eprintln!("Warning: API request failed: {}", e);
                                eprintln!("Falling back to local calculation...\n");

                                run_trajectory(&traj_config)?;
                            } else {
                                eprintln!("Error: API request failed: {}", e);
                                eprintln!("Hint: Use --offline-fallback to use local calculation on API failure");
                                std::process::exit(1);
                            }
                        }
                    }
                } else {
                    // Local calculation (default)
                    run_trajectory(&traj_config)?;
                }
            }

            // Non-online feature: just run local calculation
            #[cfg(not(feature = "online"))]
            {
                run_trajectory(&traj_config)?;
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
            // Convert inputs to metric (MBA-716)
            let velocity_metric = UnitConverter::velocity_to_metric(velocity, cli.units);
            let mass_metric = UnitConverter::mass_to_metric(bullet_mass, cli.units);
            let diameter_metric = UnitConverter::diameter_to_metric(bullet_diameter, cli.units);
            let velocity_std_metric = UnitConverter::velocity_to_metric(velocity_std, cli.units);
            let wind_std_metric = UnitConverter::wind_to_metric(wind_std, cli.units);
            let wind_speed_metric = UnitConverter::wind_to_metric(wind_speed, cli.units);
            let target_distance_metric = target_distance.map(|d| UnitConverter::distance_to_metric(d, cli.units));
            run_monte_carlo(
                velocity_metric,
                angle,
                bc,
                mass_metric,
                diameter_metric,
                num_sims,
                velocity_std_metric,
                angle_std,
                bc_std,
                wind_std_metric,
                wind_speed_metric,
                wind_direction,
                target_distance_metric,
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
            temperature,
            pressure,
            humidity,
            altitude,
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
            // Convert atmospheric conditions to metric
            let temperature_metric = UnitConverter::temperature_to_metric(temperature, cli.units);
            let pressure_metric = UnitConverter::pressure_to_metric(pressure, cli.units);
            let altitude_metric = UnitConverter::altitude_to_metric(altitude, cli.units);

            run_zero_calculation(
                velocity_metric,
                bc,
                mass_metric,
                diameter_metric,
                target_distance_metric,
                target_height_metric,
                sight_height_metric,
                temperature_metric,
                pressure_metric,
                humidity,
                altitude_metric,
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
            // Convert inputs to metric (MBA-716)
            let velocity_metric = UnitConverter::velocity_to_metric(velocity, cli.units);
            let mass_metric = UnitConverter::mass_to_metric(bullet_mass, cli.units);
            let diameter_metric = UnitConverter::diameter_to_metric(bullet_diameter, cli.units);
            let distance1_metric = UnitConverter::distance_to_metric(distance1, cli.units);
            let drop1_metric = UnitConverter::distance_to_metric(drop1, cli.units);
            let distance2_metric = UnitConverter::distance_to_metric(distance2, cli.units);
            let drop2_metric = UnitConverter::distance_to_metric(drop2, cli.units);
            run_bc_estimation(
                velocity_metric,
                mass_metric,
                diameter_metric,
                distance1_metric,
                drop1_metric,
                distance2_metric,
                drop2_metric,
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

        Commands::TrueVelocity {
            measured_drop,
            range,
            bc,
            drag_model,
            mass,
            diameter,
            chrono_velocity,
            zero_distance,
            sight_height,
            temperature,
            pressure,
            humidity,
            altitude,
            units,
            bc_table_dir,
            #[cfg(feature = "online")]
            bc_table_auto,
            #[cfg(feature = "online")]
            bc_table_url,
            offline,
            #[cfg(feature = "online")]
            offline_fallback,
            #[cfg(feature = "online")]
            api_url,
            #[cfg(feature = "online")]
            api_timeout,
            bullet_length,
            output,
        } => {
            // Convert to imperial for calculations (internal calculations use imperial)
            let range_yd = match units {
                UnitSystem::Imperial => range,
                UnitSystem::Metric => range / 0.9144, // meters to yards
            };
            let weight_gr = match units {
                UnitSystem::Imperial => mass,
                UnitSystem::Metric => mass / 0.0647989, // grams to grains
            };
            let caliber_in = match units {
                UnitSystem::Imperial => diameter,
                UnitSystem::Metric => diameter / 25.4, // mm to inches
            };
            let chrono_fps = chrono_velocity.map(|v| match units {
                UnitSystem::Imperial => v,
                UnitSystem::Metric => v / 0.3048, // m/s to fps
            });
            let zero_yd = match units {
                UnitSystem::Imperial => zero_distance,
                UnitSystem::Metric => zero_distance / 0.9144,
            };
            let sight_in = match units {
                UnitSystem::Imperial => sight_height,
                UnitSystem::Metric => sight_height / 25.4, // mm to inches
            };
            let temp_f = match units {
                UnitSystem::Imperial => temperature,
                UnitSystem::Metric => temperature * 9.0 / 5.0 + 32.0, // C to F
            };
            let press_inhg = match units {
                UnitSystem::Imperial => pressure,
                UnitSystem::Metric => pressure / 33.8639, // hPa to inHg
            };
            let alt_ft = match units {
                UnitSystem::Imperial => altitude,
                UnitSystem::Metric => altitude / 0.3048, // meters to feet
            };

            let drag_str = match drag_model {
                DragModelArg::G1 => "G1",
                DragModelArg::G7 => "G7",
            };

            // Load BC5D tables if available
            // Determine effective table directory: explicit --bc-table-dir, or auto-download cache
            #[cfg(feature = "online")]
            let effective_bc_table_dir: Option<PathBuf> = if bc_table_dir.is_some() {
                bc_table_dir.clone()
            } else if bc_table_auto {
                // Auto-download mode: ensure table is available, use cache directory
                match Bc5dDownloader::new(&bc_table_url, false) {
                    Ok(mut downloader) => {
                        match downloader.ensure_table(caliber_in) {
                            Ok(table_path) => {
                                // Return the parent directory (cache dir) as table directory
                                table_path.parent().map(|p| p.to_path_buf())
                            }
                            Err(e) => {
                                eprintln!("Warning: {}", e);
                                eprintln!("         Continuing without BC5D correction table.");
                                None
                            }
                        }
                    }
                    Err(e) => {
                        eprintln!("Warning: BC5D table download failed: {}", e);
                        eprintln!("         Continuing without BC5D correction table.");
                        None
                    }
                }
            } else {
                None
            };

            #[cfg(not(feature = "online"))]
            let effective_bc_table_dir: Option<PathBuf> = bc_table_dir.clone();

            // Load BC segments from BC5D table if available (MBA-744: uses shared helper)
            let bc_segments: Option<Vec<BCSegmentData>> = if let Some(table_dir) = &effective_bc_table_dir {
                let mut manager = Bc5dTableManager::new(table_dir);
                generate_bc5d_segments(
                    &mut manager,
                    caliber_in,
                    bc,
                    drag_str,
                    weight_gr,
                    None,
                    bullet_length,
                )
            } else {
                None
            };

            // Determine if we should use local calculation
            // Use local if: --offline flag, OR no online feature, OR (no online flag and has bc_table_dir)
            #[cfg(feature = "online")]
            let use_local = offline;
            #[cfg(not(feature = "online"))]
            let use_local = true; // Always use local when online feature not available

            if use_local {
                // Run local calculation
                eprintln!("Calculating effective velocity locally...");

                match calculate_true_velocity_local(
                    measured_drop,
                    range_yd,
                    bc,
                    drag_model,
                    weight_gr,
                    caliber_in,
                    zero_yd,
                    sight_in,
                    temp_f,
                    press_inhg,
                    humidity,
                    alt_ft,
                    &bc_segments,
                ) {
                    Ok(result) => {
                        // Calculate velocity adjustment if chrono provided
                        let velocity_adjustment = chrono_fps.map(|c| result.effective_velocity_fps - c);
                        let adjustment_percent = chrono_fps.map(|c| (result.effective_velocity_fps - c) / c * 100.0);

                        // Convert output velocity back to user's units if needed
                        let effective_vel = match units {
                            UnitSystem::Imperial => result.effective_velocity_fps,
                            UnitSystem::Metric => result.effective_velocity_fps * 0.3048,
                        };
                        let vel_unit = match units {
                            UnitSystem::Imperial => "fps",
                            UnitSystem::Metric => "m/s",
                        };

                        display_true_velocity_result(
                            effective_vel,
                            vel_unit,
                            velocity_adjustment,
                            adjustment_percent,
                            &result.confidence,
                            result.iterations,
                            result.final_error_mil,
                            result.calculated_drop_mil,
                            measured_drop,
                            units,
                            output,
                            bc_segments.is_some(),
                        )?;
                    }
                    Err(e) => {
                        eprintln!("Error: Local calculation failed: {}", e);
                        std::process::exit(1);
                    }
                }
            } else {
                // Online mode
                #[cfg(feature = "online")]
                {
                    // Check TOS acceptance first
                    if !check_tos_accepted() {
                        if !prompt_tos_acceptance()? {
                            eprintln!("Cannot use online features without accepting Terms of Service.");
                            std::process::exit(1);
                        }
                    }

                    let request = TrueVelocityRequest {
                        measured_drop_mil: measured_drop,
                        range_yd,
                        bc,
                        drag_model: drag_str.to_string(),
                        weight_gr,
                        caliber: caliber_in,
                        zero_range_yd: Some(zero_yd),
                        chrono_velocity_fps: chrono_fps,
                        altitude_ft: Some(alt_ft),
                        temperature_f: Some(temp_f),
                        pressure_inhg: Some(press_inhg),
                        humidity: Some(humidity),
                        sight_height_in: Some(sight_in),
                        use_bc_enhancement: Some(true),
                    };

                    let api_client = ApiClient::new(&api_url, api_timeout);

                    eprintln!("Calculating effective velocity via API...");

                    match api_client.true_velocity(&request) {
                        Ok(response) => {
                            // Convert output velocity back to user's units if needed
                            let effective_vel = match units {
                                UnitSystem::Imperial => response.effective_velocity_fps,
                                UnitSystem::Metric => response.effective_velocity_fps * 0.3048,
                            };
                            let vel_unit = match units {
                                UnitSystem::Imperial => "fps",
                                UnitSystem::Metric => "m/s",
                            };

                            display_true_velocity_result(
                                effective_vel,
                                vel_unit,
                                response.velocity_adjustment_fps,
                                response.adjustment_percent,
                                &response.confidence,
                                response.iterations,
                                response.final_error_mil,
                                response.calculated_drop_mil,
                                measured_drop,
                                units,
                                output,
                                false,
                            )?;
                        }
                        Err(e) => {
                            if offline_fallback {
                                eprintln!("Warning: API request failed: {}", e);
                                eprintln!("Falling back to local calculation...\n");

                                match calculate_true_velocity_local(
                                    measured_drop,
                                    range_yd,
                                    bc,
                                    drag_model,
                                    weight_gr,
                                    caliber_in,
                                    zero_yd,
                                    sight_in,
                                    temp_f,
                                    press_inhg,
                                    humidity,
                                    alt_ft,
                                    &bc_segments,
                                ) {
                                    Ok(result) => {
                                        let velocity_adjustment = chrono_fps.map(|c| result.effective_velocity_fps - c);
                                        let adjustment_percent = chrono_fps.map(|c| (result.effective_velocity_fps - c) / c * 100.0);

                                        let effective_vel = match units {
                                            UnitSystem::Imperial => result.effective_velocity_fps,
                                            UnitSystem::Metric => result.effective_velocity_fps * 0.3048,
                                        };
                                        let vel_unit = match units {
                                            UnitSystem::Imperial => "fps",
                                            UnitSystem::Metric => "m/s",
                                        };

                                        display_true_velocity_result(
                                            effective_vel,
                                            vel_unit,
                                            velocity_adjustment,
                                            adjustment_percent,
                                            &result.confidence,
                                            result.iterations,
                                            result.final_error_mil,
                                            result.calculated_drop_mil,
                                            measured_drop,
                                            units,
                                            output,
                                            bc_segments.is_some(),
                                        )?;
                                    }
                                    Err(fallback_err) => {
                                        eprintln!("Error: Fallback calculation also failed: {}", fallback_err);
                                        std::process::exit(1);
                                    }
                                }
                            } else {
                                eprintln!("API Error: {}", e);
                                eprintln!("Hint: Use --offline for local calculation or --offline-fallback for automatic fallback");
                                std::process::exit(1);
                            }
                        }
                    }
                }

                #[cfg(not(feature = "online"))]
                {
                    // This branch should never be reached since use_local is always true
                    // when online feature is not available, but add for completeness
                    eprintln!("Error: Online mode not available (compile with --features online)");
                    std::process::exit(1);
                }
            }
        }

        Commands::Mpbr {
            profile,
            velocity,
            bc,
            mass,
            diameter,
            drag_model,
            vital_zone,
            sight_height,
            temperature,
            pressure,
            humidity,
            altitude,
            output,
        } => {
            // Load profile if specified
            let profile_data = profile.as_ref().map(|name| {
                load_profile(name).unwrap_or_else(|e| {
                    eprintln!("Error: {}", e);
                    std::process::exit(1);
                })
            });

            let final_velocity = resolve_param(velocity, &profile_data, |p| p.velocity)
                .unwrap_or_else(|| { eprintln!("Error: --velocity is required (or use --profile)"); std::process::exit(1); });
            let final_bc = resolve_param(bc, &profile_data, |p| p.bc)
                .unwrap_or_else(|| { eprintln!("Error: --bc is required (or use --profile)"); std::process::exit(1); });
            let final_mass = resolve_param(mass, &profile_data, |p| p.mass)
                .unwrap_or_else(|| { eprintln!("Error: --mass is required (or use --profile)"); std::process::exit(1); });
            let final_diameter = resolve_param(diameter, &profile_data, |p| p.diameter)
                .unwrap_or_else(|| { eprintln!("Error: --diameter is required (or use --profile)"); std::process::exit(1); });
            let final_sight_height = sight_height
                .or_else(|| profile_data.as_ref().and_then(|p| p.sight_height))
                .unwrap_or(match cli.units { UnitSystem::Imperial => 2.0, UnitSystem::Metric => 50.0 });
            let final_drag_model = if profile_data.is_some() && velocity.is_none() {
                parse_drag_model_arg(&profile_data.as_ref().unwrap().drag_model)
            } else {
                drag_model
            };

            handle_mpbr(
                final_velocity, final_bc, final_mass, final_diameter,
                final_drag_model, vital_zone, final_sight_height,
                temperature, pressure, humidity, altitude,
                cli.units, output,
            )?;
        }

        Commands::ComeUps {
            profile,
            velocity,
            bc,
            mass,
            diameter,
            drag_model,
            zero_distance,
            start,
            end,
            step,
            adjustment_unit,
            sight_height,
            temperature,
            pressure,
            humidity,
            altitude,
            wind_speed,
            wind_direction,
            output,
        } => {
            let profile_data = profile.as_ref().map(|name| {
                load_profile(name).unwrap_or_else(|e| {
                    eprintln!("Error: {}", e);
                    std::process::exit(1);
                })
            });

            let final_velocity = resolve_param(velocity, &profile_data, |p| p.velocity)
                .unwrap_or_else(|| { eprintln!("Error: --velocity is required (or use --profile)"); std::process::exit(1); });
            let final_bc = resolve_param(bc, &profile_data, |p| p.bc)
                .unwrap_or_else(|| { eprintln!("Error: --bc is required (or use --profile)"); std::process::exit(1); });
            let final_mass = resolve_param(mass, &profile_data, |p| p.mass)
                .unwrap_or_else(|| { eprintln!("Error: --mass is required (or use --profile)"); std::process::exit(1); });
            let final_diameter = resolve_param(diameter, &profile_data, |p| p.diameter)
                .unwrap_or_else(|| { eprintln!("Error: --diameter is required (or use --profile)"); std::process::exit(1); });
            let final_sight_height = sight_height
                .or_else(|| profile_data.as_ref().and_then(|p| p.sight_height))
                .unwrap_or(match cli.units { UnitSystem::Imperial => 2.0, UnitSystem::Metric => 50.0 });
            let final_drag_model = if profile_data.is_some() && velocity.is_none() {
                parse_drag_model_arg(&profile_data.as_ref().unwrap().drag_model)
            } else {
                drag_model
            };

            handle_come_ups(
                final_velocity, final_bc, final_mass, final_diameter,
                final_drag_model, zero_distance, start, end, step,
                adjustment_unit, final_sight_height,
                temperature, pressure, humidity, altitude,
                wind_speed, wind_direction,
                cli.units, output,
            )?;
        }

        Commands::WindCard {
            profile,
            velocity,
            bc,
            mass,
            diameter,
            drag_model,
            zero_distance,
            wind_speeds,
            start,
            end,
            step,
            adjustment_unit,
            sight_height,
            temperature,
            pressure,
            humidity,
            altitude,
            output,
        } => {
            let profile_data = profile.as_ref().map(|name| {
                load_profile(name).unwrap_or_else(|e| {
                    eprintln!("Error: {}", e);
                    std::process::exit(1);
                })
            });

            let final_velocity = resolve_param(velocity, &profile_data, |p| p.velocity)
                .unwrap_or_else(|| { eprintln!("Error: --velocity is required (or use --profile)"); std::process::exit(1); });
            let final_bc = resolve_param(bc, &profile_data, |p| p.bc)
                .unwrap_or_else(|| { eprintln!("Error: --bc is required (or use --profile)"); std::process::exit(1); });
            let final_mass = resolve_param(mass, &profile_data, |p| p.mass)
                .unwrap_or_else(|| { eprintln!("Error: --mass is required (or use --profile)"); std::process::exit(1); });
            let final_diameter = resolve_param(diameter, &profile_data, |p| p.diameter)
                .unwrap_or_else(|| { eprintln!("Error: --diameter is required (or use --profile)"); std::process::exit(1); });
            let final_sight_height = sight_height
                .or_else(|| profile_data.as_ref().and_then(|p| p.sight_height))
                .unwrap_or(match cli.units { UnitSystem::Imperial => 2.0, UnitSystem::Metric => 50.0 });
            let final_drag_model = if profile_data.is_some() && velocity.is_none() {
                parse_drag_model_arg(&profile_data.as_ref().unwrap().drag_model)
            } else {
                drag_model
            };

            // Parse wind speeds
            let ws_vec: Vec<f64> = wind_speeds.split(',')
                .filter_map(|s| s.trim().parse::<f64>().ok())
                .collect();
            if ws_vec.is_empty() {
                eprintln!("Error: --wind-speeds must contain at least one valid number (e.g., '5,10,15,20')");
                std::process::exit(1);
            }

            handle_wind_card(
                final_velocity, final_bc, final_mass, final_diameter,
                final_drag_model, zero_distance, &ws_vec,
                start, end, step,
                adjustment_unit, final_sight_height,
                temperature, pressure, humidity, altitude,
                cli.units, output,
            )?;
        }


        Commands::Stability {
            profile,
            mass,
            diameter,
            length,
            twist_rate,
            velocity,
            temperature,
            pressure,
            altitude,
            output,
        } => {
            let profile_data = profile.as_ref().map(|name| {
                load_profile(name).unwrap_or_else(|e| {
                    eprintln!("Error: {}", e);
                    std::process::exit(1);
                })
            });

            let final_mass = resolve_param(mass, &profile_data, |p| p.mass)
                .unwrap_or_else(|| { eprintln!("Error: --mass is required (or use --profile)"); std::process::exit(1); });
            let final_diameter = resolve_param(diameter, &profile_data, |p| p.diameter)
                .unwrap_or_else(|| { eprintln!("Error: --diameter is required (or use --profile)"); std::process::exit(1); });
            let final_velocity = velocity.unwrap_or(match cli.units {
                UnitSystem::Imperial => 2700.0,
                UnitSystem::Metric => 823.0,
            });

            handle_stability(
                final_mass, final_diameter, length, twist_rate, final_velocity,
                temperature, pressure, altitude,
                cli.units, output,
            )?;
        }

        Commands::RangeTable {
            profile,
            velocity,
            bc,
            mass,
            diameter,
            drag_model,
            zero_distance,
            start,
            end,
            step,
            wind_speed,
            wind_direction,
            adjustment_unit,
            sight_height,
            temperature,
            pressure,
            humidity,
            altitude,
            output,
        } => {
            let profile_data = profile.as_ref().map(|name| {
                load_profile(name).unwrap_or_else(|e| {
                    eprintln!("Error: {}", e);
                    std::process::exit(1);
                })
            });

            let final_velocity = resolve_param(velocity, &profile_data, |p| p.velocity)
                .unwrap_or_else(|| { eprintln!("Error: --velocity is required (or use --profile)"); std::process::exit(1); });
            let final_bc = resolve_param(bc, &profile_data, |p| p.bc)
                .unwrap_or_else(|| { eprintln!("Error: --bc is required (or use --profile)"); std::process::exit(1); });
            let final_mass = resolve_param(mass, &profile_data, |p| p.mass)
                .unwrap_or_else(|| { eprintln!("Error: --mass is required (or use --profile)"); std::process::exit(1); });
            let final_diameter = resolve_param(diameter, &profile_data, |p| p.diameter)
                .unwrap_or_else(|| { eprintln!("Error: --diameter is required (or use --profile)"); std::process::exit(1); });
            let final_sight_height = sight_height
                .or_else(|| profile_data.as_ref().and_then(|p| p.sight_height))
                .unwrap_or(match cli.units { UnitSystem::Imperial => 2.0, UnitSystem::Metric => 50.0 });
            let final_drag_model = if profile_data.is_some() && velocity.is_none() {
                parse_drag_model_arg(&profile_data.as_ref().unwrap().drag_model)
            } else {
                drag_model
            };

            handle_range_table(
                final_velocity, final_bc, final_mass, final_diameter,
                final_drag_model, zero_distance, start, end, step,
                wind_speed, wind_direction,
                adjustment_unit, final_sight_height,
                temperature, pressure, humidity, altitude,
                cli.units, output,
            )?;
        }

        Commands::Profile { action } => {
            match action {
                ProfileAction::Save {
                    name, velocity, bc, mass, diameter, drag_model,
                    twist_rate, sight_height, zero_distance,
                    temperature, pressure, humidity, altitude,
                    bullet_name,
                    wind_speed, wind_direction, shooting_angle,
                    auto_zero, twist_right, use_bc_segments, bullet_length,
                } => {
                    let drag_str = match drag_model {
                        DragModelArg::G1 => "G1",
                        DragModelArg::G7 => "G7",
                    };
                    let units_str = match cli.units {
                        UnitSystem::Imperial => "imperial",
                        UnitSystem::Metric => "metric",
                    };

                    let profile = ProfileData {
                        name: name.clone(),
                        velocity,
                        bc,
                        mass,
                        diameter,
                        drag_model: drag_str.to_string(),
                        twist_rate,
                        sight_height,
                        zero_distance,
                        units: units_str.to_string(),
                        temperature,
                        pressure,
                        humidity,
                        altitude,
                        bullet_name,
                        created: Some(timestamp_string()),
                        wind_speed, wind_direction, shooting_angle,
                        auto_zero, twist_right, use_bc_segments, bullet_length,
                    };

                    let path = save_profile(&profile)?;
                    eprintln!("Profile '{}' saved to {:?}", name, path);
                }

                ProfileAction::List => {
                    let profiles = list_profiles()?;
                    if profiles.is_empty() {
                        println!("No saved profiles. Use 'ballistics profile save <name> ...' to create one.");
                    } else {
                        println!("Saved Profiles:");
                        println!("┌────────────────────┬────────┬───────┬────────┬──────────┬──────────┐");
                        println!("│ Name               │ Vel    │ BC    │ Mass   │ Diameter │ Drag     │");
                        println!("├────────────────────┼────────┼───────┼────────┼──────────┼──────────┤");
                        for p in &profiles {
                            println!("│ {:<18} │{:>7.0} │{:>6.3} │{:>7.1} │{:>9.3} │ {:<8} │",
                                     p.name, p.velocity, p.bc, p.mass, p.diameter, p.drag_model);
                        }
                        println!("└────────────────────┴────────┴───────┴────────┴──────────┴──────────┘");
                    }
                }

                ProfileAction::Show { name } => {
                    let profile = load_profile(&name)?;
                    println!();
                    println!("Profile: {}", profile.name);
                    println!("╔════════════════════════════════════════╗");
                    println!("║  Velocity:      {:>10.1} {:<10}  ║",
                        profile.velocity,
                        if profile.units == "metric" { "m/s" } else { "fps" });
                    println!("║  BC:            {:>10.4}             ║", profile.bc);
                    println!("║  Mass:          {:>10.1} {:<10}  ║",
                        profile.mass,
                        if profile.units == "metric" { "g" } else { "gr" });
                    println!("║  Diameter:      {:>10.3} {:<10}  ║",
                        profile.diameter,
                        if profile.units == "metric" { "mm" } else { "in" });
                    println!("║  Drag model:    {:>10}             ║", profile.drag_model);
                    if let Some(tw) = profile.twist_rate {
                        println!("║  Twist rate:    {:>10.1}             ║", tw);
                    }
                    if let Some(sh) = profile.sight_height {
                        println!("║  Sight height:  {:>10.2} {:<10}  ║",
                            sh,
                            if profile.units == "metric" { "mm" } else { "in" });
                    }
                    if let Some(zd) = profile.zero_distance {
                        println!("║  Zero distance: {:>10.0} {:<10}  ║",
                            zd,
                            if profile.units == "metric" { "m" } else { "yd" });
                    }
                    if let Some(ref bn) = profile.bullet_name {
                        println!("║  Bullet:        {:<24}  ║", bn);
                    }
                    if let Some(ws) = profile.wind_speed {
                        println!("║  Wind speed:    {:>10.1} {:<10}  ║",
                            ws,
                            if profile.units == "metric" { "m/s" } else { "mph" });
                    }
                    if let Some(wd) = profile.wind_direction {
                        println!("║  Wind dir:      {:>10.1}°            ║", wd);
                    }
                    if let Some(sa) = profile.shooting_angle {
                        println!("║  Shoot angle:   {:>10.1}°            ║", sa);
                    }
                    if let Some(az) = profile.auto_zero {
                        println!("║  Auto-zero:     {:>10.0} {:<10}  ║",
                            az,
                            if profile.units == "metric" { "m" } else { "yd" });
                    }
                    if let Some(tr) = profile.twist_right {
                        println!("║  Twist:         {:>10}             ║",
                            if tr { "Right" } else { "Left" });
                    }
                    if let Some(ubc) = profile.use_bc_segments {
                        println!("║  BC segments:   {:>10}             ║",
                            if ubc { "Enabled" } else { "Disabled" });
                    }
                    if let Some(bl) = profile.bullet_length {
                        println!("║  Bullet length: {:>10.3} {:<10}  ║",
                            bl,
                            if profile.units == "metric" { "mm" } else { "in" });
                    }
                    println!("║  Temperature:   {:>10.1} {:<10}  ║",
                        profile.temperature,
                        if profile.units == "metric" { "°C" } else { "°F" });
                    println!("║  Pressure:      {:>10.2} {:<10}  ║",
                        profile.pressure,
                        if profile.units == "metric" { "hPa" } else { "inHg" });
                    println!("║  Humidity:      {:>10.1}%            ║", profile.humidity);
                    println!("║  Altitude:      {:>10.0} {:<10}  ║",
                        profile.altitude,
                        if profile.units == "metric" { "m" } else { "ft" });
                    println!("╚════════════════════════════════════════╝");
                }

                ProfileAction::Delete { name } => {
                    delete_profile(&name)?;
                    eprintln!("Profile '{}' deleted.", name);
                }
            }
        }

        Commands::Completions { shell } => {
            let mut cmd = Cli::command();
            clap_complete::generate(shell, &mut cmd, "ballistics", &mut std::io::stdout());
        }
    }

    Ok(())
}

/// PDF-specific metadata for dope card generation
#[derive(Debug, Clone, Default)]
struct PdfMetadata {
    rifle_name: String,
    location_name: String,
    powder: String,
    bullet_name: String,
    target_speed_mph: f64,
    output_file: Option<PathBuf>,
    // Original imperial values for header display
    velocity_fps: f64,
    temperature_f: f64,
    pressure_inhg: f64,
    altitude_ft: f64,
    wind_speed_mph: f64,
    weight_gr: f64,
    font_scale: f32,
    bold_data: bool,
}


// ============================================================================
// BC5D Segment Generation Helper (MBA-744)
// ============================================================================

/// Generate velocity-dependent BC segments from a Bc5dTableManager.
fn generate_bc5d_segments(
    manager: &mut Bc5dTableManager,
    caliber: f64,
    base_bc: f64,
    drag_type: &str,
    weight_grains: f64,
    muzzle_velocity_fps: Option<f64>,
    bullet_length_in: Option<f64>,
) -> Option<Vec<BCSegmentData>> {
    let mut velocity_breakpoints: Vec<f64> = vec![
        4000.0, 3500.0, 3000.0, 2700.0, 2500.0, 2300.0, 2100.0, 2000.0, 1900.0, 1800.0,
        1700.0, 1600.0, 1500.0, 1400.0, 1350.0, 1300.0, 1250.0,
        1200.0, 1150.0, 1100.0, 1050.0, 1000.0, 950.0, 900.0,
        850.0, 800.0, 700.0, 600.0, 500.0,
    ];

    if let Some(mv) = muzzle_velocity_fps {
        velocity_breakpoints.push(mv);
    }

    let mut valid_velocities: Vec<f64> = velocity_breakpoints
        .into_iter()
        .filter(|&v| {
            v >= 500.0 && muzzle_velocity_fps.map_or(true, |mv| v <= mv)
        })
        .collect();
    valid_velocities.sort_by(|a, b| b.partial_cmp(a).unwrap());
    valid_velocities.dedup();

    let reference_muzzle_velocity = valid_velocities.first().copied().unwrap_or(3000.0);

    let mut segments: Vec<BCSegmentData> = Vec::new();
    let mut any_correction_found = false;

    for i in 0..valid_velocities.len().saturating_sub(1) {
        let vel_max = valid_velocities[i];
        let vel_min = valid_velocities[i + 1];
        let vel_mid = (vel_max + vel_min) / 2.0;

        if let Ok(correction) = manager.lookup(caliber, weight_grains, base_bc, reference_muzzle_velocity, vel_mid, drag_type) {
            any_correction_found = true;
            let segment_bc = base_bc * correction;

            segments.push(BCSegmentData {
                velocity_min: vel_min,
                velocity_max: vel_max,
                bc_value: segment_bc,
            });
        }
    }

    if any_correction_found && !segments.is_empty() {
        let length_display = bullet_length_in.unwrap_or_else(|| caliber * 3.5);
        eprintln!("BC5D Table: Generated {} velocity-dependent BC segments", segments.len());
        eprintln!("BC5D Table: {:.3}\" caliber, {:.1}gr, {:.3}\" length (est)", caliber, weight_grains, length_display);
        let min_bc = segments.iter().map(|s| s.bc_value).fold(f64::INFINITY, f64::min);
        let max_bc = segments.iter().map(|s| s.bc_value).fold(f64::NEG_INFINITY, f64::max);
        eprintln!("BC5D Table: BC range {:.5} - {:.5} across velocity envelope", min_bc, max_bc);
        Some(segments)
    } else {
        eprintln!("Warning: No BC5D table data found for caliber {:.3}\"", caliber);
        None
    }
}

fn run_trajectory(config: &TrajectoryConfig) -> Result<(), Box<dyn Error>> {
    // Destructure config for convenient access throughout the function
    let TrajectoryConfig {
        velocity,
        angle,
        bc,
        mass,
        diameter,
        drag_model,
        max_range,
        time_step,
        temperature,
        pressure,
        humidity,
        altitude,
        wind_speed,
        wind_direction,
        ref wind_segments,
        output,
        full,
        units,
        sight_height,
        bore_height,
        ignore_ground_impact,
        use_bc_segments,
        use_cluster_bc,
        ref bc_table_segments,
        enable_magnus,
        enable_coriolis,
        enable_spin_drift,
        enable_aerodynamic_jump,
        enable_wind_shear,
        enable_pitch_damping,
        enable_precession,
        sample_trajectory,
        sample_interval,
        use_rk4,
        use_rk45,
        twist_rate,
        twist_right,
        latitude,
        shooting_angle,
        use_powder_sensitivity,
        powder_temp_sensitivity,
        powder_temp,
        ref pdf_metadata,
    } = *config;

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
        enable_magnus,
        enable_coriolis,
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
        // Use BC segments if: 1) table generated them (parameter), 2) --use-bc-segments flag, 3) neither
        use_bc_segments: use_bc_segments || bc_table_segments.is_some(),
        bc_segments: None,
        bc_segments_data: if let Some(ref segments) = *bc_table_segments {
            // Use segments passed from BC table (highest priority)
            Some(segments.clone())
        } else if use_bc_segments {
            // Generate BC segments automatically from estimator
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
        enable_aerodynamic_jump,
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

    // Downrange-segmented wind (overrides the scalar wind when present).
    if !wind_segments.is_empty() {
        if enable_wind_shear {
            return Err("--wind-segment cannot be combined with --enable-wind-shear \
                (downrange segments + altitude shear is not yet a defined model)"
                .into());
        }
        // Note when a non-zero scalar wind is also set, since segments take precedence.
        if wind_speed != 0.0 {
            eprintln!(
                "note: --wind-segment overrides --wind-speed/--wind-direction (scalar wind ignored)"
            );
        }
        // Warn if the segments don't cover the whole trajectory (wind is zero beyond
        // the last segment's until-distance). A 1 m epsilon avoids float-noise warnings
        // when a segment is set right at max_range.
        let coverage_m = wind_segments
            .iter()
            .map(|(_, _, until)| *until)
            .fold(0.0_f64, f64::max);
        if coverage_m < max_range - 1.0 {
            let (cov_disp, max_disp, unit) = match units {
                UnitSystem::Imperial => (coverage_m / 0.9144, max_range / 0.9144, "yd"),
                UnitSystem::Metric => (coverage_m, max_range, "m"),
            };
            eprintln!(
                "warning: wind segments cover only {cov_disp:.0} {unit} of the {max_disp:.0} {unit} \
                trajectory; wind is treated as zero beyond {cov_disp:.0} {unit}"
            );
        }
        solver.set_wind_segments(wind_segments.clone());
    }

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
                            // Output contract is unchanged: the `x` field is lateral
                            // (drift), `z` is downrange. With McCoy internals these map
                            // to position.z (lateral) and position.x (downrange).
                            x: p.position.z,
                            y: p.position.y,
                            z: p.position.x,
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
                    // Drop/drift use the same small-deflection units as the table/API
                    // (inches for Imperial, meters for Metric), NOT the range unit (MBA-950).
                    let defl_unit = match units {
                        UnitSystem::Metric => "m",
                        UnitSystem::Imperial => "in",
                    };
                    println!(
                        "distance_{},drop_{},drift_{},velocity_{},energy_{},time_s",
                        dist_unit, defl_unit, defl_unit, vel_unit, energy_unit
                    );
                    for s in sampled {
                        let distance = UnitConverter::distance_from_metric(s.distance_m, units);
                        // Imperial drop/drift in inches (meters * 39.3701); Metric in meters.
                        let (drop, drift) = match units {
                            UnitSystem::Imperial => (s.drop_m * 39.3701, s.wind_drift_m * 39.3701),
                            UnitSystem::Metric => (s.drop_m, s.wind_drift_m),
                        };
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
                        // Output contract unchanged: x column is lateral, z is downrange.
                        // McCoy internals: lateral=position.z, downrange=position.x.
                        let x = UnitConverter::distance_from_metric(p.position.z, units);
                        let y = UnitConverter::distance_from_metric(p.position.y, units);
                        let z = UnitConverter::distance_from_metric(p.position.x, units);
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
                let last_range = last_point.position.x;

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
                        let x_display = UnitConverter::distance_from_metric(p.position.x, units); // X column = downrange (position.x; McCoy frame)
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

        #[cfg(feature = "pdf")]
        OutputFormat::Pdf => {
            // PDF output requires metadata and output file
            let pdf_meta = pdf_metadata.as_ref().ok_or("PDF output requires metadata (use --target-speed, --powder, --bullet-name, etc.)")?;
            let output_path = pdf_meta.output_file.as_ref().ok_or("PDF output requires --output-file")?;

            // Get sampled trajectory points (required for dope card)
            let sampled = result.sampled_points.as_ref().ok_or(
                "PDF output requires --sample-trajectory flag for trajectory data"
            )?;

            // Calculate density altitude
            let da_ft = calculate_density_altitude(
                pdf_meta.altitude_ft,
                pdf_meta.pressure_inhg,
                pdf_meta.temperature_f,
            );

            // Build PDF config
            let pdf_config = DopeCardConfig {
                rifle_name: pdf_meta.rifle_name.clone(),
                location: pdf_meta.location_name.clone(),
                density_altitude_ft: da_ft,
                pressure_inhg: pdf_meta.pressure_inhg,
                pressure_hpa: pdf_meta.pressure_inhg * 33.8639,
                temperature_f: pdf_meta.temperature_f,
                altitude_ft: pdf_meta.altitude_ft,
                wind_speed_mph: pdf_meta.wind_speed_mph,
                target_speed_mph: pdf_meta.target_speed_mph,
                solver_mode: if cfg!(feature = "online") { "online".to_string() } else { "offline".to_string() },
                powder: pdf_meta.powder.clone(),
                bullet: pdf_meta.bullet_name.clone(),
                weight_gr: pdf_meta.weight_gr,
                bc,
                drag_model: match drag_model {
                    DragModelArg::G1 => "G1".to_string(),
                    DragModelArg::G7 => "G7".to_string(),
                },
                velocity_fps: pdf_meta.velocity_fps,
                font_scale: pdf_meta.font_scale,
                bold_data: pdf_meta.bold_data,
            };

            // Convert sampled trajectory to dope card rows
            // Filter from zero range onwards (typically 100 yards)
            let rows: Vec<DopeCardRow> = sampled.iter()
                .filter(|s| s.distance_m >= UnitConverter::distance_to_metric(100.0, UnitSystem::Imperial))
                .map(|s| {
                    // Convert distance to yards for range
                    let range_yd = UnitConverter::distance_from_metric(s.distance_m, UnitSystem::Imperial);
                    // Convert drop to yards (s.drop_m is already in meters, positive = below line of sight)
                    let drop_yd = s.drop_m / 0.9144; // meters to yards
                    // Convert drift to yards
                    let drift_yd = s.wind_drift_m / 0.9144;

                    DopeCardRow {
                        range_yd: range_yd.round() as u32,
                        // Drop MIL: positive = dial up (bullet below LOS). drop_m is already
                        // positive-below-LOS (sample_trajectory: los_y - y_interp), matching the
                        // come-up / range tables, so do NOT negate it (this column was sign-flipped).
                        drop_mil: yards_to_mil(drop_yd, range_yd),
                        // Wind MIL: positive = dial right for wind from right
                        wind_mil: yards_to_mil(drift_yd, range_yd),
                        // Lead MIL for moving target
                        lead_mil: calculate_lead_mil(pdf_meta.target_speed_mph, s.time_s, range_yd),
                    }
                })
                .collect();

            if rows.is_empty() {
                return Err("No trajectory data available for PDF output. Ensure --sample-trajectory is enabled and range > 100 yards.".into());
            }

            // Generate PDF
            let pdf_bytes = pdf_dope_card::generate_dope_card_pdf(&pdf_config, &rows)?;

            // Write to file
            std::fs::write(output_path, &pdf_bytes)?;
            eprintln!("PDF dope card written to: {}", output_path.display());
            eprintln!("  {} ranges from {} to {} yards", rows.len(), rows.first().map(|r| r.range_yd).unwrap_or(0), rows.last().map(|r| r.range_yd).unwrap_or(0));
        }
        #[cfg(not(feature = "pdf"))]
        OutputFormat::Pdf => {
            return Err("PDF output requires building with the 'pdf' feature (e.g. cargo build --features pdf)".into());
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

        OutputFormat::Pdf => {
            eprintln!("Error: PDF output is not supported for API trajectory results.");
            eprintln!("Hint: Use local calculation (without --online) for PDF dope card generation.");
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
    temperature: f64,
    pressure: f64,
    humidity: f64,
    altitude: f64,
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
        temperature,
        pressure,
        humidity,
        altitude,
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

    // Calculate zero angle with atmospheric conditions
    let zero_angle =
        ballistics_engine::calculate_zero_angle_with_conditions(
            inputs.clone(),
            target_distance,
            target_height,
            WindConditions::default(),
            atmosphere.clone(),
        )?;

    // Calculate trajectory at zero angle to get additional info
    let mut zeroed_inputs = inputs;
    zeroed_inputs.muzzle_angle = zero_angle;

    let solver = TrajectorySolver::new(zeroed_inputs, Default::default(), atmosphere);
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
                    // Preserves pre-existing behavior: this read the lateral component
                    // (old position.x), which is now position.z. NOTE: this looks like a
                    // latent bug (PBR should be downrange = position.x); fix separately.
                    .map(|p| p.position.z)
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

        OutputFormat::Pdf => {
            eprintln!("Error: PDF output is not supported for zero calculation.");
            eprintln!("Hint: Use --output json, --output csv, or --output table instead.");
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

        OutputFormat::Pdf => {
            eprintln!("Error: PDF output is not supported for BC segment generation.");
            eprintln!("Hint: Use --output json, --output csv, or --output table instead.");
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

    let mut solver = TrajectorySolver::new(inputs, Default::default(), Default::default());
    // Bound the verification range like the estimator does (estimate_bc_from_trajectory caps at
    // last_point * 1.5); the solver default of 1000 m otherwise gives a bogus 100% error for
    // long-range inputs (no point reaches distance2, so calc_drop2 falls back to 0.0). Metric.
    solver.set_max_range(distance1.max(distance2) * 1.5);
    let trajectory = solver.solve()?;

    // Find drops at the specified distances (X is downrange)
    let calc_drop1 = trajectory
        .points
        .iter()
        .find(|p| p.position.x >= distance1)
        .map(|p| -p.position.y)
        .unwrap_or(0.0);

    let calc_drop2 = trajectory
        .points
        .iter()
        .find(|p| p.position.x >= distance2)
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

        OutputFormat::Pdf => {
            eprintln!("Error: PDF output is not supported for BC estimation.");
            eprintln!("Hint: Use --output json, --output csv, or --output table instead.");
        }
    }

    Ok(())
}

/// Display velocity truing results in the specified format
fn display_true_velocity_result(
    effective_vel: f64,
    vel_unit: &str,
    velocity_adjustment: Option<f64>,
    adjustment_percent: Option<f64>,
    confidence: &str,
    iterations: i32,
    final_error_mil: f64,
    calculated_drop_mil: f64,
    measured_drop: f64,
    units: UnitSystem,
    output: OutputFormat,
    used_bc_table: bool,
) -> Result<(), Box<dyn Error>> {
    match output {
        OutputFormat::Json => {
            let mut json_output = serde_json::json!({
                "effective_velocity": effective_vel,
                "velocity_unit": vel_unit,
                "confidence": confidence,
                "iterations": iterations,
                "final_error_mil": final_error_mil,
                "calculated_drop_mil": calculated_drop_mil,
                "measured_drop_mil": measured_drop,
                "used_bc_table": used_bc_table,
            });
            if let Some(adj) = velocity_adjustment {
                let adj_display = match units {
                    UnitSystem::Imperial => adj,
                    UnitSystem::Metric => adj * 0.3048,
                };
                json_output["velocity_adjustment"] = serde_json::json!(adj_display);
            }
            if let Some(pct) = adjustment_percent {
                json_output["adjustment_percent"] = serde_json::json!(pct);
            }
            println!("{}", serde_json::to_string_pretty(&json_output)?);
        }
        OutputFormat::Csv => {
            println!("effective_velocity,unit,adjustment,adjustment_pct,confidence,iterations,final_error_mil,calculated_drop_mil,used_bc_table");
            println!("{:.1},{},{},{},{},{},{:.4},{:.2},{}",
                effective_vel,
                vel_unit,
                velocity_adjustment.map(|v| {
                    let adj = match units {
                        UnitSystem::Imperial => v,
                        UnitSystem::Metric => v * 0.3048,
                    };
                    format!("{:.1}", adj)
                }).unwrap_or_default(),
                adjustment_percent.map(|v| format!("{:.2}", v)).unwrap_or_default(),
                confidence,
                iterations,
                final_error_mil,
                calculated_drop_mil,
                used_bc_table,
            );
        }
        OutputFormat::Table => {
            println!();
            println!("╔════════════════════════════════════════════════════════════╗");
            println!("║           VELOCITY TRUING RESULTS                          ║");
            println!("╠════════════════════════════════════════════════════════════╣");
            println!("║  Effective Muzzle Velocity: {:>8.1} {:>4}                 ║", effective_vel, vel_unit);
            if let Some(adj) = velocity_adjustment {
                let adj_display = match units {
                    UnitSystem::Imperial => adj,
                    UnitSystem::Metric => adj * 0.3048,
                };
                println!("║  Adjustment from Chrono:    {:>+8.1} {:>4}                 ║", adj_display, vel_unit);
                if let Some(pct) = adjustment_percent {
                    println!("║  Adjustment Percentage:     {:>+8.2}%                      ║", pct);
                }
            }
            println!("╠════════════════════════════════════════════════════════════╣");
            println!("║  Confidence:                {:>8}                        ║", confidence);
            println!("║  Iterations:                {:>8}                        ║", iterations);
            println!("║  Final Error:               {:>8.4} MIL                  ║", final_error_mil);
            println!("║  Calculated Drop:           {:>8.2} MIL                  ║", calculated_drop_mil);
            println!("║  Measured Drop:             {:>8.2} MIL                  ║", measured_drop);
            if used_bc_table {
                println!("╠════════════════════════════════════════════════════════════╣");
                println!("║  BC5D Table:                     Yes                       ║");
            }
            println!("╚════════════════════════════════════════════════════════════╝");
            println!();
        }

        OutputFormat::Pdf => {
            eprintln!("Error: PDF output is not supported for velocity truing results.");
            eprintln!("Hint: Use --output json, --output csv, or --output table instead.");
        }
    }

    Ok(())
}

/// Calculate drop at a given muzzle velocity using trajectory solver
/// Returns drop in MILs at the target range
fn calculate_drop_at_velocity(
    velocity_fps: f64,
    bc: f64,
    drag_model: DragModelArg,
    mass_gr: f64,
    diameter_in: f64,
    zero_distance_yd: f64,
    range_yd: f64,
    sight_height_in: f64,
    temperature_f: f64,
    pressure_inhg: f64,
    humidity: f64,
    altitude_ft: f64,
    bc_segments: &Option<Vec<BCSegmentData>>,
) -> Result<f64, Box<dyn Error>> {
    // Convert to SI units
    let velocity_ms = velocity_fps * 0.3048;
    let mass_kg = mass_gr * 0.0000647989;
    let diameter_m = diameter_in * 0.0254;
    let zero_m = zero_distance_yd * 0.9144;
    let range_m = range_yd * 0.9144;
    let sight_height_m = sight_height_in * 0.0254;
    let altitude_m = altitude_ft * 0.3048;
    let temperature_c = (temperature_f - 32.0) * 5.0 / 9.0;
    let pressure_hpa = pressure_inhg * 33.8639; // Convert inHg to hPa

    let drag_model_enum = match drag_model {
        DragModelArg::G1 => DragModel::G1,
        DragModelArg::G7 => DragModel::G7,
    };

    // Create base inputs - match defaults used by trajectory command
    let mut inputs = BallisticInputs {
        muzzle_velocity: velocity_ms,
        bc_value: bc,
        bc_type: drag_model_enum,
        bullet_mass: mass_kg,
        bullet_diameter: diameter_m,
        bullet_length: diameter_m * 4.5, // Approximate
        sight_height: sight_height_m,
        target_distance: range_m + 100.0, // Overshoot to ensure we have data
        use_bc_segments: bc_segments.is_some(),
        bc_segments_data: bc_segments.clone(),
        use_rk4: true,
        muzzle_angle: 0.0, // Will be set by zero angle calculation
        ..Default::default() // Uses muzzle_height: 0.0 by default
    };

    // Set up atmospheric conditions
    // AtmosphericConditions expects: temperature in Celsius, pressure in hPa, humidity 0-100, altitude in meters
    let atmosphere = AtmosphericConditions {
        temperature: temperature_c,
        pressure: pressure_hpa,
        humidity, // Already 0-100 from input
        altitude: altitude_m,
        ..Default::default()
    };

    let wind = WindConditions::default();

    // Calculate zero angle for the zero distance
    // Target height is sight_height because the bullet must cross the LOS at zero distance
    // The LOS is at y = sight_height (sight is above bore by sight_height)
    // So the bullet (starting at y = 0 = bore level) must rise to y = sight_height at zero distance
    let zero_angle = ballistics_engine::calculate_zero_angle_with_conditions(
        inputs.clone(),
        zero_m,
        sight_height_m, // target height at zero distance (LOS height)
        wind.clone(),
        atmosphere.clone(),
    )?;

    // Set the calculated zero angle
    inputs.muzzle_angle = zero_angle;

    // Create solver and solve
    let mut solver = TrajectorySolver::new(inputs, wind, atmosphere);
    solver.set_max_range(range_m + 100.0);
    solver.set_time_step(0.0001);

    let result = solver.solve()?;

    // Find the point at target range
    let target_point = result.points.iter()
        .find(|p| p.position.x >= range_m)
        .ok_or("Trajectory didn't reach target range")?;

    // Calculate drop relative to line of sight (LOS)
    // Using the same formula as the Flask API (true_velocity.py):
    //   los_height = z * tan(barrel_angle) - sight_height
    //   drop = los_height - bullet_y
    //
    // The barrel_angle is the zero_angle we calculated earlier.
    // This formula accounts for the angled barrel and scope offset.
    let barrel_angle = zero_angle;
    let z = target_point.position.x;
    let bullet_y = target_point.position.y;

    // LOS height at range z
    let los_height = z * barrel_angle.tan() - sight_height_m;

    // Drop = LOS height - bullet position (positive = below LOS)
    let drop_m = los_height - bullet_y;

    // Convert to MILs: mil = (drop_inches / 36 / range_yards) * 1000
    // Or equivalently: mil = (drop_m / range_m) * 1000
    let drop_mil = (drop_m / z) * 1000.0;

    Ok(drop_mil)
}

/// Result of local velocity truing calculation
struct TrueVelocityLocalResult {
    effective_velocity_fps: f64,
    iterations: i32,
    final_error_mil: f64,
    calculated_drop_mil: f64,
    confidence: String,
}

/// Calculate true velocity using local binary search
/// Returns (effective_velocity_fps, iterations, final_error_mil, calculated_drop_mil)
fn calculate_true_velocity_local(
    measured_drop_mil: f64,
    range_yd: f64,
    bc: f64,
    drag_model: DragModelArg,
    mass_gr: f64,
    diameter_in: f64,
    zero_distance_yd: f64,
    sight_height_in: f64,
    temperature_f: f64,
    pressure_inhg: f64,
    humidity: f64,
    altitude_ft: f64,
    bc_segments: &Option<Vec<BCSegmentData>>,
) -> Result<TrueVelocityLocalResult, Box<dyn Error>> {
    // Binary search between velocity bounds
    let mut velocity_low = 1500.0;
    let mut velocity_high = 4500.0;
    let tolerance_mil = 0.01; // 0.01 MIL tolerance
    let max_iterations = 50;

    let mut iterations = 0;
    let mut last_error = 0.0;
    let mut last_calculated_drop = 0.0;

    for i in 0..max_iterations {
        iterations = i + 1;
        let test_velocity = (velocity_low + velocity_high) / 2.0;

        // Run trajectory at test velocity
        let calculated_drop_mil = calculate_drop_at_velocity(
            test_velocity,
            bc,
            drag_model,
            mass_gr,
            diameter_in,
            zero_distance_yd,
            range_yd,
            sight_height_in,
            temperature_f,
            pressure_inhg,
            humidity,
            altitude_ft,
            bc_segments,
        )?;

        last_calculated_drop = calculated_drop_mil;
        let error = calculated_drop_mil - measured_drop_mil;
        last_error = error;

        if error.abs() < tolerance_mil {
            // Converged
            let confidence = if error.abs() < 0.005 {
                "high"
            } else {
                "medium"
            };

            return Ok(TrueVelocityLocalResult {
                effective_velocity_fps: test_velocity,
                iterations,
                final_error_mil: error,
                calculated_drop_mil,
                confidence: confidence.to_string(),
            });
        }

        // Higher calculated drop = bullet is slower = need higher velocity
        // Lower calculated drop = bullet is faster = need lower velocity
        if calculated_drop_mil > measured_drop_mil {
            // Bullet dropping more than observed = slower than actual
            // Need higher velocity
            velocity_low = test_velocity;
        } else {
            // Bullet dropping less = faster than actual
            // Need lower velocity
            velocity_high = test_velocity;
        }

        // Check for convergence issues
        if (velocity_high - velocity_low).abs() < 0.5 {
            break;
        }
    }

    // Did not converge within tolerance, return best estimate
    let final_velocity = (velocity_low + velocity_high) / 2.0;
    let confidence = if last_error.abs() < 0.1 { "medium" } else { "low" };

    Ok(TrueVelocityLocalResult {
        effective_velocity_fps: final_velocity,
        iterations,
        final_error_mil: last_error,
        calculated_drop_mil: last_calculated_drop,
        confidence: confidence.to_string(),
    })
}

// ============================================================================
// MPBR, Come-Ups, Wind Card Handler Functions
// ============================================================================

/// Shared: build BallisticInputs + atmosphere + wind from common parameters (all in metric)
fn build_trajectory_components(
    velocity: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    drag_model: DragModelArg,
    sight_height: f64,
    temperature: f64,
    pressure: f64,
    humidity: f64,
    altitude: f64,
    wind_speed: f64,
    wind_direction: f64,
    max_range: f64,
    sample_interval: f64,
) -> (BallisticInputs, WindConditions, AtmosphericConditions) {
    let drag_model_enum = match drag_model {
        DragModelArg::G1 => DragModel::G1,
        DragModelArg::G7 => DragModel::G7,
    };

    let inputs = BallisticInputs {
        bc_value: bc,
        bc_type: drag_model_enum,
        bullet_mass: mass,
        muzzle_velocity: velocity,
        bullet_diameter: diameter,
        bullet_length: diameter * 4.5,
        muzzle_angle: 0.0,
        target_distance: max_range,
        sight_height,
        altitude,
        temperature,
        pressure,
        humidity,
        wind_speed,
        wind_angle: wind_direction,
        use_rk4: true,             // Required for non-Euler solver
        use_adaptive_rk45: true,   // Use RK45 adaptive (default solver)
        enable_trajectory_sampling: true,
        sample_interval,
        caliber_inches: diameter / 0.0254,
        weight_grains: mass / 0.00006479891,
        twist_rate: 12.0,
        is_twist_right: true,
        ..Default::default()
    };

    // wind_direction is in degrees (matching BallisticInputs convention)
    let wind = WindConditions {
        speed: wind_speed,
        direction: wind_direction.to_radians(), // WindConditions expects radians
        ..Default::default()
    };

    let atmosphere = AtmosphericConditions {
        temperature,
        pressure,
        humidity,
        altitude,
        ..Default::default()
    };

    (inputs, wind, atmosphere)
}

/// Run a trajectory and return sampled points at the given zero angle
fn run_sampled_trajectory(
    velocity: f64, bc: f64, mass: f64, diameter: f64,
    drag_model: DragModelArg, sight_height: f64,
    temperature: f64, pressure: f64, humidity: f64, altitude: f64,
    wind_speed: f64, wind_direction: f64,
    max_range: f64, sample_interval: f64,
    zero_angle_rad: f64,
) -> Result<Vec<trajectory_sampling::TrajectorySample>, Box<dyn Error>> {
    let (mut inputs, wind, atmosphere) = build_trajectory_components(
        velocity, bc, mass, diameter, drag_model, sight_height,
        temperature, pressure, humidity, altitude,
        wind_speed, wind_direction, max_range, sample_interval,
    );
    inputs.muzzle_angle = zero_angle_rad;

    let mut solver = TrajectorySolver::new(inputs, wind, atmosphere);
    solver.set_max_range(max_range);
    solver.set_time_step(0.001);
    let result = solver.solve()?;

    Ok(result.sampled_points.unwrap_or_default())
}

/// Resolve bullet parameters: CLI arg overrides profile value
fn resolve_param(cli_val: Option<f64>, profile: &Option<ProfileData>, getter: fn(&ProfileData) -> f64) -> Option<f64> {
    cli_val.or_else(|| profile.as_ref().map(getter))
}

/// MPBR handler
fn handle_mpbr(
    velocity: f64, bc: f64, mass: f64, diameter: f64,
    drag_model: DragModelArg,
    vital_zone: f64, // in user units (inches or cm)
    sight_height: f64, // in user units
    temperature: f64, pressure: f64, humidity: f64, altitude: f64,
    units: UnitSystem,
    output: OutputFormat,
) -> Result<(), Box<dyn Error>> {
    // Convert everything to metric
    let velocity_m = UnitConverter::velocity_to_metric(velocity, units);
    let mass_kg = UnitConverter::mass_to_metric(mass, units);
    let diameter_m = UnitConverter::diameter_to_metric(diameter, units);
    let sight_height_m = UnitConverter::sight_height_to_metric(sight_height, units);
    let temperature_c = UnitConverter::temperature_to_metric(temperature, units);
    let pressure_hpa = UnitConverter::pressure_to_metric(pressure, units);
    let altitude_m = UnitConverter::altitude_to_metric(altitude, units);

    // Vital zone to meters
    let vital_zone_m = match units {
        UnitSystem::Imperial => vital_zone * 0.0254, // inches to meters
        UnitSystem::Metric => vital_zone * 0.01,     // cm to meters
    };
    let half_vital_m = vital_zone_m / 2.0;

    // Binary search on zero distance to find the one where max ordinate ≈ half vital zone
    let mut zero_low_m = UnitConverter::distance_to_metric(50.0, units);
    let mut zero_high_m = UnitConverter::distance_to_metric(2000.0, units);
    let tolerance_m = 0.001; // ~0.04 inches
    let max_iter = 60;

    let mut best_zero_m = 0.0;
    let mut best_max_ord_m = 0.0;
    let mut best_max_ord_dist_m = 0.0;

    for _ in 0..max_iter {
        let test_zero_m = (zero_low_m + zero_high_m) / 2.0;

        // Create inputs for zero calculation at test_zero_m
        let drag_model_enum = match drag_model {
            DragModelArg::G1 => DragModel::G1,
            DragModelArg::G7 => DragModel::G7,
        };

        let zero_inputs = BallisticInputs {
            bc_value: bc,
            bc_type: drag_model_enum,
            bullet_mass: mass_kg,
            muzzle_velocity: velocity_m,
            bullet_diameter: diameter_m,
            bullet_length: diameter_m * 4.5,
            sight_height: sight_height_m,
            use_rk4: true,
            ..Default::default()
        };

        let atmosphere = AtmosphericConditions {
            temperature: temperature_c,
            pressure: pressure_hpa,
            humidity,
            altitude: altitude_m,
            ..Default::default()
        };

        let zero_angle = match ballistics_engine::calculate_zero_angle_with_conditions(
            zero_inputs, test_zero_m, sight_height_m,
            WindConditions::default(), atmosphere,
        ) {
            Ok(a) => a,
            Err(_) => {
                zero_high_m = test_zero_m;
                continue;
            }
        };

        // Run trajectory with this zero, sampling every ~1 yard
        let samples = match run_sampled_trajectory(
            velocity_m, bc, mass_kg, diameter_m, drag_model, sight_height_m,
            temperature_c, pressure_hpa, humidity, altitude_m,
            0.0, 0.0,
            test_zero_m * 1.5, // max range past zero
            UnitConverter::distance_to_metric(1.0, UnitSystem::Imperial), // ~1 yd sample interval
            zero_angle,
        ) {
            Ok(s) => s,
            Err(_) => {
                zero_high_m = test_zero_m;
                continue;
            }
        };

        // Find max ordinate (highest point above LOS)
        // drop_m in samples is relative to LOS (negative = below LOS, but in the coordinate system
        // the y-position grows then drops, so we need to find max height above LOS)
        // Actually, drop_m is the height of the bullet relative to line of sight at that distance
        // At the zero crossing: drop_m ≈ 0. Between muzzle and zero, drop_m > 0 (above LOS).
        // Beyond zero, drop_m < 0 (below LOS).
        // NOTE: In the sampling system, drop_m is defined as bullet_y - los_y, so positive = above LOS
        let mut max_ord_m: f64 = 0.0;
        let mut max_ord_dist_m: f64 = 0.0;
        for s in &samples {
            // drop_m is the drop value; when bullet is above LOS, this value represents height above LOS
            // In the sampling code, drop_m = bullet_y - line_of_sight_y
            // When zeroed, the bullet rises above LOS then falls back.
            // A negative drop_m means below LOS. We want the max POSITIVE value (highest above LOS).
            let height_above_los = -s.drop_m; // drop_m is negative when above LOS in standard convention
            if height_above_los > max_ord_m {
                max_ord_m = height_above_los;
                max_ord_dist_m = s.distance_m;
            }
        }

        // Check: if max_ord is 0 but we have points above bore, try the other sign convention
        if max_ord_m < 0.001 {
            // Try positive drop_m as "above LOS"
            for s in &samples {
                if s.drop_m > max_ord_m {
                    max_ord_m = s.drop_m;
                    max_ord_dist_m = s.distance_m;
                }
            }
        }

        best_zero_m = test_zero_m;
        best_max_ord_m = max_ord_m;
        best_max_ord_dist_m = max_ord_dist_m;

        let diff = max_ord_m - half_vital_m;
        if diff.abs() < tolerance_m {
            break;
        }

        if max_ord_m > half_vital_m {
            // Max ordinate too high → zero closer (reduce zero distance)
            zero_high_m = test_zero_m;
        } else {
            // Max ordinate too low → zero farther (increase zero distance)
            zero_low_m = test_zero_m;
        }
    }

    // Now run the final trajectory at optimal zero to find MPBR, near zero, far zero
    let drag_model_enum = match drag_model {
        DragModelArg::G1 => DragModel::G1,
        DragModelArg::G7 => DragModel::G7,
    };

    let final_inputs = BallisticInputs {
        bc_value: bc,
        bc_type: drag_model_enum,
        bullet_mass: mass_kg,
        muzzle_velocity: velocity_m,
        bullet_diameter: diameter_m,
        bullet_length: diameter_m * 4.5,
        sight_height: sight_height_m,
        use_rk4: true,
        ..Default::default()
    };

    let atmosphere = AtmosphericConditions {
        temperature: temperature_c,
        pressure: pressure_hpa,
        humidity,
        altitude: altitude_m,
        ..Default::default()
    };

    let final_zero_angle = ballistics_engine::calculate_zero_angle_with_conditions(
        final_inputs, best_zero_m, sight_height_m,
        WindConditions::default(), atmosphere,
    )?;

    let final_samples = run_sampled_trajectory(
        velocity_m, bc, mass_kg, diameter_m, drag_model, sight_height_m,
        temperature_c, pressure_hpa, humidity, altitude_m,
        0.0, 0.0,
        best_zero_m * 2.0,
        UnitConverter::distance_to_metric(1.0, UnitSystem::Imperial),
        final_zero_angle,
    )?;

    // Find near zero crossing (first time trajectory crosses from below LOS to above)
    // Find far zero crossing (where trajectory drops back to LOS)
    // Find MPBR (where drop exceeds -half_vital_m)
    let mut near_zero_m: f64 = 0.0;
    let mut far_zero_m: f64 = 0.0;
    let mut mpbr_m: f64 = 0.0;
    let mut impact_vel_mps: f64 = 0.0;
    let mut impact_energy_j: f64 = 0.0;
    let mut found_near = false;
    let mut found_far = false;

    // Determine sign convention by checking if samples near the midpoint are positive or negative
    // The bullet should be ABOVE LOS in the middle of the trajectory
    let sign_flip = if final_samples.len() > 10 {
        let mid = &final_samples[final_samples.len() / 3];
        mid.drop_m < 0.0 // if drop_m is negative in the middle, we need to flip sign
    } else {
        false
    };

    for i in 1..final_samples.len() {
        let prev_drop = if sign_flip { -final_samples[i-1].drop_m } else { final_samples[i-1].drop_m };
        let curr_drop = if sign_flip { -final_samples[i].drop_m } else { final_samples[i].drop_m };

        // Near zero: trajectory goes from negative/zero to positive (crossing LOS upward)
        if !found_near && prev_drop <= 0.0 && curr_drop > 0.0 && final_samples[i].distance_m > 5.0 {
            // Interpolate
            let frac = (-prev_drop) / (curr_drop - prev_drop);
            near_zero_m = final_samples[i-1].distance_m + frac * (final_samples[i].distance_m - final_samples[i-1].distance_m);
            found_near = true;
        }

        // Far zero: trajectory goes from positive to negative (crossing LOS downward)
        if found_near && !found_far && prev_drop > 0.0 && curr_drop <= 0.0 {
            let frac = prev_drop / (prev_drop - curr_drop);
            far_zero_m = final_samples[i-1].distance_m + frac * (final_samples[i].distance_m - final_samples[i-1].distance_m);
            found_far = true;
        }

        // MPBR: where drop goes below -half_vital_m
        if found_far && curr_drop < -half_vital_m && mpbr_m == 0.0 {
            let frac = (-half_vital_m - prev_drop) / (curr_drop - prev_drop);
            mpbr_m = final_samples[i-1].distance_m + frac * (final_samples[i].distance_m - final_samples[i-1].distance_m);
            impact_vel_mps = final_samples[i].velocity_mps;
            impact_energy_j = final_samples[i].energy_j;
        }
    }

    // Convert outputs for display
    let mpbr_display = UnitConverter::distance_from_metric(mpbr_m, units);
    let optimal_zero_display = UnitConverter::distance_from_metric(best_zero_m, units);
    let near_zero_display = UnitConverter::distance_from_metric(near_zero_m, units);
    let far_zero_display = UnitConverter::distance_from_metric(far_zero_m, units);
    let max_ord_display = match units {
        UnitSystem::Imperial => best_max_ord_m / 0.0254, // meters to inches
        UnitSystem::Metric => best_max_ord_m * 100.0,    // meters to cm
    };
    let max_ord_dist_display = UnitConverter::distance_from_metric(best_max_ord_dist_m, units);
    let impact_vel_display = UnitConverter::velocity_from_metric(impact_vel_mps, units);
    let impact_energy_display = UnitConverter::energy_from_metric(impact_energy_j, units);

    let (dist_unit, vel_unit, energy_unit, size_unit) = match units {
        UnitSystem::Imperial => ("yd", "fps", "ft-lb", "in"),
        UnitSystem::Metric => ("m", "m/s", "J", "cm"),
    };

    let vital_zone_display = match units {
        UnitSystem::Imperial => vital_zone, // already in inches
        UnitSystem::Metric => vital_zone,   // already in cm
    };

    match output {
        OutputFormat::Json => {
            let json = serde_json::json!({
                "mpbr": mpbr_display,
                "mpbr_unit": dist_unit,
                "optimal_zero": optimal_zero_display,
                "near_zero": near_zero_display,
                "far_zero": far_zero_display,
                "max_ordinate": max_ord_display,
                "max_ordinate_unit": size_unit,
                "max_ordinate_distance": max_ord_dist_display,
                "impact_velocity": impact_vel_display,
                "impact_velocity_unit": vel_unit,
                "impact_energy": impact_energy_display,
                "impact_energy_unit": energy_unit,
                "vital_zone": vital_zone_display,
                "vital_zone_unit": size_unit,
            });
            println!("{}", serde_json::to_string_pretty(&json)?);
        }
        OutputFormat::Csv => {
            println!("metric,value,unit");
            println!("mpbr,{:.1},{}", mpbr_display, dist_unit);
            println!("optimal_zero,{:.1},{}", optimal_zero_display, dist_unit);
            println!("near_zero,{:.1},{}", near_zero_display, dist_unit);
            println!("far_zero,{:.1},{}", far_zero_display, dist_unit);
            println!("max_ordinate,{:.1},{}", max_ord_display, size_unit);
            println!("max_ordinate_distance,{:.1},{}", max_ord_dist_display, dist_unit);
            println!("impact_velocity,{:.0},{}", impact_vel_display, vel_unit);
            println!("impact_energy,{:.0},{}", impact_energy_display, energy_unit);
        }
        OutputFormat::Table | OutputFormat::Pdf => {
            println!();
            println!("MPBR Analysis (vital zone: {:.1} {})", vital_zone_display, size_unit);
            println!("╔════════════════════════════════════════╗");
            println!("║  MPBR:            {:>6.0} {:<14}║", mpbr_display, dist_unit);
            println!("║  Optimal zero:    {:>6.0} {:<14}║", optimal_zero_display, dist_unit);
            println!("║  Near zero:       {:>6.0} {:<14}║", near_zero_display, dist_unit);
            println!("║  Far zero:        {:>6.0} {:<14}║", far_zero_display, dist_unit);
            println!("║  Max ordinate:    {:>6.1} {} at {:.0} {}  ║", max_ord_display, size_unit, max_ord_dist_display, dist_unit);
            println!("║  Impact velocity: {:>6.0} {:<14}║", impact_vel_display, vel_unit);
            println!("║  Impact energy:   {:>6.0} {:<14}║", impact_energy_display, energy_unit);
            println!("╚════════════════════════════════════════╝");
        }
    }

    Ok(())
}

/// Come-ups handler
fn handle_come_ups(
    velocity: f64, bc: f64, mass: f64, diameter: f64,
    drag_model: DragModelArg,
    zero_distance: f64, start: f64, end: f64, step: f64,
    adjustment_unit: AdjustmentUnit,
    sight_height: f64,
    temperature: f64, pressure: f64, humidity: f64, altitude: f64,
    wind_speed: f64, wind_direction: f64,
    units: UnitSystem,
    output: OutputFormat,
) -> Result<(), Box<dyn Error>> {
    // Convert to metric
    let velocity_m = UnitConverter::velocity_to_metric(velocity, units);
    let mass_kg = UnitConverter::mass_to_metric(mass, units);
    let diameter_m = UnitConverter::diameter_to_metric(diameter, units);
    let sight_height_m = UnitConverter::sight_height_to_metric(sight_height, units);
    let zero_distance_m = UnitConverter::distance_to_metric(zero_distance, units);
    let temperature_c = UnitConverter::temperature_to_metric(temperature, units);
    let pressure_hpa = UnitConverter::pressure_to_metric(pressure, units);
    let altitude_m = UnitConverter::altitude_to_metric(altitude, units);
    let wind_speed_m = UnitConverter::wind_to_metric(wind_speed, units);
    let end_m = UnitConverter::distance_to_metric(end, units);
    let sample_m = UnitConverter::distance_to_metric(step, units);

    // Calculate zero angle
    let drag_model_enum = match drag_model {
        DragModelArg::G1 => DragModel::G1,
        DragModelArg::G7 => DragModel::G7,
    };

    let zero_inputs = BallisticInputs {
        bc_value: bc,
        bc_type: drag_model_enum,
        bullet_mass: mass_kg,
        muzzle_velocity: velocity_m,
        bullet_diameter: diameter_m,
        bullet_length: diameter_m * 4.5,
        sight_height: sight_height_m,
        use_rk4: true,
        ..Default::default()
    };

    let atmosphere = AtmosphericConditions {
        temperature: temperature_c,
        pressure: pressure_hpa,
        humidity,
        altitude: altitude_m,
        ..Default::default()
    };

    let zero_angle = ballistics_engine::calculate_zero_angle_with_conditions(
        zero_inputs, zero_distance_m, sight_height_m,
        WindConditions::default(), atmosphere,
    )?;

    // Run trajectory with sampling
    let samples = run_sampled_trajectory(
        velocity_m, bc, mass_kg, diameter_m, drag_model, sight_height_m,
        temperature_c, pressure_hpa, humidity, altitude_m,
        wind_speed_m, wind_direction, // degrees, converted internally
        end_m * 1.1, sample_m,
        zero_angle,
    )?;

    // Build output rows at the requested range intervals
    let adj_label = match adjustment_unit {
        AdjustmentUnit::Mil => "MIL",
        AdjustmentUnit::Moa => "MOA",
    };

    let (dist_unit, vel_unit, energy_unit) = match units {
        UnitSystem::Imperial => ("yd", "fps", "ft-lb"),
        UnitSystem::Metric => ("m", "m/s", "J"),
    };

    struct ComeUpRow {
        range: f64,
        drop_adj: f64,
        come_up: f64,
        velocity: f64,
        energy: f64,
        time: f64,
    }

    let mut rows: Vec<ComeUpRow> = Vec::new();
    let mut current_range = start;
    let mut prev_drop_adj: f64 = 0.0;

    while current_range <= end + 0.1 {
        let range_m = UnitConverter::distance_to_metric(current_range, units);

        // Find closest sampled point
        let closest = samples.iter().min_by(|a, b| {
            (a.distance_m - range_m).abs().partial_cmp(&(b.distance_m - range_m).abs()).unwrap()
        });

        if let Some(sample) = closest {
            if (sample.distance_m - range_m).abs() < sample_m * 1.5 {
                let drop_yd = UnitConverter::distance_from_metric(sample.drop_m, units);
                let range_display = UnitConverter::distance_from_metric(sample.distance_m, units);
                let drop_adj = drop_to_adjustment(drop_yd, range_display, adjustment_unit);
                let come_up = drop_adj - prev_drop_adj;

                rows.push(ComeUpRow {
                    range: current_range,
                    drop_adj,
                    come_up,
                    velocity: UnitConverter::velocity_from_metric(sample.velocity_mps, units),
                    energy: UnitConverter::energy_from_metric(sample.energy_j, units),
                    time: sample.time_s,
                });

                prev_drop_adj = drop_adj;
            }
        }

        current_range += step;
    }

    match output {
        OutputFormat::Json => {
            let json_rows: Vec<serde_json::Value> = rows.iter().map(|r| {
                serde_json::json!({
                    "range": r.range,
                    "drop": r.drop_adj,
                    "come_up": r.come_up,
                    "velocity": r.velocity,
                    "energy": r.energy,
                    "time": r.time,
                })
            }).collect();
            let json = serde_json::json!({
                "zero_distance": zero_distance,
                "adjustment_unit": adj_label,
                "distance_unit": dist_unit,
                "velocity_unit": vel_unit,
                "energy_unit": energy_unit,
                "data": json_rows,
            });
            println!("{}", serde_json::to_string_pretty(&json)?);
        }
        OutputFormat::Csv => {
            println!("range_{},drop_{},come_up_{},velocity_{},energy_{},time_s",
                     dist_unit, adj_label.to_lowercase(), adj_label.to_lowercase(), vel_unit, energy_unit);
            for r in &rows {
                println!("{:.0},{:.3},{:.3},{:.0},{:.0},{:.3}",
                         r.range, r.drop_adj, r.come_up, r.velocity, r.energy, r.time);
            }
        }
        OutputFormat::Table | OutputFormat::Pdf => {
            println!();
            println!("Come-Up Table (zero: {:.0} {}, {})", zero_distance, dist_unit, adj_label);
            println!("┌──────────┬──────────┬──────────┬──────────┬──────────┬──────────┐");
            println!("│Range ({:>2})|Drop ({:>3})|Come-Up   │ Vel ({:>3})│Energy    │ Time (s) │",
                     dist_unit, adj_label, vel_unit);
            println!("├──────────┼──────────┼──────────┼──────────┼──────────┼──────────┤");
            for (i, r) in rows.iter().enumerate() {
                let come_up_str = if i == 0 {
                    "    —     ".to_string()
                } else {
                    format!("{:>9.3} ", r.come_up)
                };
                println!("│{:>9.0} │{:>9.3} │{}│{:>9.0} │{:>9.0} │{:>9.3} │",
                         r.range, r.drop_adj, come_up_str, r.velocity, r.energy, r.time);
            }
            println!("└──────────┴──────────┴──────────┴──────────┴──────────┴──────────┘");
        }
    }

    Ok(())
}

/// Wind card handler
fn handle_wind_card(
    velocity: f64, bc: f64, mass: f64, diameter: f64,
    drag_model: DragModelArg,
    zero_distance: f64, wind_speeds: &[f64],
    start: f64, end: f64, step: f64,
    adjustment_unit: AdjustmentUnit,
    sight_height: f64,
    temperature: f64, pressure: f64, humidity: f64, altitude: f64,
    units: UnitSystem,
    output: OutputFormat,
) -> Result<(), Box<dyn Error>> {
    // Convert to metric
    let velocity_m = UnitConverter::velocity_to_metric(velocity, units);
    let mass_kg = UnitConverter::mass_to_metric(mass, units);
    let diameter_m = UnitConverter::diameter_to_metric(diameter, units);
    let sight_height_m = UnitConverter::sight_height_to_metric(sight_height, units);
    let zero_distance_m = UnitConverter::distance_to_metric(zero_distance, units);
    let temperature_c = UnitConverter::temperature_to_metric(temperature, units);
    let pressure_hpa = UnitConverter::pressure_to_metric(pressure, units);
    let altitude_m = UnitConverter::altitude_to_metric(altitude, units);
    let end_m = UnitConverter::distance_to_metric(end, units);
    let sample_m = UnitConverter::distance_to_metric(step, units);

    // Calculate zero angle (no wind)
    let drag_model_enum = match drag_model {
        DragModelArg::G1 => DragModel::G1,
        DragModelArg::G7 => DragModel::G7,
    };

    let zero_inputs = BallisticInputs {
        bc_value: bc,
        bc_type: drag_model_enum,
        bullet_mass: mass_kg,
        muzzle_velocity: velocity_m,
        bullet_diameter: diameter_m,
        bullet_length: diameter_m * 4.5,
        sight_height: sight_height_m,
        use_rk4: true,
        ..Default::default()
    };

    let atmosphere = AtmosphericConditions {
        temperature: temperature_c,
        pressure: pressure_hpa,
        humidity,
        altitude: altitude_m,
        ..Default::default()
    };

    let zero_angle = ballistics_engine::calculate_zero_angle_with_conditions(
        zero_inputs, zero_distance_m, sight_height_m,
        WindConditions::default(), atmosphere,
    )?;

    let adj_label = match adjustment_unit {
        AdjustmentUnit::Mil => "MIL",
        AdjustmentUnit::Moa => "MOA",
    };

    let (dist_unit, wind_unit) = match units {
        UnitSystem::Imperial => ("yd", "mph"),
        UnitSystem::Metric => ("m", "m/s"),
    };

    // For each wind speed, run trajectory with 90° crosswind and collect drift
    // wind_direction = 90° means full value crosswind from the right
    let crosswind_deg = 90.0; // degrees (converted to radians internally)

    // Collect data: rows = ranges, columns = wind speeds
    struct WindRow {
        range: f64,
        drifts: Vec<f64>, // drift in adjustment units per wind speed
    }

    let mut ranges: Vec<f64> = Vec::new();
    let mut current = start;
    while current <= end + 0.1 {
        ranges.push(current);
        current += step;
    }

    let mut all_drifts: Vec<Vec<f64>> = vec![Vec::new(); ranges.len()];

    for &ws in wind_speeds {
        let ws_m = UnitConverter::wind_to_metric(ws, units);

        let samples = run_sampled_trajectory(
            velocity_m, bc, mass_kg, diameter_m, drag_model, sight_height_m,
            temperature_c, pressure_hpa, humidity, altitude_m,
            ws_m, crosswind_deg,
            end_m * 1.1, sample_m,
            zero_angle,
        )?;

        for (ri, &range_display) in ranges.iter().enumerate() {
            let range_m = UnitConverter::distance_to_metric(range_display, units);

            let closest = samples.iter().min_by(|a, b| {
                (a.distance_m - range_m).abs().partial_cmp(&(b.distance_m - range_m).abs()).unwrap()
            });

            let drift_adj = if let Some(sample) = closest {
                if (sample.distance_m - range_m).abs() < sample_m * 1.5 {
                    let drift_yd = UnitConverter::distance_from_metric(sample.wind_drift_m, units);
                    drop_to_adjustment(drift_yd, range_display, adjustment_unit)
                } else {
                    0.0
                }
            } else {
                0.0
            };

            all_drifts[ri].push(drift_adj);
        }
    }

    let wind_rows: Vec<WindRow> = ranges.iter().enumerate().map(|(i, &range)| {
        WindRow { range, drifts: all_drifts[i].clone() }
    }).collect();

    match output {
        OutputFormat::Json => {
            let json_rows: Vec<serde_json::Value> = wind_rows.iter().map(|r| {
                let mut row = serde_json::json!({ "range": r.range });
                for (j, &ws) in wind_speeds.iter().enumerate() {
                    row[format!("wind_{}", ws)] = serde_json::json!(r.drifts.get(j).unwrap_or(&0.0));
                }
                row
            }).collect();
            let json = serde_json::json!({
                "zero_distance": zero_distance,
                "adjustment_unit": adj_label,
                "distance_unit": dist_unit,
                "wind_unit": wind_unit,
                "wind_speeds": wind_speeds,
                "crosswind": "full-value (90°)",
                "data": json_rows,
            });
            println!("{}", serde_json::to_string_pretty(&json)?);
        }
        OutputFormat::Csv => {
            let ws_headers: Vec<String> = wind_speeds.iter().map(|ws| format!("wind_{}_{}", ws, wind_unit)).collect();
            println!("range_{},{}", dist_unit, ws_headers.join(","));
            for r in &wind_rows {
                let drift_strs: Vec<String> = r.drifts.iter().map(|d| format!("{:.1}", d)).collect();
                println!("{:.0},{}", r.range, drift_strs.join(","));
            }
        }
        OutputFormat::Table | OutputFormat::Pdf => {
            println!();
            println!("Wind Card (zero: {:.0} {}, {}, full-value crosswind)", zero_distance, dist_unit, adj_label);

            // Header
            let col_width = 10;
            let range_header = format!("Range ({:>2})", dist_unit);
            let mut header = format!("┌{:─>w$}", "", w = col_width);
            for _ in wind_speeds { header += &format!("┬{:─>w$}", "", w = col_width); }
            header += "┐";
            println!("{}", header);

            let mut label_row = format!("│{:<w$}", range_header, w = col_width);
            for ws in wind_speeds {
                label_row += &format!("│{:>8} {} ", ws, wind_unit);
            }
            label_row += "│";
            println!("{}", label_row);

            let mut sep = format!("├{:─>w$}", "", w = col_width);
            for _ in wind_speeds { sep += &format!("┼{:─>w$}", "", w = col_width); }
            sep += "┤";
            println!("{}", sep);

            for r in &wind_rows {
                let mut row_str = format!("│{:>9.0} ", r.range);
                for d in &r.drifts {
                    row_str += &format!("│{:>9.1} ", d);
                }
                row_str += "│";
                println!("{}", row_str);
            }

            let mut footer = format!("└{:─>w$}", "", w = col_width);
            for _ in wind_speeds { footer += &format!("┴{:─>w$}", "", w = col_width); }
            footer += "┘";
            println!("{}", footer);
        }
    }

    Ok(())
}


// ============================================================================
// Stability Analysis Handler (MBA-734)
// ============================================================================

/// Calculate the Miller stability factor and minimum twist rates
fn handle_stability(
    mass: f64, diameter: f64, length: f64, twist_rate: f64, velocity: f64,
    temperature: f64, pressure: f64, altitude: f64,
    units: UnitSystem,
    output: OutputFormat,
) -> Result<(), Box<dyn Error>> {
    // Convert to metric (SI) for the stability calculation
    let mass_kg = UnitConverter::mass_to_metric(mass, units);
    let diameter_m = UnitConverter::diameter_to_metric(diameter, units);
    let length_m = UnitConverter::sight_height_to_metric(length, units); // inches/mm to meters
    let velocity_m = UnitConverter::velocity_to_metric(velocity, units);
    let temperature_c = UnitConverter::temperature_to_metric(temperature, units);
    let pressure_hpa = UnitConverter::pressure_to_metric(pressure, units);
    let altitude_m = UnitConverter::altitude_to_metric(altitude, units);

    // Twist rate: for imperial it's inches/turn, for metric it's mm/turn
    // compute_stability_coefficient expects twist_rate in inches/turn
    let twist_rate_inches = match units {
        UnitSystem::Imperial => twist_rate,
        UnitSystem::Metric => twist_rate / 25.4, // mm to inches
    };

    // Build BallisticInputs for the stability calculation
    let inputs = BallisticInputs {
        bullet_mass: mass_kg,
        bullet_diameter: diameter_m,
        bullet_length: length_m,
        muzzle_velocity: velocity_m,
        twist_rate: twist_rate_inches,
        ..Default::default()
    };

    let atmo_params = (altitude_m, temperature_c, pressure_hpa, 1.0);

    // Calculate stability factor
    let sg = ballistics_engine::stability::compute_stability_coefficient(&inputs, atmo_params);

    let status = if sg >= 1.5 {
        "Fully Stable"
    } else if sg >= 1.0 {
        "Marginally Stable"
    } else {
        "UNSTABLE"
    };

    // Calculate minimum twist rate for Sg=1.5 and Sg=1.0 by binary search.
    // Lower twist number = faster spin = higher Sg.
    let find_twist_for_sg = |target_sg: f64| -> f64 {
        let mut low = 1.0_f64;   // 1 inch/turn (extremely fast)
        let mut high = 40.0_f64; // 40 inches/turn (very slow)

        for _ in 0..100 {
            let mid = (low + high) / 2.0;
            let test_inputs = BallisticInputs {
                bullet_mass: mass_kg,
                bullet_diameter: diameter_m,
                bullet_length: length_m,
                muzzle_velocity: velocity_m,
                twist_rate: mid,
                ..Default::default()
            };
            let test_sg = ballistics_engine::stability::compute_stability_coefficient(
                &test_inputs, atmo_params,
            );
            if test_sg > target_sg {
                low = mid;
            } else {
                high = mid;
            }
            if (test_sg - target_sg).abs() < 0.001 {
                break;
            }
        }
        (low + high) / 2.0
    };

    let min_twist_1_5 = find_twist_for_sg(1.5);
    let min_twist_1_0 = find_twist_for_sg(1.0);

    // min_twist_* are in INCHES (compute_stability_coefficient treats twist_rate as inches/turn);
    // convert to mm for metric so the JSON/CSV values match their labeled unit (the Table path
    // already converts). Imperial is unchanged.
    let (min15_out, min10_out) = match units {
        UnitSystem::Imperial => (min_twist_1_5, min_twist_1_0),
        UnitSystem::Metric => (min_twist_1_5 * 25.4, min_twist_1_0 * 25.4),
    };

    // Display units
    let (len_unit, vel_unit, twist_display, min15_display, min10_display) = match units {
        UnitSystem::Imperial => (
            "\"",
            "fps",
            format!("1:{:.0}\" RH", twist_rate),
            format!("1:{:.1}\"", min_twist_1_5),
            format!("1:{:.1}\"", min_twist_1_0),
        ),
        UnitSystem::Metric => (
            "mm",
            "m/s",
            format!("1:{:.0}mm RH", twist_rate),
            format!("1:{:.1}mm", min_twist_1_5 * 25.4),
            format!("1:{:.1}mm", min_twist_1_0 * 25.4),
        ),
    };

    match output {
        OutputFormat::Json => {
            let json = serde_json::json!({
                "stability_factor": (sg * 100.0).round() / 100.0,
                "status": status,
                "twist_rate": twist_rate,
                "twist_rate_unit": if units == UnitSystem::Imperial { "in/turn" } else { "mm/turn" },
                "min_twist_sg_1_5": (min15_out * 100.0).round() / 100.0,
                "min_twist_sg_1_0": (min10_out * 100.0).round() / 100.0,
                "bullet_length": length,
                "bullet_diameter": diameter,
                "bullet_mass": mass,
                "muzzle_velocity": velocity,
                "units": if units == UnitSystem::Imperial { "imperial" } else { "metric" },
            });
            println!("{}", serde_json::to_string_pretty(&json)?);
        }
        OutputFormat::Csv => {
            println!("sg,status,twist_rate,min_twist_sg1.5,min_twist_sg1.0,length,velocity");
            println!("{:.2},{},{:.1},{:.1},{:.1},{:.3},{:.0}",
                     sg, status, twist_rate, min15_out, min10_out, length, velocity);
        }
        OutputFormat::Table | OutputFormat::Pdf => {
            println!();
            println!("╔════════════════════════════════════════╗");
            println!("║         Stability Analysis            ║");
            println!("╠════════════════════════════════════════╣");
            println!("║ Stability Factor (Sg):  {:>12.2}  ║", sg);
            println!("║ Status:              {:>14}  ║", status);
            println!("║ Twist Rate:          {:>14}  ║", twist_display);
            println!("║ Min Twist (Sg=1.5):  {:>14}  ║", min15_display);
            println!("║ Min Twist (Sg=1.0):  {:>14}  ║", min10_display);
            println!("║ Bullet Length:     {:>12.3}{:>4}  ║", length, len_unit);
            println!("║ Muzzle Velocity:     {:>10.0} {:<3}  ║", velocity, vel_unit);
            println!("╚════════════════════════════════════════╝");
        }
    }

    Ok(())
}

// ============================================================================
// Range Table Handler (MBA-733)
// ============================================================================

/// Generate a comprehensive range table combining drop, wind, velocity, energy, and ToF
fn handle_range_table(
    velocity: f64, bc: f64, mass: f64, diameter: f64,
    drag_model: DragModelArg,
    zero_distance: f64, start: f64, end: f64, step: f64,
    wind_speed: f64, wind_direction: f64,
    adjustment_unit: AdjustmentUnit,
    sight_height: f64,
    temperature: f64, pressure: f64, humidity: f64, altitude: f64,
    units: UnitSystem,
    output: OutputFormat,
) -> Result<(), Box<dyn Error>> {
    // Convert to metric
    let velocity_m = UnitConverter::velocity_to_metric(velocity, units);
    let mass_kg = UnitConverter::mass_to_metric(mass, units);
    let diameter_m = UnitConverter::diameter_to_metric(diameter, units);
    let sight_height_m = UnitConverter::sight_height_to_metric(sight_height, units);
    let zero_distance_m = UnitConverter::distance_to_metric(zero_distance, units);
    let temperature_c = UnitConverter::temperature_to_metric(temperature, units);
    let pressure_hpa = UnitConverter::pressure_to_metric(pressure, units);
    let altitude_m = UnitConverter::altitude_to_metric(altitude, units);
    let wind_speed_m = UnitConverter::wind_to_metric(wind_speed, units);
    let end_m = UnitConverter::distance_to_metric(end, units);
    let sample_m = UnitConverter::distance_to_metric(step, units);

    // Calculate zero angle (no wind for clean zero)
    let drag_model_enum = match drag_model {
        DragModelArg::G1 => DragModel::G1,
        DragModelArg::G7 => DragModel::G7,
    };

    let zero_inputs = BallisticInputs {
        bc_value: bc,
        bc_type: drag_model_enum,
        bullet_mass: mass_kg,
        muzzle_velocity: velocity_m,
        bullet_diameter: diameter_m,
        bullet_length: diameter_m * 4.5,
        sight_height: sight_height_m,
        use_rk4: true,
        ..Default::default()
    };

    let atmosphere = AtmosphericConditions {
        temperature: temperature_c,
        pressure: pressure_hpa,
        humidity,
        altitude: altitude_m,
        ..Default::default()
    };

    let zero_angle = ballistics_engine::calculate_zero_angle_with_conditions(
        zero_inputs, zero_distance_m, sight_height_m,
        WindConditions::default(), atmosphere,
    )?;

    // Run trajectory WITH wind (for wind drift values)
    let wind_samples = run_sampled_trajectory(
        velocity_m, bc, mass_kg, diameter_m, drag_model, sight_height_m,
        temperature_c, pressure_hpa, humidity, altitude_m,
        wind_speed_m, wind_direction,
        end_m * 1.1, sample_m,
        zero_angle,
    )?;

    // Run trajectory WITHOUT wind (for pure drop)
    let no_wind_samples = run_sampled_trajectory(
        velocity_m, bc, mass_kg, diameter_m, drag_model, sight_height_m,
        temperature_c, pressure_hpa, humidity, altitude_m,
        0.0, 0.0,
        end_m * 1.1, sample_m,
        zero_angle,
    )?;

    // Build output rows
    let adj_label = match adjustment_unit {
        AdjustmentUnit::Mil => "MIL",
        AdjustmentUnit::Moa => "MOA",
    };

    let (dist_unit, vel_unit, energy_unit, drop_unit, wind_unit_label) = match units {
        UnitSystem::Imperial => ("yd", "fps", "ft-lb", "in", "mph"),
        UnitSystem::Metric => ("m", "m/s", "J", "mm", "m/s"),
    };

    struct RangeRow {
        range: f64,
        drop_linear: f64,
        drop_adj: f64,
        wind_linear: f64,
        wind_adj: f64,
        velocity: f64,
        energy: f64,
        time: f64,
    }

    let mut rows: Vec<RangeRow> = Vec::new();
    let mut current_range = start;

    while current_range <= end + 0.1 {
        let range_m = UnitConverter::distance_to_metric(current_range, units);

        let nw_closest = no_wind_samples.iter().min_by(|a, b| {
            (a.distance_m - range_m).abs().partial_cmp(&(b.distance_m - range_m).abs()).unwrap()
        });

        let w_closest = wind_samples.iter().min_by(|a, b| {
            (a.distance_m - range_m).abs().partial_cmp(&(b.distance_m - range_m).abs()).unwrap()
        });

        if let (Some(nw), Some(w)) = (nw_closest, w_closest) {
            if (nw.distance_m - range_m).abs() < sample_m * 1.5 {
                let range_display = UnitConverter::distance_from_metric(nw.distance_m, units);

                let drop_linear = match units {
                    UnitSystem::Imperial => nw.drop_m / 0.0254,
                    UnitSystem::Metric => nw.drop_m * 1000.0,
                };

                let drop_yd = UnitConverter::distance_from_metric(nw.drop_m, units);
                let drop_adj = drop_to_adjustment(drop_yd, range_display, adjustment_unit);

                let wind_linear = match units {
                    UnitSystem::Imperial => w.wind_drift_m / 0.0254,
                    UnitSystem::Metric => w.wind_drift_m * 1000.0,
                };

                let drift_yd = UnitConverter::distance_from_metric(w.wind_drift_m, units);
                let wind_adj = drop_to_adjustment(drift_yd, range_display, adjustment_unit);

                rows.push(RangeRow {
                    range: current_range,
                    drop_linear,
                    drop_adj,
                    wind_linear,
                    wind_adj,
                    velocity: UnitConverter::velocity_from_metric(nw.velocity_mps, units),
                    energy: UnitConverter::energy_from_metric(nw.energy_j, units),
                    time: nw.time_s,
                });
            }
        }

        current_range += step;
    }

    match output {
        OutputFormat::Json => {
            let json_rows: Vec<serde_json::Value> = rows.iter().map(|r| {
                serde_json::json!({
                    "range": r.range,
                    "drop": r.drop_linear,
                    "drop_adj": r.drop_adj,
                    "wind_drift": r.wind_linear,
                    "wind_adj": r.wind_adj,
                    "velocity": r.velocity,
                    "energy": r.energy,
                    "time": r.time,
                })
            }).collect();
            let json = serde_json::json!({
                "zero_distance": zero_distance,
                "wind_speed": wind_speed,
                "wind_direction": wind_direction,
                "adjustment_unit": adj_label,
                "distance_unit": dist_unit,
                "velocity_unit": vel_unit,
                "energy_unit": energy_unit,
                "drop_unit": drop_unit,
                "data": json_rows,
            });
            println!("{}", serde_json::to_string_pretty(&json)?);
        }
        OutputFormat::Csv => {
            println!("range_{},drop_{},drop_{},wind_{},wind_{},velocity_{},energy_{},tof_s",
                     dist_unit, drop_unit, adj_label.to_lowercase(),
                     drop_unit, adj_label.to_lowercase(),
                     vel_unit, energy_unit);
            for r in &rows {
                println!("{:.0},{:.1},{:.3},{:.1},{:.2},{:.0},{:.0},{:.2}",
                         r.range, r.drop_linear, r.drop_adj,
                         r.wind_linear, r.wind_adj,
                         r.velocity, r.energy, r.time);
            }
        }
        OutputFormat::Table | OutputFormat::Pdf => {
            println!();
            println!("Range Table (zero: {:.0} {}, {:.0} {} crosswind)", zero_distance, dist_unit, wind_speed, wind_unit_label);
            println!("┌───────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬───────┐");
            println!("│Range  │Drop     │Drop     │Wind     │Wind     │Vel      │Energy   │ToF    │");
            println!("│({:>2})   │({:>2})     │({:>3})    │({:>2})     │({:>3})    │({:>3})    │({:>4})   │(s)    │",
                     dist_unit, drop_unit, adj_label, drop_unit, adj_label, vel_unit, energy_unit);
            println!("├───────┼─────────┼─────────┼─────────┼─────────┼─────────┼─────────┼───────┤");
            for r in &rows {
                println!("│{:>6.0} │{:>8.1} │{:>8.3} │{:>8.1} │{:>8.2} │{:>8.0} │{:>8.0} │{:>6.2} │",
                         r.range, r.drop_linear, r.drop_adj,
                         r.wind_linear, r.wind_adj,
                         r.velocity, r.energy, r.time);
            }
            println!("└───────┴─────────┴─────────┴─────────┴─────────┴─────────┴─────────┴───────┘");
        }
    }

    Ok(())
}
