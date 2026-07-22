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

mod mcp_command;
#[cfg(feature = "pdf")]
mod pdf_dope_card;
mod solve_json_command;
// MBA-1355: parse_click_value/clicks_for/ClickValue power the CLI's click-graduation
// flags (--elevation-click-value/--windage-click-value) and profile fields.
use ballistics_engine::adjustment::{adjustment_factor, clicks_for, parse_click_value, ClickBase, ClickValue};
use ballistics_engine::terminal_plot;
#[cfg(feature = "pdf")]
use pdf_dope_card::{calculate_density_altitude, DopeCardConfig, DopeCardRow, FontSizePreset};

#[cfg(feature = "online")]
use ballistics_engine::api_client::{ApiClient, TrajectoryRequestBuilder, TrueVelocityRequest};
use ballistics_engine::bc_table::BcCorrectionTable;
use ballistics_engine::bc_table_5d::Bc5dTableManager;
#[cfg(feature = "online")]
use ballistics_engine::bc_table_download::Bc5dDownloader;
use ballistics_engine::constants::{GRAINS_TO_KG, DEFAULT_POWDER_REFERENCE_TEMP_C, DEFAULT_POWDER_REFERENCE_TEMP_F, GRAMS_PER_GRAIN};
use ballistics_engine::cli_api::UnitSystem;
use ballistics_engine::truing::{
    calculate_true_velocity_local, fallback_bullet_length_m, parse_truing_observation,
    run_multi_observation_truing_core, validate_truing_observations, DragModelArg, DropUnit,
    MultiTruingReport, TruingModelInputsV1, TruingObservation,
};
use ballistics_engine::truing_plan::{
    plan_truing_experiment_v1, TruingExperimentPlanRequestV1, TruingExperimentPlanV1,
    TruingPlanModeV1,
};
use ballistics_engine::truing_uncertainty::{
    run_uncertainty_truing_v1, NormalPriorV1, TruingApproximationV1, TruingPredictionRequestV1, TruingPriorsV1,
    UncertaintyTruingReportV1, UncertaintyTruingRequestV1, WeightedTruingObservationV1,
};
use ballistics_engine::wez::{compute_wez, parse_target_size, TargetSizeMetric, WezResult, WezRow};
use ballistics_engine::{
    trajectory_sampling, AtmosphericConditions, BCSegmentData, BallisticInputs, DragModel,
    MonteCarloParams, TrajectorySolver, WindConditions,
};
use clap::{CommandFactory, Parser, Subcommand, ValueEnum};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::error::Error;
use std::fs;
use std::io::{self, Write};
use std::path::PathBuf;
use strsim::levenshtein;

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
    let mut response = ureq::get(TOS_URL)
        .config()
        .timeout_global(Some(std::time::Duration::from_secs(10)))
        .build()
        .call()
        .map_err(|e| format!("Failed to fetch Terms of Service: {}", e))?;

    response
        .body_mut()
        .read_to_string()
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
        eprintln!(
            "Terms of Service accepted. This will not be shown again unless the terms change."
        );
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

#[allow(
    clippy::large_enum_variant,
    reason = "variants intentionally mirror the stable flat Clap command shape"
)]
#[derive(Subcommand)]
enum Commands {
    /// Solve one explicit-SI v1 JSON request from stdin and write one JSON envelope to stdout
    SolveJson,

    /// Run a Model Context Protocol (MCP) server over stdio, exposing solve and engine_info tools
    Mcp,

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

        /// Vertical wind, mph (imperial) or m/s (metric); positive = updraft (raises POI)
        #[arg(
            long,
            value_parser = f64_range(-100.0, 100.0),
            default_value_t = 0.0,
            allow_hyphen_values = true
        )]
        wind_vertical: f64,

        /// Downrange wind segment "SPEED:ANGLE:UNTIL_DISTANCE" with an optional 4th
        /// ":VERTICAL" field (repeatable).
        /// SPEED/UNTIL_DISTANCE follow --units (mph & yd imperial, m/s & m metric); ANGLE is
        /// degrees, same convention as --wind-direction. The optional 4th field, VERTICAL, is
        /// ALWAYS m/s (positive = updraft, raises POI) regardless of --units — unlike SPEED it
        /// does not follow --units; omit it for the pre-existing 3-field form (vertical wind
        /// 0). Each segment applies from the previous boundary out to UNTIL_DISTANCE; wind is
        /// zero beyond the last segment. Overrides --wind-speed/-direction/-vertical. Not
        /// compatible with --enable-wind-shear.
        #[arg(long = "wind-segment", value_name = "SPEED:ANGLE:DIST[:VERTICAL]", action = clap::ArgAction::Append)]
        wind_segment: Vec<String>,

        /// Manual velocity-dependent BC segment(s): "VMIN:VMAX:BC" (repeatable). VMIN/VMAX are
        /// velocities in --units (fps imperial, m/s metric); the given BC applies while the
        /// bullet's speed is in [VMIN, VMAX). Segments are keyed to VELOCITY (independent of
        /// distance/wind). Overrides --bc-table and --bc-table-dir, and implies
        /// --use-bc-segments.
        #[arg(long = "bc-segment", value_name = "VMIN:VMAX:BC", action = clap::ArgAction::Append)]
        bc_segment: Vec<String>,

        /// Temperature (Fahrenheit or Celsius based on --units; default 59 F / 15 C)
        #[arg(long)]
        temperature: Option<f64>,

        /// Pressure (inHg or hPa based on --units; default 29.92 inHg / 1013.25 hPa)
        #[arg(long)]
        pressure: Option<f64>,

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

        /// Render an inline terminal chart (drop vs. range, then lateral drift vs. range)
        /// after the normal output. Only applies to `-o table` (the default); other
        /// `--output` formats stay machine-readable and unchanged. Bare `--plot` uses the
        /// Unicode braille-dot canvas; `--plot ascii` uses a '*'-per-cell fallback for
        /// terminals/fonts without braille glyph coverage. Monochrome (no ANSI colors —
        /// see src/terminal_plot.rs for why); fixed 72x12-cell canvas per chart, no
        /// terminal-size detection.
        #[arg(long, value_enum, num_args = 0..=1, default_missing_value = "braille")]
        plot: Option<PlotStyle>,

        /// Automatically zero to target distance (overrides --angle)
        #[arg(long)]
        auto_zero: Option<f64>,

        /// Sight height above bore (inches for imperial, mm for metric)
        #[arg(long)]
        sight_height: Option<f64>,

        /// Bore height above ground (inches for imperial, mm for metric; MBA-1339 unified
        /// this with --sight-height and the WASM --muzzle-height flag). Default: 60in/1500mm
        /// (= 5ft/1.5m). Also accepts --muzzle-height.
        #[arg(long, visible_alias = "muzzle-height")]
        bore_height: Option<f64>,

        /// Ignore ground impact detection (trajectory continues to max range)
        #[arg(long)]
        ignore_ground_impact: bool,

        /// Enable velocity-based BC segmentation
        #[arg(long)]
        use_bc_segments: bool,

        /// Path to a custom drag deck: CSV with `mach,cd` per line (Hornady CDM / Lapua radar
        /// style). Replaces the G-model + BC for drag; bc_value is ignored when set. Mach-keyed;
        /// out-of-range Mach holds the nearest tabulated Cd.
        #[arg(long, value_name = "FILE")]
        drag_table: Option<PathBuf>,

        /// Rifle cant in DEGREES, positive = clockwise from the shooter (top of scope tips
        /// right). Models "zeroed level, fired canted": POI moves right and low. 0 = level.
        #[arg(long, alias = "cant-angle", value_name = "DEGREES", default_value_t = 0.0)]
        cant: f64,

        /// Print the BC5D-generated segment ladder as ready-to-paste
        /// --bc-segment arguments (requires --bc-table-dir; velocities in
        /// the active --units)
        #[arg(long)]
        print_bc_segments: bool,

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

        /// Enable empirical Litz spin drift calculations
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
        #[arg(
            long,
            value_name = "FILE",
            help = "Use precomputed BC correction table instead of online API"
        )]
        bc_table: Option<PathBuf>,

        /// Directory containing caliber-specific BC5D tables (e.g., bc5d_308.bin)
        #[arg(
            long,
            value_name = "DIR",
            help = "Use caliber-specific 5D BC correction tables"
        )]
        bc_table_dir: Option<PathBuf>,

        /// Enable auto-download of BC5D correction tables when needed
        #[cfg(feature = "online")]
        #[arg(long, help = "Auto-download BC5D tables when needed")]
        bc_table_auto: bool,

        /// Base URL for BC5D table downloads
        #[cfg(feature = "online")]
        #[arg(
            long,
            default_value = "https://ballistics.tools/downloads/bc5d",
            help = "Base URL for BC5D table downloads"
        )]
        bc_table_url: String,

        /// Force re-download of BC5D tables even if cached
        #[cfg(feature = "online")]
        #[arg(long, help = "Force re-download even if cached")]
        bc_table_refresh: bool,

        /// Bullet length in inches (for BC table lookup; estimated from diameter if not provided)
        #[arg(long)]
        bullet_length: Option<f64>,

        /// Barrel twist rate (inches/turn for imperial, mm/turn for metric; e.g. imperial 10 = 1:10")
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

        /// Powder temperature sensitivity (fps/°F imperial, m/s/°C metric).
        /// Defaults to 1.0 fps/°F (0.54864 m/s/°C).
        #[arg(long)]
        powder_temp_sensitivity: Option<f64>,

        /// Powder temperature (°F/°C). With --powder-temp-curve, this is the powder
        /// temperature the curve is interpolated at to resolve muzzle velocity (defaults
        /// to --temperature when omitted, i.e. powder assumed at air temperature). With
        /// --powder-temp-sensitivity (linear model), it is the reference temperature the
        /// stated velocity was measured at (defaults to the 70°F equivalent).
        #[arg(long)]
        powder_temp: Option<f64>,

        /// Measured powder-temperature -> muzzle-velocity curve as comma-separated
        /// TEMP:VELOCITY points, e.g. "40:2620,70:2700,100:2760" (temps in F/C,
        /// velocities in fps / m·s⁻¹ per --units). The muzzle velocity is interpolated
        /// from this table at the powder temperature (--powder-temp, or --temperature if
        /// unset; clamped at the endpoints, no extrapolation). Data-driven, non-linear
        /// alternative that OVERRIDES --powder-temp-sensitivity when supplied.
        #[arg(long = "powder-temp-curve", value_name = "TEMP:VEL,...")]
        powder_temp_curve: Option<String>,

        // Zero-day condition overrides for --auto-zero. When supplied, the zero ANGLE is
        // solved under these conditions (the day/velocity the rifle was actually zeroed in)
        // while the trajectory runs under the current shot-day conditions. With no zero-day
        // flags, the shot-day values are reused exactly; coupled powder fallbacks are documented
        // on the individual options below.
        /// Zero-day muzzle velocity for --auto-zero (fps imperial / m·s⁻¹ metric). Use when the
        /// rifle was zeroed at a different velocity than this shot. Overrides both powder models.
        #[arg(long, value_parser = f64_range(0.0, 6000.0))]
        zero_velocity: Option<f64>,

        /// Zero-day air temperature for --auto-zero (°F imperial / °C metric). Also resolves
        /// zero-day velocity for the linear powder-sensitivity model unless --zero-velocity is set.
        #[arg(long, allow_hyphen_values = true)]
        zero_temperature: Option<f64>,

        /// Zero-day barometric pressure for --auto-zero (inHg imperial / hPa metric)
        #[arg(long)]
        zero_pressure: Option<f64>,

        /// Zero-day relative humidity for --auto-zero (percent, 0-100)
        #[arg(long, value_parser = f64_range(0.0, 100.0))]
        zero_humidity: Option<f64>,

        /// Zero-day altitude for --auto-zero (feet imperial / meters metric)
        #[arg(long, allow_hyphen_values = true)]
        zero_altitude: Option<f64>,

        /// Zero-day powder temperature for --auto-zero (°F/°C). With --powder-temp-curve,
        /// the curve is interpolated at this temperature to resolve the zero-day muzzle
        /// velocity. When omitted, the curve follows an explicit --zero-temperature; otherwise
        /// it inherits an explicit shot-day --powder-temp. An explicit --zero-velocity still
        /// takes precedence.
        #[arg(long, allow_hyphen_values = true)]
        zero_powder_temp: Option<f64>,

        // Online Mode Parameters (feature-gated)
        /// Use Flask API for ML-enhanced trajectory calculation
        #[cfg(feature = "online")]
        #[arg(
            long,
            help = "Route calculations through Flask API for ML enhancements"
        )]
        online: bool,

        /// API endpoint URL
        #[cfg(feature = "online")]
        #[arg(
            long,
            default_value = "https://api.ballistics.7.62x51mm.sh",
            help = "API endpoint URL"
        )]
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

        /// Wind shear boundary-layer model (only used with --enable-wind-shear)
        #[arg(
            long,
            value_enum,
            default_value = "power_law",
            help = "Wind shear model: none, logarithmic, power_law, ekman_spiral, custom_layers"
        )]
        wind_shear_model: WindShearModelArg,

        /// Weather zone interpolation method
        #[cfg(feature = "online")]
        #[arg(
            long,
            default_value = "linear",
            help = "Weather zone interpolation: linear, cubic, step"
        )]
        weather_zone_interpolation: String,

        // PDF Dope Card Parameters
        /// Target speed for moving-target holds (mph for imperial, m/s for metric — same
        /// convention as `lead --target-speed`). Drives two independent outputs when > 0:
        /// the PDF dope card's assumed-90-degree-crossing Lead column, and — for every
        /// output format — a per-point mover Ring column/field (MBA-1325): ring radius =
        /// target_speed x point ToF, the "fire when the mover enters the ring" technique.
        /// See CLI_USAGE.md for the full mover-ring writeup. 0 (default) disables both.
        #[arg(long, default_value = "0.0", value_parser = f64_range(0.0, 300.0))]
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

        /// Angular unit for the PDF dope card's Drop/Wind/Lead columns AND the
        /// --target-speed mover Ring table column (mil, moa, smoa, iphy, or clicks).
        /// CSV/JSON ring fields stay mil/meters — their names carry the unit contract.
        /// `clicks` requires --elevation-click-value (or a saved profile's
        /// elevation_click) and formats whole turret clicks instead of an angle.
        #[arg(long, value_enum, default_value = "mil")]
        adjustment_unit: AdjustmentUnit,

        /// Turret elevation click graduation for `--adjustment-unit clicks`, e.g. 0.1mil or
        /// 0.25moa (falls back to the saved profile's elevation_click when omitted).
        #[arg(long)]
        elevation_click_value: Option<String>,

        /// Turret windage click graduation for `--adjustment-unit clicks` (falls back to the
        /// resolved elevation graduation, then the saved profile's windage_click).
        #[arg(long)]
        windage_click_value: Option<String>,
    },

    /// Run Monte Carlo simulation
    MonteCarlo {
        /// Base velocity (fps for imperial, m/s for metric)
        #[arg(short = 'v', long, value_parser = f64_range(0.0, 6000.0))]
        velocity: f64,

        /// Launch angle (degrees)
        #[arg(short = 'a', long, default_value = "0.0")]
        angle: f64,

        /// Ballistic coefficient
        #[arg(short = 'b', long, value_parser = f64_range(0.001, 2.0))]
        bc: f64,

        /// Mass (grains for imperial, grams for metric)
        #[arg(short = 'm', long, value_parser = f64_range(0.1, 2000.0))]
        mass: f64,

        /// Diameter (inches for imperial, mm for metric)
        #[arg(short = 'd', long, value_parser = f64_range(0.01, 60.0))]
        diameter: f64,

        /// Number of simulations
        #[arg(short = 'n', long, default_value = "1000")]
        num_sims: usize,

        /// Velocity standard deviation (fps for imperial, m/s for metric)
        #[arg(long, default_value = "1.0")]
        velocity_std: f64,

        /// Angle standard deviation (degrees)
        #[arg(long, default_value = "0.1")]
        angle_std: f64,

        /// BC standard deviation
        #[arg(long, default_value = "0.01")]
        bc_std: f64,

        /// Wind speed standard deviation (mph for imperial, m/s for metric)
        #[arg(long, default_value = "1.0")]
        wind_std: f64,

        /// Wind direction standard deviation (degrees)
        #[arg(long, visible_alias = "wind-dir-std", default_value = "0.0")]
        wind_direction_std: f64,

        /// Base wind speed (mph for imperial, m/s for metric)
        #[arg(long, default_value = "0.0")]
        wind_speed: f64,

        /// Base wind direction (degrees; wind-FROM: 0=headwind, 90=from right, 180=tailwind)
        #[arg(long, default_value = "0.0")]
        wind_direction: f64,

        /// Base vertical wind, mph (imperial) or m/s (metric); positive = updraft (raises POI)
        #[arg(
            long,
            value_parser = f64_range(-100.0, 100.0),
            default_value_t = 0.0,
            allow_hyphen_values = true
        )]
        wind_vertical: f64,

        /// Target distance (yards for imperial, meters for metric)
        #[arg(long)]
        target_distance: Option<f64>,

        /// Hit-zone radius around the point of aim at the target plane (yards for imperial,
        /// meters for metric). A "hit" is a simulation landing within this radius. Default 0.3 m.
        #[arg(long)]
        target_radius: Option<f64>,

        /// Path to a custom drag deck: CSV `mach,cd` per line. Replaces the G-model + BC for drag;
        /// bc_value is ignored when set. Mach-keyed; out-of-range Mach holds the nearest tabulated Cd.
        #[arg(long, value_name = "FILE")]
        drag_table: Option<std::path::PathBuf>,

        /// Rifle cant in DEGREES, positive = clockwise from the shooter (top of scope tips
        /// right). Models "zeroed level, fired canted": POI moves right and low. 0 = level.
        #[arg(long, alias = "cant-angle", value_name = "DEGREES", default_value_t = 0.0)]
        cant: f64,

        /// WEZ (Weapon Employment Zone) sweep mode: report hit probability vs range for a
        /// target size instead of a single summary. See CLI_USAGE.md's WEZ section.
        #[arg(long)]
        wez: bool,

        /// WEZ target size: `WIDTHxHEIGHT` (inches for imperial, cm for metric; e.g. `18x30`)
        /// for a rectangle centered on the point of aim, or a single number for a circular
        /// radius (same units; same semantics as `--target-radius`, just in target-size units).
        /// Required with `--wez`.
        #[arg(long, value_name = "WxH|R", requires = "wez")]
        target_size: Option<String>,

        /// Shooter's wind CALL error (mph for imperial, m/s for metric): the uncertainty in the
        /// shooter's own estimate of wind speed, as distinct from `--wind-std`'s physical
        /// gustiness. Composed with `--wind-std` in quadrature for the underlying dispersion:
        /// effective_wind_std = sqrt(wind_std^2 + wind_call_error^2). Only used with `--wez`.
        #[arg(long, default_value = "0.0", requires = "wez")]
        wind_call_error: f64,

        /// WEZ sweep start range (yards for imperial, meters for metric).
        #[arg(long, default_value = "200.0", requires = "wez")]
        wez_start: f64,

        /// WEZ sweep end range, inclusive (yards for imperial, meters for metric).
        #[arg(long, default_value = "1000.0", requires = "wez")]
        wez_end: f64,

        /// WEZ sweep step (yards for imperial, meters for metric).
        #[arg(long, default_value = "100.0", requires = "wez")]
        wez_step: f64,

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

        /// Temperature (Fahrenheit or Celsius based on --units; default 59 F / 15 C)
        #[arg(long)]
        temperature: Option<f64>,

        /// Pressure (inHg or hPa based on --units; default 29.92 inHg / 1013.25 hPa)
        #[arg(long)]
        pressure: Option<f64>,

        /// Humidity (0-100%)
        #[arg(long, default_value = "50.0", value_parser = f64_range(0.0, 100.0))]
        humidity: f64,

        /// Altitude (feet or meters based on --units)
        #[arg(long, default_value = "0.0")]
        altitude: f64,

        /// Path to a custom drag deck: CSV `mach,cd` per line. Replaces the G-model + BC for drag;
        /// bc_value is ignored when set. Mach-keyed; out-of-range Mach holds the nearest tabulated Cd.
        #[arg(long, value_name = "FILE")]
        drag_table: Option<std::path::PathBuf>,

        /// Output format
        #[arg(short = 'o', long, default_value = "table")]
        output: OutputFormat,
    },

    /// Estimate BC from trajectory data (drop and/or velocity), for G1, G7, or both
    EstimateBC {
        /// Initial velocity (fps for imperial, m/s for metric)
        #[arg(short = 'v', long, value_parser = f64_range(0.0, 6000.0))]
        velocity: f64,

        /// Mass (grains for imperial, grams for metric)
        #[arg(short = 'm', long, value_parser = f64_range(0.1, 2000.0))]
        mass: f64,

        /// Diameter (inches for imperial, mm for metric)
        #[arg(short = 'd', long, value_parser = f64_range(0.01, 60.0))]
        diameter: f64,

        /// Drop data: "dist,drop;dist,drop;..." (distance in yd/m, drop in in/m).
        /// n-point alternative to --distance1/--drop1/--distance2/--drop2.
        #[arg(long)]
        data: Option<String>,

        /// Velocity data: "dist,vel;dist,vel;..." (distance in yd/m, velocity in fps/mps).
        /// Enables the velocity-retention fit variants.
        #[arg(long)]
        velocity_data: Option<String>,

        /// Drag model to estimate: g1, g7, or both
        #[arg(long, default_value = "both")]
        drag_model: String,

        /// Zero range of the drop data, in yd/m (dope cards are zeroed). Given → drop is
        /// fit as drop-below-line-of-sight from a rifle zeroed here. Omitted → drop is
        /// treated as bore-referenced (flat-fire drop below the extended bore axis).
        #[arg(long)]
        zero_range: Option<f64>,

        /// Sight height above bore (inches/mm) for the zeroed drop fit [default 2 in / 50 mm]
        #[arg(long)]
        sight_height: Option<f64>,

        /// Air temperature the data was measured at (°F imperial / °C metric) [default 59 F / 15 C]
        #[arg(long)]
        temperature: Option<f64>,

        /// Barometric pressure the data was measured at (inHg imperial / hPa metric) [default 29.92 / 1013.25]
        #[arg(long)]
        pressure: Option<f64>,

        /// Relative humidity the data was measured at (percent, 0–100)
        #[arg(long, default_value = "50.0", value_parser = f64_range(0.0, 100.0))]
        humidity: f64,

        /// Altitude the data was measured at (feet imperial / meters metric)
        #[arg(long, default_value = "0.0")]
        altitude: f64,

        /// Distance 1 (yd/m) — legacy 2-point drop input
        #[arg(long)]
        distance1: Option<f64>,

        /// Drop at distance 1 (inches for imperial, meters for metric) — legacy
        #[arg(long, allow_hyphen_values = true)]
        drop1: Option<f64>,

        /// Distance 2 (yd/m) — legacy 2-point drop input
        #[arg(long)]
        distance2: Option<f64>,

        /// Drop at distance 2 (inches for imperial, meters for metric) — legacy
        #[arg(long, allow_hyphen_values = true)]
        drop2: Option<f64>,

        /// Output format
        #[arg(short = 'o', long, default_value = "table")]
        output: OutputFormat,
    },

    /// Generate BC segments for velocity-dependent BC
    GenerateBCSegments {
        /// Base ballistic coefficient
        #[arg(short = 'b', long, value_parser = f64_range(0.001, 2.0))]
        bc: f64,

        /// Projectile mass (grains for imperial, grams for metric)
        #[arg(short = 'm', long, value_parser = f64_range(0.1, 2000.0))]
        mass: f64,

        /// Projectile diameter (inches for imperial, mm for metric)
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
        /// Measured drop at the target range. Unit is MIL by default; in
        /// multi-observation mode (see --observed) it follows --drop-unit and
        /// counts as the first observation.
        #[arg(long)]
        measured_drop: f64,

        /// Range at which drop was measured (yards for imperial, meters for metric)
        #[arg(long)]
        range: f64,

        /// Additional observed impact "RANGE:DROP" for joint MV+BC calibration
        /// (MBA-1316). Repeatable. RANGE follows --units (yd/m); DROP follows
        /// --drop-unit. With >=1 --observed the command fits muzzle velocity and
        /// (when the observations constrain it) ballistic coefficient jointly via
        /// the real forward model. With zero --observed it behaves exactly as the
        /// classic single-observation velocity truing.
        #[arg(long = "observed", value_name = "RANGE:DROP", action = clap::ArgAction::Append)]
        observed: Vec<String>,

        /// Known one-standard-deviation error for an impact reading, in
        /// --drop-unit. Supplying this enables uncertainty-aware joint MV/BC
        /// truing; an optional third field in --observed (RANGE:DROP:SIGMA)
        /// overrides it for that observation.
        #[arg(long, value_name = "SIGMA")]
        observation_sigma: Option<f64>,

        /// Explicit normal prior for muzzle velocity, MEAN:SIGMA. Values follow
        /// --units (fps for imperial, m/s for metric). Never inferred from
        /// --chrono-velocity.
        #[arg(long, value_name = "MEAN:SIGMA")]
        mv_prior: Option<String>,

        /// Explicit normal prior for the scalar ballistic coefficient,
        /// MEAN:SIGMA (dimensionless).
        #[arg(long, value_name = "MEAN:SIGMA")]
        bc_prior: Option<String>,

        /// Range at which to report a propagated model-drop interval. Repeatable;
        /// follows --units (yards/meters).
        #[arg(long, value_name = "RANGE", action = clap::ArgAction::Append)]
        predict_range: Vec<f64>,

        /// Optional one-standard-deviation error for a future impact reading, in
        /// --drop-unit. Adds a future-observation band alongside the latent model
        /// band at every --predict-range.
        #[arg(long, value_name = "SIGMA")]
        prediction_sigma: Option<f64>,

        /// Unit for --measured-drop and every --observed drop in multi-observation
        /// mode. Ignored in single-observation mode (which is always MIL).
        #[arg(long = "drop-unit", default_value = "mil")]
        drop_unit: DropUnit,

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

        /// Sight height above bore (inches for imperial, mm for metric) [default: 2 in / 50 mm]
        #[arg(long)]
        sight_height: Option<f64>,

        /// Temperature (Fahrenheit for imperial, Celsius for metric; default 59 F / 15 C)
        #[arg(long)]
        temperature: Option<f64>,

        /// Pressure (inHg for imperial, hPa for metric; default 29.92 inHg / 1013.25 hPa)
        #[arg(long)]
        pressure: Option<f64>,

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
        #[arg(
            long,
            value_name = "DIR",
            help = "Use caliber-specific 5D BC correction tables"
        )]
        bc_table_dir: Option<PathBuf>,

        /// Enable auto-download of BC5D correction tables when needed
        #[cfg(feature = "online")]
        #[arg(long, help = "Auto-download BC5D tables when needed")]
        bc_table_auto: bool,

        /// Base URL for BC5D table downloads
        #[cfg(feature = "online")]
        #[arg(
            long,
            default_value = "https://ballistics.tools/downloads/bc5d",
            help = "Base URL for BC5D table downloads"
        )]
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

    /// Recommend target ranges for an identifiable MV + BC truing experiment
    PlanTruing {
        /// Load nominal load, rifle, and atmosphere values from a saved profile.
        /// Explicit flags override profile values.
        #[arg(long, value_name = "NAME")]
        profile: Option<String>,

        /// Nominal muzzle velocity (fps or m/s based on global --units)
        #[arg(short = 'v', long, value_parser = f64_range(0.0, 6000.0))]
        velocity: Option<f64>,

        /// Nominal scalar ballistic coefficient
        #[arg(short = 'b', long, value_parser = f64_range(0.001, 2.0))]
        bc: Option<f64>,

        /// Scalar drag model (G1 or G7)
        #[arg(long)]
        drag_model: Option<DragModelArg>,

        /// Bullet weight (grains or grams based on global --units)
        #[arg(short = 'm', long, value_parser = f64_range(0.1, 2000.0))]
        mass: Option<f64>,

        /// Bullet diameter (inches or millimeters based on global --units)
        #[arg(short = 'd', long, value_parser = f64_range(0.01, 60.0))]
        diameter: Option<f64>,

        /// Comma-separated available target ranges in global --units
        #[arg(
            long,
            value_name = "R1,R2,...",
            value_delimiter = ',',
            num_args = 1..,
            conflicts_with = "range_grid",
            required_unless_present = "range_grid"
        )]
        candidate_ranges: Vec<f64>,

        /// Explicit candidate grid START:END:STEP in global --units
        #[arg(
            long,
            value_name = "START:END:STEP",
            conflicts_with = "candidate_ranges",
            required_unless_present = "candidate_ranges"
        )]
        range_grid: Option<String>,

        /// Number of impact observations to collect
        #[arg(long, default_value_t = 3)]
        observation_count: usize,

        /// Required separation between selected target ranges in global --units
        #[arg(long, default_value_t = 0.0, value_parser = f64_range(0.0, 100000.0))]
        minimum_separation: f64,

        /// Assumed independent one-standard-deviation impact-reading error in
        /// --drop-unit
        #[arg(long, value_parser = f64_range(f64::MIN_POSITIVE, 1000000.0))]
        measurement_resolution: f64,

        /// Unit used for the assumed impact-reading error and predicted drop
        #[arg(long, default_value = "mil")]
        drop_unit: DropUnit,

        /// Zero distance in global --units
        #[arg(long)]
        zero_distance: Option<f64>,

        /// Sight height (inches or millimeters based on global --units)
        #[arg(long)]
        sight_height: Option<f64>,

        /// Temperature (Fahrenheit or Celsius based on global --units)
        #[arg(long)]
        temperature: Option<f64>,

        /// Pressure (inHg or hPa based on global --units)
        #[arg(long)]
        pressure: Option<f64>,

        /// Relative humidity (0-100 percent)
        #[arg(long, value_parser = f64_range(0.0, 100.0))]
        humidity: Option<f64>,

        /// Altitude (feet or meters based on global --units)
        #[arg(long)]
        altitude: Option<f64>,

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

        /// Temperature (Fahrenheit or Celsius based on --units; default 59 F / 15 C)
        #[arg(long)]
        temperature: Option<f64>,

        /// Pressure (inHg or hPa based on --units; default 29.92 inHg / 1013.25 hPa)
        #[arg(long)]
        pressure: Option<f64>,

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
        #[arg(long, default_value = "50.0", value_parser = f64_range(0.001, 100000.0))]
        step: f64,

        /// Adjustment unit (mil, moa, smoa, iphy, or clicks). `clicks` requires
        /// --elevation-click-value (or a saved profile's elevation_click) and formats
        /// whole turret clicks instead of an angle.
        #[arg(long, default_value = "mil")]
        adjustment_unit: AdjustmentUnit,

        /// Turret elevation click graduation for `--adjustment-unit clicks`, e.g. 0.1mil or
        /// 0.25moa (falls back to the saved profile's elevation_click when omitted).
        #[arg(long)]
        elevation_click_value: Option<String>,

        /// Turret windage click graduation for `--adjustment-unit clicks` (accepted for parity
        /// with `trajectory`; come-ups has no windage column to apply it to).
        #[arg(long)]
        windage_click_value: Option<String>,

        /// Sight height above bore (inches for imperial, mm for metric)
        #[arg(long)]
        sight_height: Option<f64>,

        /// Temperature (Fahrenheit or Celsius based on --units; default 59 F / 15 C)
        #[arg(long)]
        temperature: Option<f64>,

        /// Pressure (inHg or hPa based on --units; default 29.92 inHg / 1013.25 hPa)
        #[arg(long)]
        pressure: Option<f64>,

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

    /// Moving-target lead table (hold in the direction of target travel)
    Lead {
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

        /// Sight height above bore (inches for imperial, mm for metric)
        #[arg(long)]
        sight_height: Option<f64>,

        /// Temperature (Fahrenheit or Celsius based on --units; default 59 F / 15 C)
        #[arg(long)]
        temperature: Option<f64>,

        /// Pressure (inHg or hPa based on --units; default 29.92 inHg / 1013.25 hPa)
        #[arg(long)]
        pressure: Option<f64>,

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

        /// Enable powder temperature sensitivity
        #[arg(long)]
        use_powder_sensitivity: bool,

        /// Powder temperature sensitivity (fps/°F imperial, m/s/°C metric).
        /// Defaults to 1.0 fps/°F (0.54864 m/s/°C).
        #[arg(long)]
        powder_temp_sensitivity: Option<f64>,

        /// Powder temperature (°F/°C). With --powder-temp-curve, this is the powder
        /// temperature the curve is interpolated at to resolve muzzle velocity (defaults
        /// to --temperature when omitted, i.e. powder assumed at air temperature). With
        /// --powder-temp-sensitivity (linear model), it is the reference temperature the
        /// stated velocity was measured at (defaults to the 70°F equivalent).
        #[arg(long)]
        powder_temp: Option<f64>,

        /// Measured powder-temperature -> muzzle-velocity curve as comma-separated
        /// TEMP:VELOCITY points, e.g. "40:2620,70:2700,100:2760" (temps in F/C,
        /// velocities in fps / m·s⁻¹ per --units). The muzzle velocity is interpolated
        /// from this table at the powder temperature (--powder-temp, or --temperature if
        /// unset; clamped at the endpoints, no extrapolation). Data-driven, non-linear
        /// alternative that OVERRIDES --powder-temp-sensitivity when supplied.
        #[arg(long = "powder-temp-curve", value_name = "TEMP:VEL,...")]
        powder_temp_curve: Option<String>,

        /// Target speed (mph for imperial, m/s for metric)
        #[arg(long, value_parser = f64_range(0.0, 300.0))]
        target_speed: f64,

        /// Direction of target travel relative to the line of sight, degrees:
        /// 0 = directly away, 90 = crossing left-to-right, 180 = directly toward,
        /// 270 = crossing right-to-left
        #[arg(long, default_value = "90.0")]
        target_angle: f64,

        /// Target length for body-length holds (inches for imperial, mm for metric)
        #[arg(long, value_parser = f64_range(0.001, 10000.0))]
        target_length: Option<f64>,

        /// Start range (yards or meters)
        #[arg(long, default_value = "100.0")]
        start: f64,

        /// End range (yards or meters)
        #[arg(long, default_value = "600.0")]
        end: f64,

        /// Range step (yards or meters)
        #[arg(long, default_value = "100.0", value_parser = f64_range(0.001, 100000.0))]
        step: f64,

        /// Adjustment unit (mil or moa)
        #[arg(long, default_value = "mil")]
        adjustment_unit: AdjustmentUnit,

        /// Output format
        #[arg(short = 'o', long, default_value = "table")]
        output: OutputFormat,
    },

    /// Resolve the powder-temperature-adjusted muzzle velocity without running a
    /// trajectory (MBA-737): the linear sensitivity model or a measured
    /// temperature->velocity curve, with an optional sweep table across temperatures.
    Powder {
        /// Nominal muzzle velocity (fps or m/s based on --units) at the reference
        /// temperature. Required for the linear model. With --powder-temp-curve the
        /// curve supplies the velocity and this only anchors the reported shift.
        #[arg(short = 'v', long, value_parser = f64_range(0.0, 6000.0))]
        velocity: Option<f64>,

        /// Shot-day air temperature (°F/°C based on --units). Drives the linear
        /// model's delta from the reference temperature; with --powder-temp-curve it
        /// is the curve's lookup temperature when --powder-temp is not given (powder
        /// assumed at air temperature).
        #[arg(long)]
        temperature: Option<f64>,

        /// Powder temperature sensitivity (fps/°F imperial, m/s/°C metric).
        /// Defaults to 1.0 fps/°F (0.54864 m/s/°C).
        #[arg(long)]
        powder_temp_sensitivity: Option<f64>,

        /// Powder temperature (°F/°C). With --powder-temp-curve, this is the powder
        /// temperature the curve is interpolated at to resolve muzzle velocity (defaults
        /// to --temperature when omitted, i.e. powder assumed at air temperature). With
        /// the linear model, it is the reference temperature the stated velocity was
        /// measured at (defaults to the 70°F equivalent).
        #[arg(long)]
        powder_temp: Option<f64>,

        /// Measured powder-temperature -> muzzle-velocity curve as comma-separated
        /// TEMP:VELOCITY points, e.g. "40:2620,70:2700,100:2760" (temps in F/C,
        /// velocities in fps / m·s⁻¹ per --units). The muzzle velocity is interpolated
        /// from this table at the powder temperature (--powder-temp, or --temperature if
        /// unset; clamped at the endpoints, no extrapolation). Data-driven, non-linear
        /// alternative that OVERRIDES --powder-temp-sensitivity when supplied.
        #[arg(long = "powder-temp-curve", value_name = "TEMP:VEL,...")]
        powder_temp_curve: Option<String>,

        /// Print a velocity table across a temperature range instead of a single
        /// resolution: START:END:STEP in °F (imperial) or °C (metric), e.g.
        /// "20:110:10". With a curve, the sweep temperatures are powder temperatures.
        #[arg(long, value_name = "START:END:STEP")]
        sweep: Option<String>,

        /// Bullet mass (grains imperial / grams metric). When given, muzzle energy
        /// (ft·lb / J) is reported alongside each resolved velocity.
        #[arg(short = 'm', long, value_parser = f64_range(0.1, 2000.0))]
        mass: Option<f64>,

        /// Output format (table, json, or csv)
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

        /// Wind angle in DEGREES, wind-FROM convention (matches --wind-direction):
        /// 0 = headwind, 90 = from the right (full value), 180 = tailwind,
        /// 270 = from the left. Positive drift = dial right. Default 90.
        #[arg(long, value_parser = f64_range(0.0, 360.0), conflicts_with = "wind_angles")]
        wind_angle: Option<f64>,

        /// Comma-separated wind angles (degrees, wind-FROM) — emits one complete
        /// card per angle. Example: --wind-angles 30,60,90
        #[arg(long, conflicts_with = "wind_angle")]
        wind_angles: Option<String>,

        /// Start range (yards or meters)
        #[arg(long, default_value = "100.0")]
        start: f64,

        /// End range (yards or meters)
        #[arg(long, default_value = "1000.0")]
        end: f64,

        /// Range step (yards or meters)
        #[arg(long, default_value = "100.0", value_parser = f64_range(0.001, 100000.0))]
        step: f64,

        /// Adjustment unit (mil or moa)
        #[arg(long, default_value = "mil")]
        adjustment_unit: AdjustmentUnit,

        /// Sight height above bore (inches for imperial, mm for metric)
        #[arg(long)]
        sight_height: Option<f64>,

        /// Temperature (Fahrenheit or Celsius based on --units; default 59 F / 15 C)
        #[arg(long)]
        temperature: Option<f64>,

        /// Pressure (inHg or hPa based on --units; default 29.92 inHg / 1013.25 hPa)
        #[arg(long)]
        pressure: Option<f64>,

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

        /// Temperature (Fahrenheit or Celsius based on --units; default 59 F / 15 C)
        #[arg(long)]
        temperature: Option<f64>,

        /// Pressure (inHg or hPa based on --units; default 29.92 inHg / 1013.25 hPa)
        #[arg(long)]
        pressure: Option<f64>,

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
        #[arg(long, default_value = "50.0", value_parser = f64_range(0.001, 100000.0))]
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

        /// Temperature (Fahrenheit or Celsius based on --units; default 59 F / 15 C)
        #[arg(long)]
        temperature: Option<f64>,

        /// Pressure (inHg or hPa based on --units; default 29.92 inHg / 1013.25 hPa)
        #[arg(long)]
        pressure: Option<f64>,

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

    /// Compare multiple loads side-by-side at the same conditions (MBA-735)
    Compare {
        /// A load as NAME:DRAG:BC:MASS:VELOCITY with an optional sixth :DIAMETER
        /// field (repeat 2-8 times).
        /// DRAG is g1|g7; MASS is grains/grams, VELOCITY fps|m/s, DIAMETER in|mm
        /// per --units (diameter defaults to .308 in / 7.82 mm). NAME must not
        /// contain ':'. Mixable with --profile.
        #[arg(long = "load", value_name = "SPEC")]
        loads: Vec<String>,

        /// A saved profile as a load (repeatable; mixable with --load). Like
        /// range-table, a profile's bc_segments/drag_curve are not yet consumed
        /// here — the scalar BC is used.
        #[arg(long = "profile", value_name = "NAME")]
        profiles: Vec<String>,

        /// Shared zero distance (yards for imperial, meters for metric); every
        /// load is zeroed independently at this distance
        #[arg(long)]
        zero_distance: f64,

        /// Start range (yards or meters)
        #[arg(long, default_value = "100.0")]
        start: f64,

        /// End range (yards or meters)
        #[arg(long, default_value = "500.0")]
        end: f64,

        /// Range step (yards or meters)
        #[arg(long, default_value = "100.0", value_parser = f64_range(0.001, 100000.0))]
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
        #[arg(long, default_value = "1.5")]
        sight_height: f64,

        /// Temperature (F or C based on --units)
        #[arg(long)]
        temperature: Option<f64>,

        /// Pressure (inHg or hPa based on --units)
        #[arg(long)]
        pressure: Option<f64>,

        /// Humidity (percent)
        #[arg(long, default_value = "50.0")]
        humidity: f64,

        /// Altitude (feet or meters based on --units)
        #[arg(long, default_value = "0.0")]
        altitude: f64,

        /// Output format (table, json, csv)
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

#[allow(
    clippy::large_enum_variant,
    reason = "variants intentionally mirror the stable flat Clap command shape"
)]
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

        /// Barrel twist rate (inches per turn for imperial, mm per turn for metric)
        #[arg(long)]
        twist_rate: Option<f64>,

        /// Sight height above bore (inches for imperial, mm for metric)
        #[arg(long)]
        sight_height: Option<f64>,

        /// Default zero distance (yards for imperial, meters for metric)
        #[arg(long)]
        zero_distance: Option<f64>,

        /// Default temperature (Fahrenheit for imperial, Celsius for metric; default 59 F / 15 C)
        #[arg(long)]
        temperature: Option<f64>,

        /// Default pressure (inHg for imperial, hPa for metric; default 29.92 inHg / 1013.25 hPa)
        #[arg(long)]
        pressure: Option<f64>,

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

        /// Bullet length (inches for imperial, mm for metric; used for BC table lookup)
        #[arg(long)]
        bullet_length: Option<f64>,

        /// Turret elevation click graduation for `--adjustment-unit clicks` (MBA-1355), e.g. 0.1mil
        /// or 0.25moa. Validated at save time so a profile can't store garbage.
        #[arg(long)]
        elevation_click: Option<String>,

        /// Turret windage click graduation for `--adjustment-unit clicks` (MBA-1355); falls back to
        /// the elevation graduation when omitted, e.g. 0.1mil or 0.25moa.
        #[arg(long)]
        windage_click: Option<String>,
    },

    /// Import a profile from a third-party file (.a7p — ArcherBC2 format)
    #[cfg(feature = "profile-import")]
    Import {
        /// Path to the .a7p file
        file: PathBuf,

        /// Profile name to save under (defaults to the file's profile name, sanitized)
        #[arg(long)]
        name: Option<String>,

        /// Parse and show the full mapping report without saving anything
        #[arg(long)]
        dry_run: bool,

        /// Treat a checksum mismatch as a fatal error instead of a warning
        #[arg(long)]
        strict: bool,
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

// UnitSystem and DragModelArg moved to ballistics_engine::truing (MBA-1343).

#[derive(Debug, Clone, Copy, ValueEnum)]
enum OutputFormat {
    Json,
    Csv,
    Table,
    Pdf,
}

// DropUnit moved to ballistics_engine::truing (MBA-1343).

#[derive(Debug, Clone, Copy, ValueEnum)]
enum MonteCarloOutput {
    Summary,
    Full,
    Statistics,
}

/// Wind-shear boundary-layer profile model. Parsed as a clap enum so invalid
/// values are rejected at the command line (MBA-965). Both snake_case (the
/// historical/engine spelling) and kebab-case are accepted for each model.
#[derive(Debug, Clone, Copy, PartialEq, ValueEnum)]
enum WindShearModelArg {
    None,
    #[value(name = "logarithmic")]
    Logarithmic,
    #[value(
        name = "power_law",
        alias = "power-law",
        alias = "powerlaw",
        alias = "exponential"
    )]
    PowerLaw,
    #[value(name = "ekman_spiral", alias = "ekman-spiral", alias = "ekman")]
    EkmanSpiral,
    #[value(name = "custom_layers", alias = "custom-layers", alias = "custom")]
    CustomLayers,
}

impl WindShearModelArg {
    /// Canonical lower-snake string understood by the engine
    /// (`cli_api`/`wind_shear`).
    fn as_engine_str(self) -> &'static str {
        match self {
            WindShearModelArg::None => "none",
            WindShearModelArg::Logarithmic => "logarithmic",
            WindShearModelArg::PowerLaw => "power_law",
            WindShearModelArg::EkmanSpiral => "ekman_spiral",
            WindShearModelArg::CustomLayers => "custom_layers",
        }
    }
}

#[derive(Debug, Clone, Copy, ValueEnum, Default, PartialEq)]
enum AdjustmentUnit {
    /// Milliradians (1 MIL = 3.6 inches at 100 yards)
    #[default]
    Mil,
    /// Minutes of Angle, true MOA (1 MOA = 1.047 inches at 100 yards)
    Moa,
    /// Shooter's MOA (exactly 1 inch per 100 yards)
    Smoa,
    /// Inches per hundred yards (numerically identical to SMOA)
    Iphy,
    /// Whole turret clicks; requires an elevation click graduation
    /// (--elevation-click-value or the profile's elevation_click)
    Clicks,
}

/// Terminal chart renderer for `trajectory --plot` (MBA-1320). Bare `--plot` (no value)
/// resolves to `Braille` via clap's `default_missing_value`; see the `plot` field on
/// `Commands::Trajectory`.
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
enum PlotStyle {
    /// Unicode braille-dot canvas (2x4 dots per character cell). The default terminal
    /// chart renderer.
    Braille,
    /// ASCII '*'-per-cell canvas, for terminals/fonts without braille glyph coverage.
    Ascii,
}

impl PlotStyle {
    fn as_canvas_style(self) -> terminal_plot::CanvasStyle {
        match self {
            PlotStyle::Braille => terminal_plot::CanvasStyle::Braille,
            PlotStyle::Ascii => terminal_plot::CanvasStyle::Ascii,
        }
    }
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
    /// Turret elevation click graduation for `--adjustment-unit clicks` (MBA-1355), e.g. "0.1mil" or
    /// "0.25moa" — parsed by `parse_click_value` at both save-time (validation) and
    /// resolve-time (`resolve_click_values`). Unit-invariant (an angular graduation, not a
    /// linear measurement), so `converted_to` leaves it untouched — see the `bc_segments`/
    /// `drag_curve` comment below for why.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    elevation_click: Option<String>,
    /// Turret windage click graduation for `--adjustment-unit clicks` (MBA-1355). Falls back to the
    /// resolved elevation graduation when unset — see `resolve_click_values`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    windage_click: Option<String>,
    /// Velocity-banded BC breakpoints (MBA-1323 Phase 2: multi-row `.a7p` G1/G7 import).
    /// `velocity_mps` in each entry is ALWAYS meters/second regardless of this profile's
    /// `units` field — see [`ProfileBcSegment`]. The scalar `bc` field above is retained as
    /// the highest-velocity row for tools that only understand a single BC; this list is the
    /// authoritative full schedule once present. `None` (the pre-Phase-2 shape) means "no
    /// velocity-banded schedule was captured" and callers fall back to the scalar `bc`.
    ///
    /// FORWARD-COMPAT WARNING (one-way): `#[serde(default)]` means this field round-trips
    /// safely through readers that predate Phase 2, but "safely" only means *deserialization*
    /// doesn't error — a pre-Phase-2 (or otherwise un-updated, e.g. stale WASM/bindings) reader
    /// silently drops this key and solves with only the scalar `bc` above. That is a materially
    /// different, unwarned answer whenever the schedule's non-fastest bands matter (empirically
    /// confirmed: ~639 m/s vs. ~411 m/s impact velocity for the same imported profile — see
    /// CLI_USAGE.md's "Multi-BC and CUSTOM drag curves" section). There is no sentinel trick
    /// available here the way there is for `drag_curve`/CUSTOM below (a real, plausible-looking
    /// `bc` value is unavoidable for back-compat with single-BC tools), so this direction of
    /// version skew degrades silently by design and must stay documented rather than "fixed".
    #[serde(default, skip_serializing_if = "Option::is_none")]
    bc_segments: Option<Vec<ProfileBcSegment>>,
    /// Full Mach/Cd drag curve (MBA-1323 Phase 2: `.a7p` `bc_type == CUSTOM` import). When
    /// present, the scalar `bc`/`drag_model` fields are not physically meaningful for the
    /// solve (`drag_model` reads "CUSTOM"; see `map_a7p_to_profile`'s CUSTOM handling for why
    /// `bc` is an inert `0.0` sentinel rather than a real coefficient).
    ///
    /// FORWARD-COMPAT NOTE: unlike `bc_segments` above, a reader that predates Phase 2 (or
    /// otherwise doesn't consume this field) is safe by construction, not just by omission: it
    /// still sees `bc == 0.0` and `drag_model == "CUSTOM"`, so `BallisticInputs::validate_for_solve`
    /// rejects the solve loudly ("bc_value must be finite and greater than zero") instead of
    /// silently running physics under a fabricated coefficient.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    drag_curve: Option<Vec<ProfileDragPoint>>,
}

/// One velocity-banded BC breakpoint (profile schema v2, MBA-1323 Phase 2). Stored as a raw
/// breakpoint, NOT a pre-computed `velocity_min`/`velocity_max` band — banding into the
/// engine's [`BCSegmentData`] shape happens at solve time (`bc_segments_from_profile`), so the
/// stored JSON stays a simple, order-independent list that round-trips cleanly.
///
/// `velocity_mps` is ALWAYS meters/second, independent of [`ProfileData::units`]. This is
/// intentional (matches the engine's internal BC-segment plumbing, which is also always in
/// engine units) and is why [`ProfileData::converted_to`] leaves this field untouched — see
/// the comment there.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct ProfileBcSegment {
    bc: f64,
    velocity_mps: f64,
}

/// One Mach/Cd point of a saved custom drag curve (profile schema v2, MBA-1323 Phase 2).
/// Both fields are unit-invariant (Mach is dimensionless; Cd is dimensionless), so
/// [`ProfileData::converted_to`] leaves `drag_curve` untouched too.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct ProfileDragPoint {
    mach: f64,
    cd: f64,
}

fn default_unit_system() -> String {
    "imperial".to_string()
}
fn default_temperature() -> f64 {
    59.0
}
fn default_pressure() -> f64 {
    29.92
}
fn default_humidity() -> f64 {
    50.0
}

#[derive(Debug, Serialize, Deserialize)]
struct TrajectoryPoint {
    time: f64,
    x: f64,
    y: f64,
    z: f64,
    velocity: f64,
    energy: f64,
    /// Mover-ring linear radius, meters (MBA-1325): `target_speed_mps * time`. Present
    /// only when `--target-speed` was supplied (> 0); unit is fixed at meters regardless
    /// of `--units` (units-in-the-name, see MBA-1315).
    #[serde(skip_serializing_if = "Option::is_none")]
    mover_ring_m: Option<f64>,
    /// Mover-ring angular radius, milliradians: `mover_ring_m / downrange_m * 1000`.
    /// Omitted at the muzzle (downrange = 0, ratio undefined) and whenever
    /// `mover_ring_m` itself is absent.
    #[serde(skip_serializing_if = "Option::is_none")]
    mover_ring_mil: Option<f64>,
}

#[derive(Debug, Serialize, Deserialize)]
struct TrajectoryResult {
    /// Unit system the numeric fields are expressed in ("imperial" or "metric").
    units: String,
    max_range: f64,
    max_height: f64,
    time_of_flight: f64,
    impact_velocity: f64,
    impact_energy: f64,
    stability_coefficient: Option<f64>,
    spin_drift: Option<f64>,
    /// Minimum pitch-damping coefficient (only when --enable-pitch-damping).
    #[serde(skip_serializing_if = "Option::is_none")]
    min_pitch_damping: Option<f64>,
    /// Mach number when entering the transonic regime (pitch-damping diagnostic).
    #[serde(skip_serializing_if = "Option::is_none")]
    transonic_mach: Option<f64>,
    /// Maximum yaw angle during flight, radians (only when --enable-precession).
    #[serde(skip_serializing_if = "Option::is_none")]
    max_yaw_angle: Option<f64>,
    /// Maximum precession angle during flight, radians (only when --enable-precession).
    #[serde(skip_serializing_if = "Option::is_none")]
    max_precession_angle: Option<f64>,
    /// Final pitch angle, radians (only when --enable-precession).
    #[serde(skip_serializing_if = "Option::is_none")]
    final_pitch_angle: Option<f64>,
    /// Final yaw angle, radians (only when --enable-precession).
    #[serde(skip_serializing_if = "Option::is_none")]
    final_yaw_angle: Option<f64>,
    trajectory: Vec<TrajectoryPoint>,
    /// Self-describing metadata for this document (MBA-1315), additive only: every key/value
    /// above and every `trajectory[]` point key are unchanged. Appended last so the
    /// pre-existing fields serialize byte-identically ahead of it. See `TrajectoryLegend` for
    /// how the unit/axis text was derived (empirically, from controlled crosswind/no-wind
    /// runs) and CLI_USAGE.md's JSON Format section for the same text with worked examples.
    legend: TrajectoryLegend,
}

/// Concrete unit labels and `trajectory[]` axis semantics for one [`TrajectoryResult`]
/// document (MBA-1315). A field tester's tooling once misread `trajectory[].x`/`z` as feet
/// and swapped their meaning because the legacy JSON never stated either; this block states
/// both. Nested under one `legend` key rather than top-level `units`/`axes` keys because
/// `units` already names the pre-existing "imperial"/"metric" string field on
/// [`TrajectoryResult`].
#[derive(Debug, Serialize, Deserialize)]
struct TrajectoryLegend {
    /// Concrete unit abbreviation for each numeric quantity family in this document.
    units: TrajectoryLegendUnits,
    /// `trajectory[]` point axis semantics, verified empirically (not assumed): a pure
    /// crosswind run and a no-wind run, comparing which point component moved, in which sign
    /// direction, against the table output (see `tests/legacy_trajectory_json_contract.rs`).
    axes: TrajectoryLegendAxes,
}

/// Unit abbreviation for each numeric quantity family. `distance` applies to
/// `trajectory[].x`/`y`/`z`, `max_range`, `max_height`, and `spin_drift`; `velocity` to
/// `trajectory[].velocity` and `impact_velocity`; `energy` to `trajectory[].energy` and
/// `impact_energy`. Angle fields (`max_yaw_angle`, `max_precession_angle`,
/// `final_pitch_angle`, `final_yaw_angle`) and time fields are unit-invariant (always
/// radians / seconds) and are not covered here.
#[derive(Debug, Serialize, Deserialize)]
struct TrajectoryLegendUnits {
    distance: String,
    velocity: String,
    energy: String,
}

/// Axis semantics for `trajectory[]` points. `x` is lateral and `z` is downrange in THIS
/// output — the reverse of the x=lateral/z=depth convention some 3D tooling assumes, and the
/// original source of the field-tester misdiagnosis this ticket fixes.
#[derive(Debug, Serialize, Deserialize)]
struct TrajectoryLegendAxes {
    x: String,
    y: String,
    z: String,
}

impl TrajectoryLegend {
    fn for_units(units: UnitSystem) -> Self {
        let (distance, velocity, energy) = match units {
            UnitSystem::Imperial => ("yd", "fps", "ft-lb"),
            UnitSystem::Metric => ("m", "m/s", "J"),
        };
        TrajectoryLegend {
            units: TrajectoryLegendUnits {
                distance: distance.to_string(),
                velocity: velocity.to_string(),
                energy: energy.to_string(),
            },
            axes: TrajectoryLegendAxes {
                x: "lateral offset from the muzzle's initial aiming direction; positive means \
                    to the shooter's right (e.g. a crosswind FROM the left, --wind-direction \
                    270, drifts x positive; FROM the right, --wind-direction 90, drifts x \
                    negative). Zero at the muzzle."
                    .to_string(),
                y: "height above the ground in the world frame; positive means up. Starts near \
                    bore height and reaches 0 at ground impact. This is NOT height above the \
                    line of sight (compare solve-json v1's LOS-relative drop_m)."
                    .to_string(),
                z: "downrange distance from the muzzle; zero at the muzzle, always increasing."
                    .to_string(),
            },
        }
    }
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
    bullet_length: f64,
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
    wind_direction: f64, // degrees at the CLI boundary
    wind_vertical: f64,  // m/s, positive = updraft
    // Downrange-segmented wind (engine units: speed km/h, angle deg, distance m).
    // When non-empty, overrides the scalar wind above.
    wind_segments: Vec<ballistics_engine::wind::WindSegment>,

    // Output
    output: OutputFormat,
    full: bool,
    // Inline terminal chart (MBA-1320); None = no chart (the default). Only rendered for
    // OutputFormat::Table — see run_trajectory's OutputFormat::Table arm.
    plot: Option<PlotStyle>,
    units: UnitSystem,

    // Sight / bore geometry (metric)
    sight_height: f64,
    bore_height: f64,
    ignore_ground_impact: bool,
    // Rifle cant, degrees (CLI boundary; converted to radians for BallisticInputs.cant_angle).
    cant: f64,

    // BC options
    use_bc_segments: bool,
    use_cluster_bc: bool,
    bc_segments_data: Option<Vec<BCSegmentData>>,
    // Custom drag deck (--drag-table), takes precedence over the G-model + BC when set.
    custom_drag_table: Option<ballistics_engine::drag::DragTable>,

    // Advanced physics toggles
    enable_magnus: bool,
    enable_coriolis: bool,
    enable_spin_drift: bool,
    enable_aerodynamic_jump: bool,
    enable_wind_shear: bool,
    wind_shear_model: WindShearModelArg,
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
    shot_direction: Option<f64>, // compass bearing of the shot, degrees, 0=N (for Coriolis)
    shooting_angle: f64,

    // Powder sensitivity
    use_powder_sensitivity: bool,
    powder_temp_sensitivity: f64,
    powder_temp: f64,
    // Optional measured (temperature_celsius, muzzle_velocity_m_s) curve; SI, sorted.
    powder_temp_curve: Option<Vec<(f64, f64)>>,
    // Powder temperature (Celsius) to interpolate the curve at; None = ambient temp.
    powder_curve_temp_c: Option<f64>,

    // Mover ring (MBA-1325): target speed in display units (mph imperial / m/s metric,
    // same convention as `lead --target-speed`). 0.0 (default) disables the per-point
    // Ring column/fields in every output format; also feeds the PDF Lead column.
    target_speed: f64,
    // Angular unit for the ring TABLE column (mil or moa) — from --adjustment-unit,
    // which trajectory already exposes for the PDF dope card. CSV keeps ring_mil and
    // JSON keeps mover_ring_m/mover_ring_mil regardless: those carry the unit in the
    // name (the contract), the table is the human display surface.
    adjustment_unit: AdjustmentUnit,
    // Resolved turret click graduations for `--adjustment-unit clicks` (MBA-1355): Some(...) iff
    // adjustment_unit == Clicks (resolved once, eagerly, in the Trajectory command's
    // dispatch — see resolve_click_values); None for every other adjustment_unit. Drives
    // the mover Ring column and the PDF dope-card rows.
    elevation_click: Option<ClickValue>,
    windage_click: Option<ClickValue>,

    // PDF metadata
    pdf_metadata: Option<PdfMetadata>,
}

#[derive(Debug, Serialize, Deserialize)]
struct MonteCarloResult {
    mean_range: f64,
    std_range: f64,
    mean_impact_velocity: f64,
    std_impact_velocity: f64,
    /// Median radial miss among simulations that reached the target plane; `None` when no
    /// simulation reached it.
    cep: Option<f64>,
    target_shortfall_fraction: f64,
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
            UnitSystem::Imperial => val * GRAINS_TO_KG, // grains to kg
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

    /// Convert a temperature *difference* (delta), not an absolute point.
    /// A 1 F delta equals a 5/9 C delta (no -32 offset). Used to convert
    /// per-degree quantities such as powder temperature sensitivity (MBA-963).
    fn temperature_delta_to_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val * 5.0 / 9.0, // Fahrenheit delta to Celsius delta
        }
    }

    fn pressure_to_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val * 33.8639, // inHg to hPa
        }
    }

    /// Resolve the pressure CLI argument AFTER --units is known.
    ///
    /// `None` -> the per-unit standard atmosphere (29.92 inHg / 1013.25 hPa).
    /// `Some(v)` -> validated against a plausible range in the user's own unit
    /// and returned IN USER UNITS (the caller still applies `pressure_to_metric`).
    fn resolve_pressure(val: Option<f64>, units: UnitSystem) -> Result<f64, String> {
        match (val, units) {
            (None, UnitSystem::Imperial) => Ok(29.92),
            (None, UnitSystem::Metric) => Ok(1013.25),
            (Some(v), UnitSystem::Imperial) => {
                if !(8.0..=33.0).contains(&v) {
                    Err(format!(
                        "--pressure {v} inHg is out of range (expected ~8..33 inHg for imperial units)"
                    ))
                } else {
                    Ok(v)
                }
            }
            (Some(v), UnitSystem::Metric) => {
                if !(250.0..=1100.0).contains(&v) {
                    Err(format!(
                        "--pressure {v} hPa is out of range (expected ~250..1100 hPa for metric units)"
                    ))
                } else {
                    Ok(v)
                }
            }
        }
    }

    /// Resolve the temperature CLI argument AFTER --units is known.
    ///
    /// `None` -> the per-unit standard temperature (59 F / 15 C).
    /// `Some(v)` -> validated against a plausible range in the user's own unit
    /// and returned IN USER UNITS (the caller still applies `temperature_to_metric`).
    fn resolve_temperature(val: Option<f64>, units: UnitSystem) -> Result<f64, String> {
        match (val, units) {
            (None, UnitSystem::Imperial) => Ok(59.0),
            (None, UnitSystem::Metric) => Ok(15.0),
            (Some(v), UnitSystem::Imperial) => {
                if !(-148.0..=392.0).contains(&v) {
                    Err(format!(
                        "--temperature {v} F is out of range (expected ~-148..392 F for imperial units)"
                    ))
                } else {
                    Ok(v)
                }
            }
            (Some(v), UnitSystem::Metric) => {
                if !(-100.0..=200.0).contains(&v) {
                    Err(format!(
                        "--temperature {v} C is out of range (expected ~-100..200 C for metric units)"
                    ))
                } else {
                    Ok(v)
                }
            }
        }
    }

    fn altitude_to_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val * 0.3048, // feet to meters
        }
    }

    // target_size_to_metric moved to ballistics_engine::wez (MBA-1343 Phase B).

    // Output conversions (from metric)
    fn velocity_from_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val / 0.3048, // m/s to fps
        }
    }

    fn mass_from_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val * 1000.0,          // kg to grams
            UnitSystem::Imperial => val / GRAINS_TO_KG, // kg to grains
        }
    }

    fn diameter_from_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val * 1000.0,   // meters to mm
            UnitSystem::Imperial => val / 0.0254, // meters to inches
        }
    }

    fn sight_height_from_metric(val: f64, units: UnitSystem) -> f64 {
        Self::diameter_from_metric(val, units)
    }

    fn distance_from_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val / 0.9144, // meters to yards
        }
    }

    fn wind_from_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val / 0.44704, // m/s to mph
        }
    }

    fn temperature_from_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val * 9.0 / 5.0 + 32.0, // Celsius to Fahrenheit
        }
    }

    fn pressure_from_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val / 33.8639, // hPa to inHg
        }
    }

    fn altitude_from_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val / 0.3048, // meters to feet
        }
    }

    fn energy_from_metric(val: f64, units: UnitSystem) -> f64 {
        match units {
            UnitSystem::Metric => val,
            UnitSystem::Imperial => val * 0.737562, // Joules to ft-lbs
        }
    }
}

impl ProfileData {
    /// Convert a loaded profile's display-unit values into the active CLI unit system before
    /// merging them with explicit arguments. The persisted profile remains unchanged.
    fn converted_to(mut self, target: UnitSystem) -> Result<Self, Box<dyn Error>> {
        let source = match self.units.trim().to_ascii_lowercase().as_str() {
            "imperial" => UnitSystem::Imperial,
            "metric" => UnitSystem::Metric,
            other => {
                return Err(format!(
                    "Profile '{}' has unsupported units '{}'; expected 'imperial' or 'metric'",
                    self.name, other
                )
                .into())
            }
        };

        // Avoid needless round trips (and preserve every stored bit) for the common same-unit load.
        if source == target {
            return Ok(self);
        }

        self.velocity = UnitConverter::velocity_from_metric(
            UnitConverter::velocity_to_metric(self.velocity, source),
            target,
        );
        self.mass = UnitConverter::mass_from_metric(
            UnitConverter::mass_to_metric(self.mass, source),
            target,
        );
        self.diameter = UnitConverter::diameter_from_metric(
            UnitConverter::diameter_to_metric(self.diameter, source),
            target,
        );
        self.twist_rate = self.twist_rate.map(|value| {
            UnitConverter::diameter_from_metric(
                UnitConverter::diameter_to_metric(value, source),
                target,
            )
        });
        self.sight_height = self.sight_height.map(|value| {
            UnitConverter::sight_height_from_metric(
                UnitConverter::sight_height_to_metric(value, source),
                target,
            )
        });
        self.zero_distance = self.zero_distance.map(|value| {
            UnitConverter::distance_from_metric(
                UnitConverter::distance_to_metric(value, source),
                target,
            )
        });
        self.temperature = UnitConverter::temperature_from_metric(
            UnitConverter::temperature_to_metric(self.temperature, source),
            target,
        );
        self.pressure = UnitConverter::pressure_from_metric(
            UnitConverter::pressure_to_metric(self.pressure, source),
            target,
        );
        self.altitude = UnitConverter::altitude_from_metric(
            UnitConverter::altitude_to_metric(self.altitude, source),
            target,
        );
        self.wind_speed = self.wind_speed.map(|value| {
            UnitConverter::wind_from_metric(UnitConverter::wind_to_metric(value, source), target)
        });
        self.auto_zero = self.auto_zero.map(|value| {
            UnitConverter::distance_from_metric(
                UnitConverter::distance_to_metric(value, source),
                target,
            )
        });
        self.bullet_length = self.bullet_length.map(|value| {
            UnitConverter::diameter_from_metric(
                UnitConverter::diameter_to_metric(value, source),
                target,
            )
        });

        // `bc_segments`, `drag_curve`, `elevation_click`, and `windage_click` are
        // INTENTIONALLY left untouched here — do not "complete" this by adding a
        // conversion pass for them:
        //   * `ProfileBcSegment::velocity_mps` is pinned to SI regardless of `units` (see its
        //     doc comment), so there is nothing to convert.
        //   * `ProfileDragPoint`'s `mach` and `cd` are both dimensionless, so there is nothing
        //     to convert either.
        //   * `elevation_click`/`windage_click` (MBA-1355) are angular turret graduations
        //     (mil/moa/smoa/iphy, chosen by their own suffix), not linear distances, so
        //     imperial/metric has no bearing on them.
        // Converting them would silently rescale physically meaningful data (e.g. treating an
        // m/s breakpoint as if it were fps) rather than fix anything.
        self.units = match target {
            UnitSystem::Imperial => "imperial",
            UnitSystem::Metric => "metric",
        }
        .to_string();

        Ok(self)
    }
}

// ============================================================================
// Fuzzy Column Name Matching (MBA-614)
// ============================================================================

/// Known column names for gun profiles and locations
const KNOWN_COLUMNS: &[&str] = &[
    // Gun profile columns
    "RIFLE_NAME",
    "BULLET_NAME",
    "CALIBER",
    "BULLET_LENGTH",
    "VELOCITY",
    "MUZZLE_VELOCITY",
    "MV",
    "ZERO_TEMP",
    "ZERO_ALT",
    "ZERO_RANGE",
    "ZERO_DISTANCE",
    "ZERO",
    "BC",
    "BC_TYPE",
    "BC_ADJ",
    "BC_ADJUSTMENT",
    "BULLET_WEIGHT",
    "WIND_SPEED",
    "TARGET_SPEED",
    "WIND_ANGLE",
    "SIGHT_HEIGHT",
    "BARREL_TWIST",
    "TWIST_RATE",
    "TWIST",
    "VELOCITY_ADJ",
    "VEL_ADJ",
    "RANGE_MIN",
    "RANGE_MAX",
    "RANGE_INCR",
    "GUN_MPH",
    "COA",
    "POWDER_NAME",
    "POWDER_CHARGE",
    "PLASTIC_TIP_LENGTH",
    "JBM_BULLET_ID",
    "SPIN_DRIFT",
    "CHART_TYPE",
    "REGENERATE",
    "V_OFFSET_MIL",
    "H_OFFSET_MIL",
    "LATITUDE",
    "LONGITUDE",
    "DRAG_MODEL",
    // Location columns
    "LOCATION_NAME",
    "ALTITUDE",
    "ALT",
    "PRESSURE",
    "PRESSURE(HPA OR INHG)",
    "TARGET_TEMP",
    "TEMPERATURE",
    "TEMP",
    "HUMIDITY",
    "DA",
    "DENSITY_ALTITUDE",
    "WIND_DIR",
    "WIND_DIRECTION",
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
    let headers: Vec<String> = reader
        .headers()?
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
        eprintln!(
            "Warning: Column '{}' not recognized - did you mean '{}'? Using {}.",
            original.trim(),
            corrected,
            corrected
        );
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
        return Err(format!(
            "Profile '{}' not found. Use 'ballistics profile list' to see available profiles.",
            name
        )
        .into());
    }
    let content = fs::read_to_string(&path)?;
    let profile: ProfileData = serde_json::from_str(&content)?;
    Ok(profile)
}

/// Load a profile for calculation, converting its stored display units to the active CLI units.
fn load_profile_for_units(name: &str, units: UnitSystem) -> Result<ProfileData, Box<dyn Error>> {
    load_profile(name)?.converted_to(units)
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

// ============================================================================
// .a7p profile import: A7pDocument -> ProfileData mapping + import report
// ============================================================================
//
// Gated on the `profile-import` feature (default-on; off for WASM builds via
// `--no-default-features`) because the signatures below name types from
// `ballistics_engine::profile_import`, which does not exist without the
// feature. Mirrors the `#[cfg(feature = "pdf")]` gating pattern used for the
// PDF dope-card path elsewhere in this file.

/// Everything `profile import` produced: the profile to save plus the honest
/// account of what mapped, what did not, and why.
#[cfg(feature = "profile-import")]
#[derive(Debug)]
struct ImportReport {
    /// (source field, raw value, converted value, destination field)
    mapped: Vec<[String; 4]>,
    /// (source field, human explanation) — data the profile store cannot hold.
    unmapped: Vec<(String, String)>,
    warnings: Vec<String>,
}

#[cfg(feature = "profile-import")]
#[derive(Debug)]
struct A7pImportOutcome {
    profile: ProfileData,
    report: ImportReport,
}

// Re-export of the shared constant under the name this module's a7p-mapping
// call sites already use (MBA-1327: single source of truth for grain<->gram).
#[cfg(feature = "profile-import")]
use ballistics_engine::constants::GRAMS_PER_GRAIN as GRAIN_TO_GRAM;
#[cfg(feature = "profile-import")]
const IN_TO_MM: f64 = 25.4;

/// Restrict imported profile names to characters that are safe as file names
/// in the profile store (`~/.ballistics/profiles/<name>.json`).
#[cfg(feature = "profile-import")]
fn sanitize_profile_name(raw: &str) -> String {
    let cleaned: String = raw
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, ' ' | '.' | '_' | '-') {
                c
            } else {
                '_'
            }
        })
        .collect();
    let trimmed = cleaned.trim().to_string();
    if trimmed.is_empty() {
        "imported-a7p".to_string()
    } else {
        trimmed
    }
}

#[cfg(feature = "profile-import")]
fn map_a7p_to_profile(
    doc: &ballistics_engine::profile_import::A7pDocument,
    name_override: Option<&str>,
) -> Result<A7pImportOutcome, String> {
    use ballistics_engine::profile_import::{A7pBcType, EnvelopeStatus};
    let src = &doc.profile;

    let mut report = ImportReport {
        mapped: Vec::new(),
        unmapped: Vec::new(),
        warnings: Vec::new(),
    };
    if let EnvelopeStatus::Mismatch { expected, actual } = &doc.envelope {
        report.warnings.push(format!(
            "checksum mismatch (file says {expected}, payload hashes to {actual}) — file may be corrupted"
        ));
    }

    let name = match name_override {
        Some(n) => n.to_string(),
        None => sanitize_profile_name(&src.profile_name),
    };

    let mut push = |field: &str, raw: String, converted: String, dest: &str| {
        report
            .mapped
            .push([field.to_string(), raw, converted, dest.to_string()]);
    };
    push(
        "profile_name",
        src.profile_name.clone(),
        name.clone(),
        "name",
    );
    push(
        "c_muzzle_velocity",
        format!("{:.1} m/s", src.muzzle_velocity_mps),
        format!("{:.1} m/s", src.muzzle_velocity_mps),
        "velocity (muzzle velocity)",
    );

    // Resolve drag model + BC-related profile fields. Branches by bc_type because G1/G7 and
    // CUSTOM disagree on what `coef_rows_raw` even means (velocity-BC rows vs Mach-Cd rows —
    // see A7pProfile::bc_rows()/custom_rows()) and on what the scalar `bc` field should hold.
    let (drag_model, bc, bc_segments, drag_curve): (
        &str,
        f64,
        Option<Vec<ProfileBcSegment>>,
        Option<Vec<ProfileDragPoint>>,
    ) = match src.bc_type {
        A7pBcType::G1 | A7pBcType::G7 => {
            let drag_model = if matches!(src.bc_type, A7pBcType::G1) {
                "G1"
            } else {
                "G7"
            };
            let rows = src.bc_rows();
            // The row measured at the highest velocity is the muzzle-regime BC, retained as
            // the scalar `bc` for back-compat with tools that only understand one BC.
            let (bc, bc_row_velocity) = rows.iter().copied().max_by(|a, b| a.1.total_cmp(&b.1)).ok_or_else(
                || "no BC rows in file — cannot build a profile without a BC".to_string(),
            )?;
            push(
                "coef_rows[fastest]",
                format!("BC {bc:.3} @ {bc_row_velocity:.0} m/s"),
                format!("{bc:.3} ({drag_model})"),
                "bc + drag_model",
            );

            let bc_segments = if rows.len() > 1 {
                // Descending by velocity: matches bc_segments_from_profile's/the engine's
                // "fastest row governs the muzzle regime" convention, and puts the back-compat
                // scalar `bc` above (= sorted[0].bc) in visible agreement with this list.
                let mut sorted = rows.clone();
                sorted.sort_by(|a, b| b.1.total_cmp(&a.1));
                push(
                    "coef_rows[all]",
                    format!("{} row(s), fastest {bc:.3} @ {bc_row_velocity:.0} m/s", rows.len()),
                    format!("{} bc_segments (velocity-banded, descending)", sorted.len()),
                    "bc_segments",
                );
                Some(
                    sorted
                        .into_iter()
                        .map(|(bc, velocity_mps)| ProfileBcSegment { bc, velocity_mps })
                        .collect(),
                )
            } else {
                None
            };
            (drag_model, bc, bc_segments, None)
        }
        A7pBcType::Custom => {
            let mut pairs = src.custom_rows(); // (Cd, Mach), file order
            pairs.sort_by(|a, b| a.1.total_cmp(&b.1)); // ascending by Mach (DragTable requirement)
            let mach_values: Vec<f64> = pairs.iter().map(|&(_, mach)| mach).collect();
            let cd_values: Vec<f64> = pairs.iter().map(|&(cd, _)| cd).collect();
            // Validate now (at import time) rather than only at first solve, so a malformed
            // curve is a clear import-time error instead of a confusing later failure.
            ballistics_engine::drag::DragTable::try_new(mach_values.clone(), cd_values.clone())
                .map_err(|e| format!("CUSTOM drag curve is invalid: {e}"))?;
            push(
                "coef_rows[custom]",
                format!("{} Cd/Mach row(s)", pairs.len()),
                format!(
                    "{} drag_curve point(s), Mach {:.3}-{:.3}",
                    pairs.len(),
                    mach_values.first().copied().unwrap_or(0.0),
                    mach_values.last().copied().unwrap_or(0.0)
                ),
                "drag_curve",
            );
            // No scalar BC applies to a full drag curve — see map_a7p_to_profile's module-level
            // rationale in the report below. `bc: 0.0` is an intentionally-invalid sentinel:
            //   * it is physically inert once drag_curve is consumed (custom_drag_table divides
            //     by sectional density, not `bc_value` — BallisticInputs::custom_drag_denominator);
            //   * commands that do not YET consume drag_curve (see CLI_USAGE.md's a7p import
            //     section) will fail loudly (`bc_value must be finite and greater than zero`)
            //     instead of silently running the wrong physics under an assumed G1 model. That
            //     loud failure is the honest outcome for an unwired path, not a bug to paper over.
            report.warnings.push(
                "bc_type CUSTOM: no scalar BC applies to a full drag curve; the profile's 'bc' \
                 field is set to 0.0 as an inert sentinel. It is unused once drag_curve is \
                 consumed (a custom drag table replaces the BC-based retardation model \
                 entirely); commands that do not yet consume drag_curve will fail loudly \
                 (bc_value must be > 0) rather than silently assuming a G1 model — see \
                 CLI_USAGE.md."
                    .to_string(),
            );
            let drag_curve = Some(
                pairs
                    .into_iter()
                    .map(|(cd, mach)| ProfileDragPoint { mach, cd })
                    .collect(),
            );
            ("CUSTOM", 0.0, None, drag_curve)
        }
        A7pBcType::Other(v) => {
            return Err(format!("unknown bc_type {v} — file newer than this importer"))
        }
    };

    push(
        "b_weight",
        format!("{:.1} gr", src.bullet_weight_gr),
        format!("{:.4} g", src.bullet_weight_gr * GRAIN_TO_GRAM),
        "mass",
    );
    push(
        "b_diameter",
        format!("{:.3} in", src.bullet_diameter_in),
        format!("{:.3} mm", src.bullet_diameter_in * IN_TO_MM),
        "diameter",
    );
    push(
        "b_length",
        format!("{:.3} in", src.bullet_length_in),
        format!("{:.2} mm", src.bullet_length_in * IN_TO_MM),
        "bullet_length",
    );
    push(
        "r_twist / twist_dir",
        format!(
            "{:.2} in/turn, {}",
            src.twist_in_per_turn,
            if src.twist_right { "RIGHT" } else { "LEFT" }
        ),
        format!("{:.1} mm/turn", src.twist_in_per_turn * IN_TO_MM),
        "twist_rate + twist_right",
    );
    push(
        "sc_height",
        format!("{:.0} mm", src.sight_height_mm),
        format!("{:.0} mm", src.sight_height_mm),
        "sight_height",
    );
    if let Some(zd) = src.zero_distance_m {
        push(
            "distances[c_zero_distance_idx]",
            format!("{zd:.1} m"),
            format!("{zd:.1} m"),
            "zero_distance + auto_zero",
        );
    }
    push(
        "c_zero_air_temperature",
        format!("{:.1} C", src.air_temperature_c),
        format!("{:.1} C", src.air_temperature_c),
        "temperature",
    );
    push(
        "c_zero_air_pressure",
        format!("{:.1} hPa", src.air_pressure_hpa),
        format!("{:.1} hPa", src.air_pressure_hpa),
        "pressure",
    );
    push(
        "c_zero_air_humidity",
        format!("{:.0} %", src.air_humidity_pct),
        format!("{:.0} %", src.air_humidity_pct),
        "humidity",
    );
    if !src.bullet_name.is_empty() {
        push(
            "bullet_name",
            src.bullet_name.clone(),
            src.bullet_name.clone(),
            "bullet_name",
        );
    }

    // Honest non-mapping: things the profile store cannot hold today.
    let tcoeff_mps_per_c =
        src.muzzle_velocity_mps * (src.temp_coeff_pct_per_15c / 100.0) / 15.0;
    report.unmapped.push((
        "c_t_coeff".to_string(),
        format!(
            "{:.3} %/15C = {:.3} m/s per C powder sensitivity — profile schema does not model \
             powder sensitivity",
            src.temp_coeff_pct_per_15c, tcoeff_mps_per_c
        ),
    ));
    report.unmapped.push((
        "c_zero_p_temperature".to_string(),
        format!("{:.0} C powder temperature at zeroing", src.powder_temperature_c),
    ));
    report.unmapped.push((
        "c_zero_temperature".to_string(),
        format!(
            "{:.0} C ambient temperature when the zero was established",
            src.zero_temperature_c
        ),
    ));
    if src.w_pitch_raw != 0 {
        report.unmapped.push((
            "c_zero_w_pitch".to_string(),
            format!(
                "zeroing pitch value ({}) — sight-pitch state is not modeled",
                src.w_pitch_raw
            ),
        ));
    }
    if src.zero_x_raw != 0 || src.zero_y_raw != 0 {
        report.unmapped.push((
            "zero_x / zero_y".to_string(),
            format!(
                "scope zeroing click offsets ({}, {}) — click state is not modeled",
                src.zero_x_raw, src.zero_y_raw
            ),
        ));
    }
    if !src.distances_m.is_empty() {
        report.unmapped.push((
            "distances".to_string(),
            format!("{} range-card entries (device UI list)", src.distances_m.len()),
        ));
    }
    if src.switches_count > 0 {
        report.unmapped.push((
            "switches".to_string(),
            format!("{} device UI switch entries", src.switches_count),
        ));
    }
    for (field, value) in [
        ("cartridge_name", &src.cartridge_name),
        ("caliber", &src.caliber),
        ("short_name_top", &src.short_name_top),
        ("short_name_bot", &src.short_name_bot),
        ("device_uuid", &src.device_uuid),
    ] {
        if !value.is_empty() {
            report
                .unmapped
                .push((field.to_string(), format!("\"{}\"", value.trim())));
        }
    }
    if !src.user_note.trim().is_empty() {
        report
            .unmapped
            .push(("user_note".to_string(), format!("\"{}\"", src.user_note.trim())));
    }
    for unknown in &doc.unknown_fields {
        report.unmapped.push((
            format!("{} field #{}", unknown.context, unknown.number),
            "unknown field (file newer than this importer)".to_string(),
        ));
    }

    let profile = ProfileData {
        name,
        velocity: src.muzzle_velocity_mps,
        bc,
        mass: src.bullet_weight_gr * GRAIN_TO_GRAM,
        diameter: src.bullet_diameter_in * IN_TO_MM,
        drag_model: drag_model.to_string(),
        twist_rate: Some(src.twist_in_per_turn * IN_TO_MM),
        sight_height: Some(src.sight_height_mm),
        zero_distance: src.zero_distance_m,
        units: "metric".to_string(),
        temperature: src.air_temperature_c,
        pressure: src.air_pressure_hpa,
        humidity: src.air_humidity_pct,
        altitude: 0.0,
        bullet_name: if src.bullet_name.is_empty() {
            None
        } else {
            Some(src.bullet_name.clone())
        },
        created: Some(timestamp_string()),
        wind_speed: None,
        wind_direction: None,
        shooting_angle: None,
        auto_zero: src.zero_distance_m,
        twist_right: Some(src.twist_right),
        // Left unset (not forced to Some(true)): consumers already treat "bc_segments
        // present" as implying velocity-segment behavior is active (the existing
        // `effective_use_bc_segments = use_bc_segments || bc_segments_data.is_some()`
        // pattern), so this boolean keeps its original meaning rather than being
        // overloaded by import.
        use_bc_segments: None,
        bullet_length: Some(src.bullet_length_in * IN_TO_MM),
        // .a7p has no click-graduation concept; the profile leaves these for
        // `profile save --elevation-click/--windage-click` (or a later edit) to fill in.
        elevation_click: None,
        windage_click: None,
        bc_segments,
        drag_curve,
    };

    Ok(A7pImportOutcome { profile, report })
}

#[cfg(feature = "profile-import")]
fn render_import_report(report: &ImportReport) -> String {
    let mut out = String::new();
    out.push_str("Imported fields:\n");
    out.push_str(&format!(
        "  {:<32} {:<26} {:<22} {}\n",
        "SOURCE (.a7p)", "VALUE", "CONVERTED", "DESTINATION"
    ));
    for [field, raw, converted, dest] in &report.mapped {
        out.push_str(&format!(
            "  {field:<32} {raw:<26} {converted:<22} {dest}\n"
        ));
    }
    // muzzle velocity appears in the header row name for the g7 render test
    if !report.unmapped.is_empty() {
        out.push_str("\nNOT imported (no destination in the profile store):\n");
        for (field, why) in &report.unmapped {
            out.push_str(&format!("  {field:<32} {why}\n"));
        }
    }
    if !report.warnings.is_empty() {
        out.push_str("\nWarnings:\n");
        for w in &report.warnings {
            out.push_str(&format!("  ! {w}\n"));
        }
    }
    out
}

/// Convert drop to adjustment unit (MIL, MOA, or SMOA/IPHY). `unit` maps onto the shared
/// `ballistics_engine::adjustment` module's `ClickBase`/`adjustment_factor` (MBA-1355) so
/// the CLI and the WASM terminal share one conversion table.
///
/// `AdjustmentUnit::Clicks` has no fixed factor of its own — clicks depend on a graduation
/// (`ClickValue`) supplied separately — so callers MUST resolve clicks via `clicks_for()`
/// before ever reaching this function with `unit == Clicks`. That arm exists only as a
/// release-safe fallback (falls back to MIL) guarded by a `debug_assert!`.
fn drop_to_adjustment(drop_yd: f64, range_yd: f64, unit: AdjustmentUnit) -> f64 {
    if range_yd < 1.0 {
        return 0.0;
    }
    let factor = match unit {
        AdjustmentUnit::Mil => adjustment_factor(ClickBase::Mil),
        AdjustmentUnit::Moa => adjustment_factor(ClickBase::Moa),
        AdjustmentUnit::Smoa | AdjustmentUnit::Iphy => adjustment_factor(ClickBase::Smoa),
        AdjustmentUnit::Clicks => {
            // Callers resolve clicks via clicks_for() BEFORE display; reaching this
            // arm is a scoping bug. Fall back to MIL so release builds stay sane.
            debug_assert!(false, "Clicks must be resolved via clicks_for()");
            adjustment_factor(ClickBase::Mil)
        }
    };
    (drop_yd / range_yd) * factor
}

/// SMOA/IPHY-per-MIL ratio (3600/1000 = 3.6, exact) for sites that already hold a
/// mil-based value (e.g. `moving_target::LeadSolution::lead_mil`) and need to rescale it
/// to SMOA/IPHY display units without recomputing from raw drop/range (MBA-1355).
fn smoa_per_mil() -> f64 {
    adjustment_factor(ClickBase::Smoa) / adjustment_factor(ClickBase::Mil)
}

/// Resolves the per-axis turret click graduations for `--adjustment-unit clicks` (MBA-1355).
///
/// Precedence per axis: an explicit CLI flag beats the saved profile's field. Windage
/// falls back to the resolved elevation graduation when neither `flag_wind` nor
/// `profile.windage_click` is set (most turrets share one graduation between the two
/// knobs). Elevation MUST resolve from at least one source — clicks output has nowhere
/// else to get a graduation from — so a missing elevation source is the only failure
/// mode, and the error names `--elevation-click-value` so the fix is obvious from the
/// message alone.
///
/// Only `trajectory` and `come-ups` call this; every other command that still accepts
/// `--adjustment-unit clicks` rejects it before reaching here (see
/// `reject_clicks_out_of_scope`).
fn resolve_click_values(
    flag_elev: Option<&str>,
    flag_wind: Option<&str>,
    profile: Option<&ProfileData>,
) -> Result<Option<(ClickValue, ClickValue)>, String> {
    let elev_str = flag_elev
        .map(str::to_string)
        .or_else(|| profile.and_then(|p| p.elevation_click.clone()));
    let Some(elev_str) = elev_str else {
        return Err(
            "--adjustment-unit clicks requires a turret elevation graduation: pass \
             --elevation-click-value <SIZE><UNIT> (e.g. 0.25moa or 0.1mil), or save one on \
             the profile with `profile save --elevation-click`"
                .to_string(),
        );
    };
    let elevation = parse_click_value(&elev_str)?;

    let wind_str = flag_wind
        .map(str::to_string)
        .or_else(|| profile.and_then(|p| p.windage_click.clone()));
    let windage = match wind_str {
        Some(s) => parse_click_value(&s)?,
        None => elevation,
    };

    Ok(Some((elevation, windage)))
}

/// Elevation/windage adjustment display value for one axis (MBA-1355): whole clicks
/// (rounded via `clicks_for`) when a click graduation has been resolved for this axis
/// (`--adjustment-unit clicks`), else the pre-existing `drop_to_adjustment` angular value. Shared by
/// the come-ups table and the trajectory PDF dope-card rows so both format Clicks
/// identically.
fn adjustment_display(
    drop_yd: f64,
    range_yd: f64,
    unit: AdjustmentUnit,
    click: Option<ClickValue>,
) -> f64 {
    match click {
        Some(c) => clicks_for(drop_yd, range_yd, &c) as f64,
        None => drop_to_adjustment(drop_yd, range_yd, unit),
    }
}

/// `--adjustment-unit clicks` real click resolution exists only for `trajectory`/`come-ups`
/// (MBA-1355 Task 2 scope). Every other command that still reads `AdjustmentUnit` (wind
/// card, lead/moving-target, range-table, compare) must reject it up front rather than
/// silently falling back to MIL through the `debug_assert!`-guarded arms further down —
/// call this as the first statement of each such handler.
fn reject_clicks_out_of_scope(unit: AdjustmentUnit) {
    if matches!(unit, AdjustmentUnit::Clicks) {
        eprintln!(
            "error: --adjustment-unit clicks is currently supported for trajectory and come-ups only (MBA-1355)"
        );
        std::process::exit(1);
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

/// Warning for drag-model strings the CLI silently coerces to G1. The native
/// CLI supports G1/G7 only (DragModelArg); anything else — including families
/// the library enum knows (G2/G5/G6/G8/GI/GS) and plain typos — became G1
/// with no feedback until MBA-1386's fix-half added this warning.
fn drag_model_arg_warning(s: &str) -> Option<String> {
    match s.to_uppercase().as_str() {
        "G1" | "G7" => None,
        other => Some(format!(
            "warning: drag model '{other}' is not supported by the CLI (G1/G7 only); using G1. Full-family support is tracked in MBA-1386."
        )),
    }
}

/// Parse the drag model string from a profile
fn parse_drag_model_arg(s: &str) -> DragModelArg {
    if let Some(w) = drag_model_arg_warning(s) {
        eprintln!("{w}");
    }
    match s.to_uppercase().as_str() {
        "G7" => DragModelArg::G7,
        _ => DragModelArg::G1,
    }
}

/// Parse a `--powder-temp-curve` value `"TEMP:VEL,TEMP:VEL,..."` into engine-canonical
/// `(temperature_celsius, muzzle_velocity_m_s)` points, converting from the CLI display
/// units. Points are sorted ascending by temperature; at least two distinct-temperature
/// points are required.
fn parse_powder_temp_curve(s: &str, units: UnitSystem) -> Result<Vec<(f64, f64)>, String> {
    let mut pts: Vec<(f64, f64)> = Vec::new();
    for (i, raw) in s.split(',').enumerate() {
        let part = raw.trim();
        if part.is_empty() {
            continue;
        }
        let (t_str, v_str) = part.split_once(':').ok_or_else(|| {
            format!(
                "--powder-temp-curve point {} '{}' must be TEMP:VELOCITY",
                i + 1,
                part
            )
        })?;
        let t: f64 = t_str
            .trim()
            .parse()
            .map_err(|_| format!("invalid temperature '{}' in --powder-temp-curve", t_str.trim()))?;
        let v: f64 = v_str
            .trim()
            .parse()
            .map_err(|_| format!("invalid velocity '{}' in --powder-temp-curve", v_str.trim()))?;
        if v.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
            return Err(format!(
                "--powder-temp-curve velocity must be positive (got {})",
                v
            ));
        }
        // Convert to SI: absolute temperature F->C and velocity fps->m/s (metric passes through).
        pts.push((
            UnitConverter::temperature_to_metric(t, units),
            UnitConverter::velocity_to_metric(v, units),
        ));
    }
    if pts.len() < 2 {
        return Err(format!(
            "--powder-temp-curve needs at least 2 TEMP:VELOCITY points (got {})",
            pts.len()
        ));
    }
    pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    for w in pts.windows(2) {
        if (w[0].0 - w[1].0).abs() < f64::EPSILON {
            return Err("--powder-temp-curve has duplicate temperatures".to_string());
        }
    }
    Ok(pts)
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

/// Parse a `--bc-segment` value `"VMIN:VMAX:BC"` into a velocity-keyed `BCSegmentData`.
/// VMIN/VMAX are in the CLI display velocity units (fps imperial, m/s metric) and are
/// converted to fps (the engine's segment unit); BC is dimensionless. The BC applies while
/// the bullet's current velocity is in `[VMIN, VMAX)`.
fn parse_bc_segment(s: &str, units: UnitSystem) -> Result<BCSegmentData, String> {
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() != 3 {
        return Err(format!(
            "--bc-segment expects VMIN:VMAX:BC (e.g. 1500:1800:0.243), got '{}'",
            s
        ));
    }
    let vmin: f64 = parts[0]
        .trim()
        .parse()
        .map_err(|_| format!("--bc-segment: invalid VMIN in '{}'", s))?;
    let vmax: f64 = parts[1]
        .trim()
        .parse()
        .map_err(|_| format!("--bc-segment: invalid VMAX in '{}'", s))?;
    let bc: f64 = parts[2]
        .trim()
        .parse()
        .map_err(|_| format!("--bc-segment: invalid BC in '{}'", s))?;
    if vmin.partial_cmp(&vmax) != Some(std::cmp::Ordering::Less) {
        return Err(format!("--bc-segment: VMIN must be < VMAX in '{}'", s));
    }
    if bc <= 0.0 {
        return Err(format!("--bc-segment: BC must be > 0 in '{}'", s));
    }
    // Display velocity -> fps (the unit the solver compares segment bounds against).
    let to_fps = if matches!(units, UnitSystem::Imperial) {
        1.0
    } else {
        3.280_839_895
    };
    Ok(BCSegmentData {
        velocity_min: vmin * to_fps,
        velocity_max: vmax * to_fps,
        bc_value: bc,
    })
}

/// Load a user drag deck for the CLI, exiting with a clear message on any validation error.
fn load_drag_table_or_exit(path: &std::path::Path) -> ballistics_engine::drag::DragTable {
    match ballistics_engine::drag::DragTable::from_file(path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("Error: {e}");
            std::process::exit(1);
        }
    }
}

/// Build an engine [`DragTable`](ballistics_engine::drag::DragTable) from a profile's stored
/// Mach/Cd curve ([`ProfileData::drag_curve`]). Returns a descriptive error (never panics) if
/// the stored points fail `DragTable::try_new`'s validation — `map_a7p_to_profile` already
/// validates on import, so this should only fire for a hand-edited or foreign-tool-written
/// profile JSON.
fn drag_table_from_profile(
    points: &[ProfileDragPoint],
) -> Result<ballistics_engine::drag::DragTable, String> {
    let mach: Vec<f64> = points.iter().map(|p| p.mach).collect();
    let cd: Vec<f64> = points.iter().map(|p| p.cd).collect();
    ballistics_engine::drag::DragTable::try_new(mach, cd)
}

/// Convert a profile's stored velocity-BC breakpoints ([`ProfileData::bc_segments`]) into the
/// engine's velocity-banded [`BCSegmentData`] schedule.
///
/// Mirrors the ArcherBC2 "coefficient row" convention preserved by `A7pProfile::bc_rows()`,
/// which matches standard multi-BC velocity tables (Applied Ballistics/Kestrel-style): each
/// row's OWN velocity is the point *at and above which* that row's BC applies, down to the
/// next (lower-velocity) row's own velocity, which is where the next row's BC takes over.
/// Concretely, after sorting descending by velocity, row `i`'s band is
/// `[breakpoint(i), breakpoint(i-1))` — the fastest row's band is open-ended upward (its own
/// velocity is a floor, not a ceiling: it still governs a muzzle velocity measured higher than
/// its own characterization point) and the SLOWEST row's band is forced open-ended down to
/// zero (overriding its own stated velocity as a lower bound) so the full velocity axis is
/// covered with no gap. This is deliberately non-overlapping and gapless, unlike naively
/// pairing each row with its neighbor's velocity on both sides.
///
/// Velocities convert m/s -> fps (`BCSegmentData`'s unit, matching `parse_bc_segment`). The
/// fastest row's upper bound uses the same 5000 fps top sentinel as
/// `BCSegmentEstimator::estimate_bc_segments` (comfortably above any realistic muzzle
/// velocity), floored to strictly exceed the row's own velocity so an unrealistically fast
/// import still gets a valid (non-empty) top band.
fn bc_segments_from_profile(rows: &[ProfileBcSegment]) -> Vec<BCSegmentData> {
    const MPS_TO_FPS: f64 = 3.280_839_895;
    const TOP_SENTINEL_FPS: f64 = 5000.0;

    let mut sorted: Vec<(f64, f64)> = rows
        .iter()
        .map(|r| (r.bc, r.velocity_mps * MPS_TO_FPS))
        .collect();
    sorted.sort_by(|a, b| b.1.total_cmp(&a.1)); // descending by velocity
    let n = sorted.len();
    sorted
        .iter()
        .enumerate()
        .map(|(i, &(bc, v))| {
            let velocity_max = if i == 0 {
                v.max(TOP_SENTINEL_FPS) + 1.0
            } else {
                sorted[i - 1].1
            };
            let velocity_min = if i + 1 == n { 0.0 } else { v };
            BCSegmentData {
                velocity_min,
                velocity_max,
                bc_value: bc,
            }
        })
        .collect()
}

/// Resolve the velocity-keyed BC schedule once so auto-zero, native flight, and compare flight
/// cannot synthesize different drag inputs. Explicit manual/table segments take precedence over
/// the characteristic-based estimator.
fn resolve_velocity_bc_segments(
    provided: Option<&[BCSegmentData]>,
    use_estimator: bool,
    bc: f64,
    diameter_m: f64,
    mass_kg: f64,
    drag_model: DragModelArg,
) -> Option<Vec<BCSegmentData>> {
    if let Some(segments) = provided {
        return Some(segments.to_vec());
    }
    if !use_estimator {
        return None;
    }

    Some(
        ballistics_engine::bc_estimation::BCSegmentEstimator::estimate_bc_segments(
            bc,
            diameter_m / 0.0254,
            mass_kg / GRAINS_TO_KG,
            "",
            match drag_model {
                DragModelArg::G1 => "G1",
                DragModelArg::G7 => "G7",
            },
        ),
    )
}

fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();

    match cli.command {
        Commands::SolveJson => std::process::exit(solve_json_command::run_stdio()),

        Commands::Mcp => std::process::exit(mcp_command::run_stdio()),

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
            wind_vertical,
            wind_segment,
            bc_segment,
            temperature,
            pressure,
            humidity,
            altitude,
            output,
            full,
            plot,
            auto_zero,
            sight_height,
            bore_height,
            ignore_ground_impact,
            use_bc_segments,
            drag_table,
            cant,
            print_bc_segments,
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
            powder_temp_curve,
            zero_velocity,
            zero_temperature,
            zero_pressure,
            zero_humidity,
            zero_altitude,
            zero_powder_temp,
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
            adjustment_unit,
            elevation_click_value,
            windage_click_value,
        } => {
            // Load profile from CSV if specified
            let profile_data: HashMap<String, String> =
                if let (Some(path), Some(row)) = (&profile, &profile_row) {
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
                load_profile_for_units(name, cli.units).unwrap_or_else(|e| {
                    eprintln!("Error loading saved profile: {}", e);
                    std::process::exit(1);
                })
            });
            if let Some(ref sp) = saved_profile_data {
                eprintln!("Loaded saved profile '{}'", sp.name);
            }

            // MBA-1355: resolve turret click graduations FIRST, eagerly — before any of
            // the Ring column / PDF dope card display work below — so a missing
            // graduation fails fast with a clear flag name, regardless of whether this
            // particular run would even reach a display site that needs it.
            let (elevation_click, windage_click): (Option<ClickValue>, Option<ClickValue>) =
                if matches!(adjustment_unit, AdjustmentUnit::Clicks) {
                    match resolve_click_values(
                        elevation_click_value.as_deref(),
                        windage_click_value.as_deref(),
                        saved_profile_data.as_ref(),
                    ) {
                        Ok(Some((el, wi))) => (Some(el), Some(wi)),
                        Ok(None) => (None, None),
                        Err(e) => {
                            eprintln!("error: {e}");
                            std::process::exit(1);
                        }
                    }
                } else {
                    (None, None)
                };

            // Load location from CSV if specified
            let location_data: HashMap<String, String> =
                if let (Some(path), Some(site_name)) = (&location, &site) {
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
                .or_else(|| {
                    let v = csv_get_f64(&profile_data, &["VELOCITY", "MV", "MUZZLE_VELOCITY"], 0.0);
                    if v > 0.0 {
                        Some(v)
                    } else {
                        None
                    }
                })
                .or_else(|| saved_profile_data.as_ref().map(|p| p.velocity))
                .unwrap_or(2800.0);
            let final_velocity_adj = if velocity_adjustment != 0.0 {
                velocity_adjustment
            } else {
                csv_get_f64(&profile_data, &["VELOCITY_ADJ", "VEL_ADJ"], 0.0)
            };
            let final_bc = bc
                .or_else(|| {
                    let v = csv_get_f64(&profile_data, &["BC"], 0.0);
                    if v > 0.0 {
                        Some(v)
                    } else {
                        None
                    }
                })
                .or_else(|| saved_profile_data.as_ref().map(|p| p.bc))
                .unwrap_or(0.5);
            let final_bc_adj = if bc_adjustment != 1.0 {
                bc_adjustment
            } else {
                csv_get_f64(&profile_data, &["BC_ADJ", "BC_ADJUSTMENT"], 1.0)
            };
            let final_mass = mass
                .or_else(|| saved_profile_data.as_ref().map(|p| p.mass))
                .unwrap_or_else(|| {
                    eprintln!("Error: --mass is required (or use --saved-profile)");
                    std::process::exit(1);
                });
            let final_diameter = diameter
                .or_else(|| saved_profile_data.as_ref().map(|p| p.diameter))
                .unwrap_or_else(|| {
                    eprintln!("Error: --diameter is required (or use --saved-profile)");
                    std::process::exit(1);
                });
            let final_max_range = max_range;
            let final_wind_speed = if wind_speed != 0.0 {
                wind_speed
            } else {
                saved_profile_data
                    .as_ref()
                    .and_then(|p| p.wind_speed)
                    .unwrap_or(0.0)
            };
            let final_wind_direction = if wind_direction != 0.0 {
                wind_direction
            } else {
                csv_get_f64(&location_data, &["WIND_DIR", "WIND_DIRECTION"], 0.0)
            };

            // Location overrides (environmental conditions).
            // Resolve the per-unit standard once so CSV-less runs use the right
            // standard atmosphere for --units (MBA-960/961).
            let std_temperature = UnitConverter::resolve_temperature(None, cli.units)?;
            let std_pressure = UnitConverter::resolve_pressure(None, cli.units)?;
            let final_temperature = match temperature {
                Some(t) => UnitConverter::resolve_temperature(Some(t), cli.units)?,
                None => csv_get_f64(
                    &location_data,
                    &["TARGET_TEMP", "TEMPERATURE", "TEMP"],
                    csv_get_f64(&profile_data, &["ZERO_TEMP"], std_temperature),
                ),
            };
            let final_pressure = match pressure {
                Some(p) => UnitConverter::resolve_pressure(Some(p), cli.units)?,
                None => csv_get_f64(
                    &location_data,
                    &["PRESSURE", "PRESSURE(HPA OR INHG)"],
                    std_pressure,
                ),
            };
            let final_humidity = if humidity != 50.0 {
                humidity
            } else {
                csv_get_f64(&location_data, &["HUMIDITY"], 50.0)
            };
            let final_altitude = if altitude != 0.0 {
                altitude
            } else {
                csv_get_f64(
                    &location_data,
                    &["ALTITUDE", "ALT"],
                    csv_get_f64(&profile_data, &["ZERO_ALT"], 0.0),
                )
            };

            // Get zero range: CLI --auto-zero overrides profile ZERO_RANGE
            let final_auto_zero: Option<f64> = auto_zero.or_else(|| {
                let zero_from_csv =
                    csv_get_f64(&profile_data, &["ZERO_RANGE", "ZERO_DISTANCE", "ZERO"], 0.0);
                if zero_from_csv > 0.0 {
                    Some(zero_from_csv)
                } else {
                    saved_profile_data
                        .as_ref()
                        .and_then(|p| p.auto_zero.or(p.zero_distance))
                }
            });

            // Resolve additional params from saved profile (if not explicitly set via CLI)
            let drag_model = if saved_profile_data.is_some() && velocity.is_none() && bc.is_none() {
                saved_profile_data
                    .as_ref()
                    .map(|p| parse_drag_model_arg(&p.drag_model))
                    .unwrap_or(drag_model)
            } else {
                drag_model
            };
            let use_bc_segments = use_bc_segments
                || saved_profile_data
                    .as_ref()
                    .and_then(|p| p.use_bc_segments)
                    .unwrap_or(false);
            let twist_right = saved_profile_data
                .as_ref()
                .and_then(|p| p.twist_right)
                .unwrap_or(twist_right);
            let shooting_angle = if shooting_angle != 0.0 {
                shooting_angle
            } else {
                saved_profile_data
                    .as_ref()
                    .and_then(|p| p.shooting_angle)
                    .unwrap_or(0.0)
            };
            let bullet_length =
                bullet_length.or_else(|| saved_profile_data.as_ref().and_then(|p| p.bullet_length));
            let sight_height =
                sight_height.or_else(|| saved_profile_data.as_ref().and_then(|p| p.sight_height));
            let twist_rate =
                twist_rate.or_else(|| saved_profile_data.as_ref().and_then(|p| p.twist_rate));
            // --twist-rate is mm/turn in metric and inches/turn in imperial, matching the
            // `stability` subcommand (MBA-970). The engine and all downstream paths
            // (local solver, TrajectoryConfig, the --compare API request) treat twist as
            // inches/turn, so convert once here. Previously metric twist was used raw as
            // inches (~25x off), so e.g. a 254 mm (1:10") twist read as 254"/turn and
            // collapsed stability/spin-drift to ~0.
            let twist_rate = twist_rate.map(|t| match cli.units {
                UnitSystem::Imperial => t,
                UnitSystem::Metric => t / 25.4, // mm/turn -> inches/turn
            });
            // --powder-temp serves two roles (the models are mutually exclusive):
            //   * with --powder-temp-curve: the POWDER temperature to interpolate at.
            //     Only set when the user supplied it; otherwise the curve falls back to
            //     the ambient --temperature (powder assumed at air temperature).
            //   * with --powder-temp-sensitivity (linear): the reference temperature the
            //     stated velocity was measured at (default 70°F / 21°C).
            // powder_curve_temp_c is the curve override in Celsius (Some iff provided);
            // powder_temp is the linear reference in display units (always resolved).
            let powder_curve_temp_c: Option<f64> =
                powder_temp.map(|t| UnitConverter::temperature_to_metric(t, cli.units));
            let powder_temp = powder_temp.unwrap_or(match cli.units {
                UnitSystem::Imperial => DEFAULT_POWDER_REFERENCE_TEMP_F,
                UnitSystem::Metric => DEFAULT_POWDER_REFERENCE_TEMP_C,
            });
            let powder_temp_sensitivity = powder_temp_sensitivity.unwrap_or(match cli.units {
                UnitSystem::Imperial => 1.0,
                UnitSystem::Metric => 0.3048 / (5.0 / 9.0),
            });

            // Parse the optional measured powder-temperature -> velocity curve into
            // canonical SI points. When present it supersedes the linear sensitivity
            // model at solve time (see cli_api::TrajectorySolver::new).
            let powder_temp_curve_si: Option<Vec<(f64, f64)>> = match powder_temp_curve.as_deref()
            {
                Some(s) => Some(parse_powder_temp_curve(s, cli.units)?),
                None => None,
            };

            // Apply truing adjustments
            let trued_velocity = final_velocity + final_velocity_adj;
            let pre_table_bc = final_bc * final_bc_adj;
            let mut trued_bc = pre_table_bc;
            let trued_velocity_fps = match cli.units {
                UnitSystem::Imperial => trued_velocity,
                UnitSystem::Metric => trued_velocity / 0.3048,
            };
            let caliber_in = match cli.units {
                UnitSystem::Imperial => final_diameter,
                UnitSystem::Metric => final_diameter / 25.4,
            };

            // Apply BC correction from table if provided
            // Generate velocity-dependent BC segments from the table for accurate trajectory
            let mut bc_table_segments: Option<Vec<BCSegmentData>> = None;
            let bc_table_correction: Option<f64> = if let Some(table_path) = &bc_table {
                match BcCorrectionTable::load(table_path) {
                    Ok(table) => {
                        // Get bullet length: CLI arg > CSV profile > estimate from diameter
                        let bullet_length_in = bullet_length
                            .map(|l| match cli.units {
                                UnitSystem::Imperial => l,
                                UnitSystem::Metric => l / 25.4,
                            })
                            .or_else(|| {
                                let len =
                                    csv_get_f64(&profile_data, &["BULLET_LENGTH", "LENGTH"], 0.0);
                                if len > 0.0 {
                                    Some(match cli.units {
                                        UnitSystem::Imperial => len,
                                        UnitSystem::Metric => len / 25.4,
                                    })
                                } else {
                                    None
                                }
                            })
                            .unwrap_or_else(|| {
                                // MBA-1135: mass-based length estimate (was a mass-blind
                                // caliber*3.5 heuristic). Length is a real axis of this BC
                                // correction table. Fall back to caliber*3.5 only if mass<=0.
                                let mass_gr = match cli.units {
                                    UnitSystem::Imperial => final_mass,
                                    UnitSystem::Metric => final_mass * 15.4324,
                                };
                                let est_m = ballistics_engine::stability::estimate_bullet_length_m(
                                    caliber_in * 0.0254,
                                    mass_gr * GRAINS_TO_KG,
                                );
                                if est_m > 0.0 {
                                    est_m / 0.0254
                                } else {
                                    caliber_in * 3.5
                                }
                            });

                        // Get bullet mass in grains
                        let mass_grains = match cli.units {
                            UnitSystem::Imperial => final_mass,         // already in grains
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
                            trued_velocity_fps, // Muzzle velocity
                            2700.0,
                            2500.0,
                            2300.0,
                            2100.0,
                            2000.0,
                            1900.0,
                            1800.0,
                            1700.0,
                            1600.0,
                            1500.0,
                            1400.0,
                            1350.0,
                            1300.0,
                            1250.0,
                            1200.0,
                            1150.0,
                            1100.0,
                            1050.0,
                            1000.0,
                            950.0,
                            900.0,
                            850.0,
                            800.0,
                            700.0,
                            600.0,
                            500.0,
                        ];

                        // Filter breakpoints to be below muzzle velocity and sort descending
                        let mut valid_velocities: Vec<f64> = velocity_breakpoints
                            .into_iter()
                            .filter(|&v| v <= trued_velocity_fps && v >= 500.0)
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
                            let segment_bc = trued_bc * correction;

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
                            trued_velocity_fps,
                        );

                        eprintln!(
                            "BC Table: Loaded {} (v{}, {})",
                            table_path.display(),
                            table.version(),
                            table.dimensions_str()
                        );
                        eprintln!(
                            "BC Table: Generated {} velocity-dependent BC segments",
                            segments.len()
                        );
                        eprintln!("BC Table: Muzzle correction={:.4} for BC={:.3} {} {}gr {:.3}\" @ {:.0}fps",
                            muzzle_correction, final_bc, bc_type_str, mass_grains, bullet_length_in, trued_velocity_fps);

                        // Show BC range if segments were generated
                        if !segments.is_empty() {
                            let min_bc = segments
                                .iter()
                                .map(|s| s.bc_value)
                                .fold(f64::INFINITY, f64::min);
                            let max_bc = segments
                                .iter()
                                .map(|s| s.bc_value)
                                .fold(f64::NEG_INFINITY, f64::max);
                            eprintln!(
                                "BC Table: BC range {:.5} - {:.5} across velocity envelope",
                                min_bc, max_bc
                            );
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
                    if let Ok(table) = manager.get_table(caliber_in) {
                        eprintln!(
                            "BC5D Table: Loaded caliber {:.3} (v{}, API {}, {})",
                            table.caliber(),
                            table.version(),
                            table.api_version(),
                            table.dimensions_str()
                        );
                    }

                    let bullet_length_user: Option<f64> = bullet_length
                        .map(|l| match cli.units {
                            UnitSystem::Imperial => l,
                            UnitSystem::Metric => l / 25.4,
                        })
                        .or_else(|| {
                            let csv_len =
                                csv_get_f64(&profile_data, &["BULLET_LENGTH", "LENGTH"], 0.0);
                            if csv_len > 0.0 {
                                Some(match cli.units {
                                    UnitSystem::Imperial => csv_len,
                                    UnitSystem::Metric => csv_len / 25.4,
                                })
                            } else {
                                None
                            }
                        });
                    let length_is_user = bullet_length_user.is_some();
                    let bullet_length_in = bullet_length_user.unwrap_or_else(|| {
                        // MBA-1135: mass-based length estimate (was a mass-blind caliber*3.5
                        // heuristic). Informational for the v2 BC5D table (length is not a
                        // lookup axis); fall back to caliber*3.5 only if mass<=0.
                        let est_m = ballistics_engine::stability::estimate_bullet_length_m(
                            caliber_in * 0.0254,
                            mass_grains * GRAINS_TO_KG,
                        );
                        if est_m > 0.0 {
                            est_m / 0.0254
                        } else {
                            caliber_in * 3.5
                        }
                    });

                    // BC5D segment generation is keyed by the manually adjusted BC. Preserve
                    // that pre-table value so the muzzle lookup samples the same BC axis.
                    let table_base_bc = trued_bc;

                    if let Some(segments) = generate_bc5d_segments(
                        &mut manager,
                        caliber_in,
                        table_base_bc,
                        bc_type_str,
                        mass_grains,
                        Some(trued_velocity_fps),
                        Some(bullet_length_in),
                        length_is_user,
                        if print_bc_segments { Some(cli.units) } else { None },
                    ) {
                        bc_table_segments = Some(segments);

                        let muzzle_correction = manager
                            .lookup(
                                caliber_in,
                                mass_grains,
                                table_base_bc,
                                trued_velocity_fps,
                                trued_velocity_fps,
                                bc_type_str,
                            )
                            .unwrap_or(1.0);
                        eprintln!(
                            "BC5D Table: Muzzle correction={:.4} for BC={:.3} {} {}gr @ {:.0}fps",
                            muzzle_correction,
                            table_base_bc,
                            bc_type_str,
                            mass_grains,
                            trued_velocity_fps
                        );
                        trued_bc = table_base_bc * muzzle_correction;
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
            if !profile_data.is_empty()
                || !location_data.is_empty()
                || saved_profile_data.is_some()
                || combined_bc_correction.is_some()
            {
                let bc_info = if let Some(corr) = combined_bc_correction {
                    format!(
                        "BC={:.3} (table-corrected={:.4}, factor={:.4})",
                        final_bc, trued_bc, corr
                    )
                } else {
                    format!("BC={:.3} (trued={:.4})", final_bc, trued_bc)
                };
                eprintln!(
                    "Effective values: velocity={:.1} (trued={:.1}), {}",
                    final_velocity, trued_velocity, bc_info
                );
                if !location_data.is_empty() {
                    eprintln!(
                        "                  altitude={:.0}, temp={:.1}, pressure={:.2}",
                        final_altitude, final_temperature, final_pressure
                    );
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
            let wind_vertical_metric = UnitConverter::wind_to_metric(wind_vertical, cli.units);
            let temperature_metric =
                UnitConverter::temperature_to_metric(final_temperature, cli.units);
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

            // Bore height above ground. MBA-1339 unified this flag onto the same inches/mm
            // units as --sight-height and the WASM --muzzle-height flag (previously feet/
            // meters). Defaults are unchanged: 60 in = 5 ft = 1.524 m; 1500 mm = 1.5 m.
            let bore_height_default = match cli.units {
                UnitSystem::Imperial => 60.0,  // inches
                UnitSystem::Metric => 1500.0,  // mm
            };
            let bore_height_value = bore_height.unwrap_or(bore_height_default);
            let bore_height_metric = match cli.units {
                UnitSystem::Imperial => bore_height_value * 0.0254, // inches to meters
                UnitSystem::Metric => bore_height_value * 0.001,    // mm to meters
            };
            // Trajectory altitude feeds local air density, so an absurd bore height
            // silently thins the air over the whole flight (a 25 km "muzzle" flies in
            // ~2% density) or blows up the integrator. Warn — never clamp.
            if bore_height_metric > 1000.0 {
                eprintln!(
                    "warning: --bore-height puts the bore {bore_height_metric:.0} m above ground; \
                     air density is computed at trajectory altitude, so this thins the air for \
                     the entire flight. If you meant site elevation use --altitude; to disable \
                     ground-impact truncation use --ignore-ground-impact."
                );
            }

            // Construct PDF metadata if PDF output is requested
            let pdf_metadata = if matches!(output, OutputFormat::Pdf) {
                // Get rifle name from profile_row or use a default
                let rifle_name = profile_row.clone().unwrap_or_else(|| "Rifle".to_string());
                // Get location name: CLI > site > default
                let loc_name = location_name
                    .clone()
                    .or_else(|| site.clone())
                    .unwrap_or_else(|| "Field".to_string());
                // Get powder and bullet from CLI, or from profile if available
                let powder_name = powder
                    .clone()
                    .or_else(|| profile_data.get("POWDER_NAME").cloned())
                    .unwrap_or_default();
                let bullet_display = bullet_name
                    .clone()
                    .or_else(|| profile_data.get("BULLET_NAME").cloned())
                    .unwrap_or_default();

                // Resolve effective font scale: --font-scale overrides --font-preset
                #[cfg(feature = "pdf")]
                let effective_font_scale = if font_scale != 1.0 {
                    font_scale
                } else if let Some(ref preset) = font_preset {
                    FontSizePreset::from_str(preset)
                        .map(|p| p.scale())
                        .unwrap_or_else(|| {
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
                    // MBA-1325: was a bare passthrough (always mph, ignoring --units — the
                    // only field in this block that didn't follow the sibling unit-aware
                    // conversions below). --target-speed follows the documented mph
                    // imperial / m/s metric convention like every other speed flag, so
                    // metric callers now get their m/s value converted here instead of it
                    // being silently misread as mph.
                    target_speed_mph: match cli.units {
                        UnitSystem::Imperial => target_speed,
                        UnitSystem::Metric => target_speed / 0.44704,
                    },
                    output_file: output_file.clone(),
                    velocity_fps: match cli.units {
                        UnitSystem::Imperial => trued_velocity,
                        UnitSystem::Metric => trued_velocity / 0.3048,
                    },
                    temperature_f: match cli.units {
                        UnitSystem::Imperial => final_temperature,
                        UnitSystem::Metric => final_temperature * 9.0 / 5.0 + 32.0,
                    },
                    pressure_inhg: match cli.units {
                        UnitSystem::Imperial => final_pressure,
                        UnitSystem::Metric => final_pressure / 33.8639,
                    },
                    altitude_ft: match cli.units {
                        UnitSystem::Imperial => final_altitude,
                        UnitSystem::Metric => final_altitude / 0.3048,
                    },
                    wind_speed_mph: match cli.units {
                        UnitSystem::Imperial => final_wind_speed,
                        UnitSystem::Metric => final_wind_speed / 0.44704,
                    },
                    weight_gr: match cli.units {
                        UnitSystem::Imperial => bullet_mass,
                        UnitSystem::Metric => bullet_mass * 15.4324,
                    },
                    font_scale: effective_font_scale,
                    bold_data,
                    adjustment_unit,
                })
            } else {
                None
            };

            // Manual velocity:BC segments (--bc-segment "VMIN:VMAX:BC", display units
            // -> fps) take precedence over any --bc-table-dir / BC5D segments resolved
            // above. Resolved BEFORE the auto-zero solve below so the zero angle is
            // computed with the SAME velocity-keyed BC the trajectory uses (otherwise a
            // segment that changes early-flight drag would mis-zero the shot); via
            // bc_table_segments.is_some() this also implies --use-bc-segments.
            let manual_bc_segments: Vec<BCSegmentData> = bc_segment
                .iter()
                .map(|s| parse_bc_segment(s, cli.units))
                .collect::<Result<Vec<_>, String>>()?;
            if !manual_bc_segments.is_empty() {
                bc_table_segments = Some(manual_bc_segments.clone());
                // Manual segments fully override either correction-table format. Restore
                // the manually adjusted, pre-table BC directly so interior gaps use the
                // same scalar fallback with or without a table, keeping parity with WASM.
                trued_bc = pre_table_bc;
            }

            // Saved-profile velocity-BC breakpoints (MBA-1323 Phase 2: `.a7p` multi-row G1/G7
            // import) fill the schedule only when nothing more specific for THIS run already
            // did (no --bc-segment, no --bc-table-dir correction table) — a saved default,
            // not an override of explicit CLI intent.
            if bc_table_segments.is_none() {
                if let Some(rows) = saved_profile_data.as_ref().and_then(|p| p.bc_segments.as_ref())
                {
                    bc_table_segments = Some(bc_segments_from_profile(rows));
                }
            }

            // Resolve the final velocity-keyed BC schedule before auto-zero. This is the single
            // schedule cloned into the zero solve, native flight, and compare-local flight.
            let bc_segments_data = resolve_velocity_bc_segments(
                bc_table_segments.as_deref(),
                use_bc_segments,
                trued_bc,
                diameter_metric,
                mass_metric,
                drag_model,
            );
            let effective_use_bc_segments = use_bc_segments || bc_segments_data.is_some();

            // Resolve --drag-table once, shared by the zero solve, native flight, and
            // compare-local flight, same as the BC schedule above. A saved profile's Mach/Cd
            // drag curve (MBA-1323 Phase 2: `.a7p` CUSTOM import) is the fallback when
            // --drag-table was not given for THIS run — same "saved default, not an override
            // of explicit CLI intent" precedence as the BC segments above.
            let custom_drag_table = drag_table.as_deref().map(load_drag_table_or_exit).or_else(|| {
                saved_profile_data
                    .as_ref()
                    .and_then(|p| p.drag_curve.as_ref())
                    .map(|pts| {
                        drag_table_from_profile(pts).unwrap_or_else(|e| {
                            eprintln!("Error: saved profile's drag curve is invalid: {e}");
                            std::process::exit(1);
                        })
                    })
            });
            if custom_drag_table.is_some() && effective_use_bc_segments {
                eprintln!(
                    "Warning: --drag-table and BC segments were both provided; the drag table takes \
                     precedence and BC segments are ignored."
                );
            }

            // Calculate zero angle if auto-zero is specified (from CLI or profile)
            let muzzle_angle = if let Some(zero_distance) = final_auto_zero {
                let zero_distance_metric =
                    UnitConverter::distance_to_metric(zero_distance, cli.units);

                // Resolve the condition set the zero ANGLE is solved under: shot-day values
                // overridden by any --zero-* flags supplied. This models a rifle zeroed on a
                // different day (different air density) and/or at a different muzzle velocity
                // (e.g. from a powder-temp/velocity table), while the trajectory below still
                // runs under the current shot-day conditions. Omitting all --zero-* flags
                // falls through to the shot-day values, reproducing the prior behavior.
                let zero_velocity_metric = match zero_velocity {
                    Some(v) => UnitConverter::velocity_to_metric(v, cli.units),
                    None => velocity_metric,
                };
                let zero_temperature_metric = match zero_temperature {
                    Some(t) => UnitConverter::temperature_to_metric(t, cli.units),
                    None => temperature_metric,
                };
                let zero_pressure_metric = match zero_pressure {
                    Some(p) => UnitConverter::pressure_to_metric(p, cli.units),
                    None => pressure_metric,
                };
                let zero_humidity_value = zero_humidity.unwrap_or(final_humidity);
                let zero_altitude_metric = match zero_altitude {
                    Some(a) => UnitConverter::altitude_to_metric(a, cli.units),
                    None => altitude_metric,
                };
                if zero_velocity.is_some()
                    || zero_temperature.is_some()
                    || zero_pressure.is_some()
                    || zero_humidity.is_some()
                    || zero_altitude.is_some()
                {
                    eprintln!(
                        "Solving zero angle under supplied zero-day conditions (velocity/atmosphere \
                         may differ from the shot-day trajectory)."
                    );
                }

                // An explicit zero velocity is absolute. Otherwise let the solver apply the
                // linear powder model at the zero-day air temperature; a curve, when present,
                // retains precedence inside TrajectorySolver::new.
                let apply_zero_linear = zero_velocity.is_none() && use_powder_sensitivity;

                // Create inputs for zero calculation
                let zero_inputs = BallisticInputs {
                    muzzle_velocity: zero_velocity_metric,
                    bc_value: trued_bc,
                    bc_type: match drag_model {
                        DragModelArg::G1 => DragModel::G1,
                        DragModelArg::G7 => DragModel::G7,
                    },
                    bullet_mass: mass_metric,
                    bullet_diameter: diameter_metric,
                    caliber_inches: diameter_metric / 0.0254,
                    weight_grains: mass_metric / GRAINS_TO_KG,
                    bullet_length: bullet_length
                        .map(|l| match cli.units {
                            UnitSystem::Imperial => l * 0.0254,
                            UnitSystem::Metric => l * 0.001,
                        })
                        .unwrap_or_else(|| fallback_bullet_length_m(diameter_metric, mass_metric)),
                    sight_height: sight_height_metric,
                    muzzle_height: bore_height_metric,
                    ground_threshold: 0.0,
                    altitude: zero_altitude_metric,
                    temperature: zero_temperature_metric,
                    pressure: zero_pressure_metric,
                    humidity: zero_humidity_value,
                    use_powder_sensitivity: apply_zero_linear,
                    powder_temp_sensitivity: if apply_zero_linear {
                        UnitConverter::velocity_to_metric(powder_temp_sensitivity, cli.units)
                            / UnitConverter::temperature_delta_to_metric(1.0, cli.units)
                    } else {
                        0.0
                    },
                    // For the linear model this is the reference temperature at which the
                    // stated velocity was measured; --zero-powder-temp remains curve-only.
                    powder_temp: UnitConverter::temperature_to_metric(powder_temp, cli.units),
                    // Let a powder curve set the zero-day velocity at the zero-day
                    // temperature — UNLESS an explicit --zero-velocity was given, which
                    // takes precedence (the constructor would otherwise override it).
                    powder_temp_curve: if zero_velocity.is_some() {
                        None
                    } else {
                        powder_temp_curve_si.clone()
                    },
                    // Explicit zero-day powder temperature wins; otherwise an explicit zero-day
                    // air temperature retains the existing curve-at-zero-air behavior. With
                    // neither override, inherit the shot-day powder lookup temperature.
                    powder_curve_temp_c: if let Some(t) = zero_powder_temp {
                        Some(UnitConverter::temperature_to_metric(t, cli.units))
                    } else if zero_temperature.is_some() {
                        None
                    } else {
                        powder_curve_temp_c
                    },
                    // Zero with the exact BC configuration the subsequent flight uses. The
                    // resolved schedule covers manual/BC5D data and estimator-generated data;
                    // cluster correction remains a scalar-BC fallback under MBA-1175.
                    use_bc_segments: effective_use_bc_segments,
                    bc_segments_data: bc_segments_data.clone(),
                    use_cluster_bc,
                    custom_drag_table: custom_drag_table.clone(),
                    ..Default::default()
                };

                let zero_wind = WindConditions {
                    speed: wind_speed_metric,
                    direction: final_wind_direction.to_radians(),
                    vertical_speed: wind_vertical_metric,
                };
                let zero_atmosphere = AtmosphericConditions {
                    temperature: zero_temperature_metric,
                    pressure: zero_pressure_metric,
                    humidity: zero_humidity_value,
                    altitude: zero_altitude_metric,
                };

                // Target height is the line of sight's ground-referenced height.
                let zero_angle = ballistics_engine::calculate_zero_angle_with_conditions(
                    zero_inputs,
                    zero_distance_metric,
                    bore_height_metric + sight_height_metric,
                    zero_wind,
                    zero_atmosphere,
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
                bullet_length: bullet_length
                    .map(|l| match cli.units {
                        UnitSystem::Imperial => l * 0.0254,
                        UnitSystem::Metric => l * 0.001,
                    })
                    .unwrap_or_else(|| fallback_bullet_length_m(diameter_metric, mass_metric)),
                drag_model,
                max_range: max_range_metric,
                time_step,
                temperature: temperature_metric,
                pressure: pressure_metric,
                humidity: final_humidity,
                altitude: altitude_metric,
                wind_speed: wind_speed_metric,
                wind_direction: final_wind_direction,
                wind_vertical: wind_vertical_metric,
                wind_segments,
                output,
                full,
                plot,
                units: cli.units,
                sight_height: sight_height_metric,
                bore_height: bore_height_metric,
                ignore_ground_impact,
                cant,
                use_bc_segments: effective_use_bc_segments,
                use_cluster_bc,
                bc_segments_data: bc_segments_data.clone(),
                custom_drag_table: custom_drag_table.clone(),
                enable_magnus,
                enable_coriolis,
                enable_spin_drift,
                enable_aerodynamic_jump,
                enable_wind_shear,
                wind_shear_model,
                enable_pitch_damping,
                enable_precession,
                sample_trajectory,
                sample_interval,
                use_rk4: !use_euler,
                use_rk45: !use_rk4_fixed,
                twist_rate,
                twist_right,
                latitude,
                shot_direction,
                shooting_angle,
                use_powder_sensitivity,
                powder_temp_sensitivity,
                powder_temp,
                powder_temp_curve: powder_temp_curve_si.clone(),
                powder_curve_temp_c,
                target_speed,
                adjustment_unit,
                elevation_click,
                windage_click,
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
                            return Err(
                                "Cannot use --online without accepting Terms of Service.".into()
                            );
                        }
                        Err(e) => {
                            eprintln!("Warning: Could not verify TOS acceptance: {}", e);
                            eprintln!("Proceeding with online mode...");
                        }
                    }

                    // Build API request
                    // With no --auto-zero, use 100 in the active CLI distance units.
                    // Convert that fallback before the API request so imperial means
                    // 100 yd (91.44 m), not the previous 100 m / 109.4 yd (MBA-1158).
                    let zero_range_metric = final_auto_zero
                        .map(|d| UnitConverter::distance_to_metric(d, cli.units))
                        .unwrap_or_else(|| UnitConverter::distance_to_metric(100.0, cli.units));

                    let api_request = TrajectoryRequestBuilder::new()
                        .bc_value(trued_bc)
                        .bc_type(match drag_model {
                            DragModelArg::G1 => "G1",
                            DragModelArg::G7 => "G7",
                        })
                        .bullet_mass(mass_metric * 1000.0) // kg to grams
                        .muzzle_velocity(velocity_metric)
                        .target_distance(max_range_metric)
                        .zero_range(zero_range_metric)
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
                            0.0
                        })
                        .enable_weather_zones(enable_weather_zones)
                        .enable_3d_weather(enable_3d_weather)
                        .wind_shear_model(wind_shear_model.as_engine_str())
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
                    // `twist_rate` was already normalized to inches/turn above (metric mm
                    // converted by /25.4, MBA-970), and api_client / the local solver both
                    // consume inches/turn — so forward it verbatim. Both the local and API
                    // legs of --compare see the same converted value, so they stay in
                    // agreement.
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
                                let local_wind_direction_rad = final_wind_direction.to_radians();

                                let mut local_inputs = BallisticInputs {
                                    bc_value: trued_bc,
                                    bc_type: drag_model_enum,
                                    bullet_mass: mass_metric,
                                    muzzle_velocity: velocity_metric,
                                    bullet_diameter: diameter_metric,
                                    bullet_length: bullet_length
                                        .map(|l| match cli.units {
                                            UnitSystem::Imperial => l * 0.0254,
                                            UnitSystem::Metric => l * 0.001,
                                        })
                                        .unwrap_or_else(|| fallback_bullet_length_m(diameter_metric, mass_metric)),
                                    muzzle_angle: muzzle_angle.to_radians(),
                                    target_distance: max_range_metric,
                                    azimuth_angle: 0.0,
                                    shot_azimuth: shot_direction.map(|d| d.to_radians()).unwrap_or(0.0),
                                    shooting_angle: shooting_angle.to_radians(),
                                    cant_angle: cant.to_radians(),
                                    sight_height: sight_height_metric,
                                    muzzle_height: bore_height_metric,
                                    target_height: 0.0,
                                    ground_threshold: if ignore_ground_impact { f64::NEG_INFINITY } else { 0.0 },
                                    altitude: altitude_metric,
                                    temperature: temperature_metric,
                                    pressure: pressure_metric,
                                    humidity: final_humidity,
                                    latitude,
                                    wind_speed: wind_speed_metric,
                                    wind_angle: local_wind_direction_rad,
                                    // MBA-1135: caliber/weight-aware default twist (Miller-inverse)
                                    // instead of a fixed 1:12" when the shooter omits --twist-rate.
                                    twist_rate: twist_rate.unwrap_or_else(|| {
                                        ballistics_engine::stability::default_twist_inches(
                                            diameter_metric,
                                            mass_metric,
                                            velocity_metric,
                                        )
                                    }),
                                    is_twist_right: twist_right,
                                    caliber_inches: diameter_metric / 0.0254,
                                    weight_grains: mass_metric / GRAINS_TO_KG,
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
                                        // Per-degree DELTA conversion, not absolute point (MBA-963).
                                        UnitConverter::velocity_to_metric(powder_temp_sensitivity, cli.units)
                                            / UnitConverter::temperature_delta_to_metric(1.0, cli.units)
                                    } else { 0.0 },
                                    powder_temp: UnitConverter::temperature_to_metric(powder_temp, cli.units),
                                    powder_temp_curve: powder_temp_curve_si.clone(),
                                    powder_curve_temp_c,
                                    tipoff_yaw: 0.0,
                                    tipoff_decay_distance: 50.0,
                                    use_bc_segments: effective_use_bc_segments,
                                    bc_segments: None,
                                    bc_segments_data: bc_segments_data.clone(),
                                    use_enhanced_spin_drift: enable_spin_drift,
                                    enable_aerodynamic_jump,
                                    use_form_factor: false,
                                    enable_wind_shear,
                                    wind_shear_model: if enable_wind_shear { wind_shear_model.as_engine_str().to_string() } else { "none".to_string() },
                                    enable_trajectory_sampling: sample_trajectory,
                                    sample_interval,
                                    enable_pitch_damping,
                                    enable_precession_nutation: enable_precession,
                                    use_cluster_bc,
                                    bc_type_str: None,
                                    custom_drag_table: custom_drag_table.clone(),
                                };

                                let local_wind = WindConditions {
                                    speed: wind_speed_metric,
                                    direction: local_wind_direction_rad,
                                    vertical_speed: wind_vertical_metric,
                                };
                                let local_atmo = AtmosphericConditions {
                                    temperature: temperature_metric,
                                    pressure: pressure_metric,
                                    humidity: final_humidity,
                                    altitude: altitude_metric,
                                };

                                // The API always zeroes the shot at `zero_range_metric`. If the
                                // user did not supply --auto-zero, zero the local comparison leg
                                // at that same default range instead of firing it at the raw
                                // --angle while comparing it to a zeroed API trajectory.
                                if final_auto_zero.is_none() {
                                    local_inputs.muzzle_angle =
                                        ballistics_engine::calculate_zero_angle_with_conditions(
                                            local_inputs.clone(),
                                            zero_range_metric,
                                            bore_height_metric + sight_height_metric,
                                            local_wind.clone(),
                                            local_atmo.clone(),
                                        )?;
                                }

                                let mut local_solver =
                                    TrajectorySolver::new(local_inputs, local_wind, local_atmo);
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
                                println!(
                                    "║ BC Confidence:     {:>8.1}%           ║",
                                    confidence * 100.0
                                );
                                println!("╠════════════════════════════════════════╣");
                            }

                            let (dist_unit, _drop_unit, _vel_unit) = match cli.units {
                                UnitSystem::Metric => ("m", "m", "m/s"),
                                UnitSystem::Imperial => ("yd", "in", "fps"),
                            };

                            println!("║ Range {} │ API Drop │Local Drop│  Δ Drop  ║", dist_unit);
                            println!("╠══════════╪══════════╪══════════╪══════════╣");

                            for point in api_response.trajectory.iter().take(10) {
                                let range_display =
                                    UnitConverter::distance_from_metric(point.range, cli.units);
                                let api_drop_display = if cli.units == UnitSystem::Imperial {
                                    point.drop * 39.3701
                                } else {
                                    point.drop
                                };

                                if let Some(ref local_res) = local_result {
                                    if let Some(pos) = local_res.position_at_range(point.range) {
                                        // API drop is signed height relative to the horizontal
                                        // line of sight (negative below). Local `pos.y` is an
                                        // absolute height above ground, so remove both bore height
                                        // and the sight-over-bore offset before comparing.
                                        let local_drop =
                                            pos.y - (bore_height_metric + sight_height_metric);
                                        let local_drop_display =
                                            if cli.units == UnitSystem::Imperial {
                                                local_drop * 39.3701
                                            } else {
                                                local_drop
                                            };
                                        let delta = api_drop_display - local_drop_display;
                                        println!(
                                            "║ {:>8.1} │ {:>8.2} │ {:>8.2} │ {:>+8.2} ║",
                                            range_display,
                                            api_drop_display,
                                            local_drop_display,
                                            delta
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
            wind_direction_std,
            wind_speed,
            wind_direction,
            wind_vertical,
            target_distance,
            target_radius,
            drag_table,
            cant,
            wez,
            target_size,
            wind_call_error,
            wez_start,
            wez_end,
            wez_step,
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
            let wind_vertical_metric = UnitConverter::wind_to_metric(wind_vertical, cli.units);
            let target_distance_metric =
                target_distance.map(|d| UnitConverter::distance_to_metric(d, cli.units));
            // Hit-zone radius in meters: default 0.3 m, or the user's --target-radius (in the
            // selected unit). MBA-971.
            let target_radius_metric = target_radius
                .map(|r| UnitConverter::distance_to_metric(r, cli.units))
                .unwrap_or(ballistics_engine::DEFAULT_HIT_RADIUS_M);
            // Resolve --drag-table; no bc-segments flag exists on this subcommand, so there is
            // nothing to conflict-warn against (mirrors the trajectory resolve at line ~3338).
            let custom_drag_table = drag_table.as_deref().map(load_drag_table_or_exit);

            if wez {
                // MBA-1317: WEZ (Weapon Employment Zone) sweep mode.
                let target_size_spec = target_size
                    .as_deref()
                    .ok_or("--target-size is required with --wez (e.g. --target-size 18x30)")?;
                let target_size_parsed = parse_target_size(target_size_spec)
                    .map_err(|e| format!("--target-size: {e}"))?;
                let target_size_metric = target_size_parsed.to_metric(cli.units);
                let wind_call_error_metric =
                    UnitConverter::wind_to_metric(wind_call_error, cli.units);
                let wez_start_metric = UnitConverter::distance_to_metric(wez_start, cli.units);
                let wez_end_metric = UnitConverter::distance_to_metric(wez_end, cli.units);
                let wez_step_metric = UnitConverter::distance_to_metric(wez_step, cli.units);

                run_monte_carlo_wez(
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
                    wind_direction_std,
                    wind_speed_metric,
                    wind_direction,
                    wind_vertical_metric,
                    wind_call_error_metric,
                    target_size_metric,
                    wez_start_metric,
                    wez_end_metric,
                    wez_step_metric,
                    custom_drag_table,
                    cant,
                    output,
                    cli.units,
                )?;
            } else {
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
                    wind_direction_std,
                    wind_speed_metric,
                    wind_direction,
                    wind_vertical_metric,
                    target_distance_metric,
                    target_radius_metric,
                    custom_drag_table,
                    cant,
                    output,
                )?;
            }
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
            drag_table,
            output,
        } => {
            let temperature = UnitConverter::resolve_temperature(temperature, cli.units)?;
            let pressure = UnitConverter::resolve_pressure(pressure, cli.units)?;
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
            // Resolve --drag-table; no bc-segments flag exists on this subcommand, so there is
            // nothing to conflict-warn against (mirrors the trajectory resolve at line ~3338).
            let custom_drag_table = drag_table.as_deref().map(load_drag_table_or_exit);

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
                custom_drag_table,
                output,
                cli.units,
            )?;
        }

        Commands::EstimateBC {
            velocity,
            mass,
            diameter,
            data,
            velocity_data,
            drag_model,
            zero_range,
            sight_height,
            temperature,
            pressure,
            humidity,
            altitude,
            distance1,
            drop1,
            distance2,
            drop2,
            output,
        } => {
            // Convert scalar inputs to metric (MBA-716).
            let velocity_metric = UnitConverter::velocity_to_metric(velocity, cli.units);
            let mass_metric = UnitConverter::mass_to_metric(mass, cli.units);
            let diameter_metric = UnitConverter::diameter_to_metric(diameter, cli.units);

            // Atmosphere the data was measured at — BC is only meaningful relative to air
            // density, so this must match the dope card's conditions (defaults = ICAO std).
            let atmosphere = AtmosphericConditions {
                temperature: temperature
                    .map(|t| UnitConverter::temperature_to_metric(t, cli.units))
                    .unwrap_or(15.0),
                pressure: pressure
                    .map(|p| UnitConverter::pressure_to_metric(p, cli.units))
                    .unwrap_or(1013.25),
                humidity,
                altitude: match cli.units {
                    UnitSystem::Imperial => altitude * 0.3048,
                    UnitSystem::Metric => altitude,
                },
            };
            // Zero range (dope-card frame) and sight height, in meters.
            let zero_range_m = zero_range.map(|z| UnitConverter::distance_to_metric(z, cli.units));
            let sight_height_m = sight_height
                .map(|s| match cli.units {
                    UnitSystem::Imperial => s * 0.0254,
                    UnitSystem::Metric => s / 1000.0,
                })
                .unwrap_or(0.05);

            // Drop series: from --data (n-point) or the legacy --distance1/--drop1/... pair.
            let drop_raw: Vec<(f64, f64)> = if let Some(s) = &data {
                parse_data_pairs(s)?
            } else if let (Some(d1), Some(dr1), Some(d2), Some(dr2)) =
                (distance1, drop1, distance2, drop2)
            {
                vec![(d1, dr1), (d2, dr2)]
            } else {
                Vec::new()
            };
            // Velocity series: from --velocity-data (n-point).
            let vel_raw: Vec<(f64, f64)> = match &velocity_data {
                Some(s) => parse_data_pairs(s)?,
                None => Vec::new(),
            };

            if drop_raw.is_empty() && vel_raw.is_empty() {
                return Err("No data provided. Supply drop data (--data \"d,drop;...\" or \
                     --distance1/--drop1/--distance2/--drop2) and/or velocity data \
                     (--velocity-data \"d,vel;...\")."
                    .into());
            }

            // Convert both series to metric (distance -> m, drop in -> m, velocity fps -> m/s).
            let drop_metric: Vec<(f64, f64)> = drop_raw
                .iter()
                .map(|(d, drop)| {
                    (
                        UnitConverter::distance_to_metric(*d, cli.units),
                        match cli.units {
                            UnitSystem::Imperial => *drop * 0.0254,
                            UnitSystem::Metric => *drop,
                        },
                    )
                })
                .collect();
            let vel_metric: Vec<(f64, f64)> = vel_raw
                .iter()
                .map(|(d, v)| {
                    (
                        UnitConverter::distance_to_metric(*d, cli.units),
                        UnitConverter::velocity_to_metric(*v, cli.units),
                    )
                })
                .collect();

            let models = parse_drag_models(&drag_model)?;

            // Guard the most common mistake: feeding a zeroed dope card (a point with ~0
            // drop) without --zero-range, which makes the fit treat it as bore-referenced
            // and returns a wrong (often maxed-out) BC.
            if zero_range_m.is_none()
                && drop_metric.iter().any(|(_, dr)| dr.abs() < 0.05)
                && drop_metric.iter().any(|(_, dr)| dr.abs() > 0.25)
            {
                let zd = drop_metric
                    .iter()
                    .min_by(|a, b| a.1.abs().partial_cmp(&b.1.abs()).unwrap())
                    .map(|(d, _)| match cli.units {
                        UnitSystem::Imperial => *d / 0.9144,
                        UnitSystem::Metric => *d,
                    })
                    .unwrap_or(0.0);
                let unit = if cli.units == UnitSystem::Imperial { "yd" } else { "m" };
                eprintln!(
                    "Warning: your drop data looks zeroed near {zd:.0} {unit} (a point has ~0 drop), \
                     but --zero-range was not given."
                );
                eprintln!(
                    "         Dope-card drops are below line of sight; pass --zero-range {zd:.0} \
                     for an accurate BC. Without it, drop is treated as bore-referenced (flat fire)."
                );
            }

            run_bc_estimation_multi(
                velocity_metric,
                mass_metric,
                diameter_metric,
                &drop_metric,
                &vel_metric,
                &models,
                atmosphere,
                zero_range_m,
                sight_height_m,
                cli.units,
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
            observed,
            observation_sigma,
            mv_prior,
            bc_prior,
            predict_range,
            prediction_sigma,
            drop_unit,
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
            let temperature = UnitConverter::resolve_temperature(temperature, units)?;
            let pressure = UnitConverter::resolve_pressure(pressure, units)?;
            // Convert to imperial for calculations (internal calculations use imperial)
            let range_yd = match units {
                UnitSystem::Imperial => range,
                UnitSystem::Metric => range / 0.9144, // meters to yards
            };
            let weight_gr = match units {
                UnitSystem::Imperial => mass,
                UnitSystem::Metric => mass / GRAMS_PER_GRAIN, // grams to grains
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
            let sight_height_default = match units {
                UnitSystem::Imperial => 2.0,
                UnitSystem::Metric => 50.0,
            };
            let sight_height = sight_height.unwrap_or(sight_height_default);
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
            let bc_segments: Option<Vec<BCSegmentData>> =
                if let Some(table_dir) = &effective_bc_table_dir {
                    let mut manager = Bc5dTableManager::new(table_dir);
                    generate_bc5d_segments(
                        &mut manager,
                        caliber_in,
                        bc,
                        drag_str,
                        weight_gr,
                        None,
                        bullet_length,
                        bullet_length.is_some(),
                        None,
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

            let uncertainty_requested = observation_sigma.is_some()
                || mv_prior.is_some()
                || bc_prior.is_some()
                || !predict_range.is_empty()
                || prediction_sigma.is_some()
                || observed.iter().any(|token| token.matches(':').count() == 2);

            if uncertainty_requested {
                let default_sigma = observation_sigma.ok_or(
                    "--observation-sigma is required for uncertainty-aware truing (it supplies the primary observation's known 1-sigma error)",
                )?;
                if !default_sigma.is_finite() || default_sigma <= 0.0 {
                    return Err("--observation-sigma must be positive and finite".into());
                }
                if observed.is_empty() {
                    return Err(
                        "uncertainty-aware joint MV+BC truing requires at least one --observed RANGE:DROP[:SIGMA] in addition to the primary observation"
                            .into(),
                    );
                }
                if bc_segments.as_ref().is_some_and(|segments| !segments.is_empty()) {
                    return Err(
                        "uncertainty-aware scalar-BC truing does not support velocity-banded BC tables; omit --bc-table-dir/--bc-table-auto or fit an explicit table scale model"
                            .into(),
                    );
                }

                let mut weighted_observations = Vec::with_capacity(observed.len() + 1);
                weighted_observations.push(WeightedTruingObservationV1 {
                    range_yd,
                    drop: measured_drop,
                    sigma: default_sigma,
                });
                for token in &observed {
                    let (observation, sigma) =
                        parse_uncertain_truing_observation(token, units, default_sigma)?;
                    weighted_observations.push(WeightedTruingObservationV1 {
                        range_yd: observation.range_yd,
                        drop: observation.drop,
                        sigma,
                    });
                }

                let mv_prior = mv_prior
                    .as_deref()
                    .map(|token| {
                        let (mean, sigma) = parse_normal_prior(token, "--mv-prior")?;
                        let to_fps = match units {
                            UnitSystem::Imperial => 1.0,
                            UnitSystem::Metric => 1.0 / 0.3048,
                        };
                        Ok::<NormalPriorV1, String>(NormalPriorV1 {
                            mean: mean * to_fps,
                            sigma: sigma * to_fps,
                        })
                    })
                    .transpose()?;
                let bc_prior = bc_prior
                    .as_deref()
                    .map(|token| {
                        let (mean, sigma) = parse_normal_prior(token, "--bc-prior")?;
                        Ok::<NormalPriorV1, String>(NormalPriorV1 { mean, sigma })
                    })
                    .transpose()?;
                if prediction_sigma.is_some_and(|sigma| !sigma.is_finite() || sigma <= 0.0) {
                    return Err("--prediction-sigma must be positive and finite".into());
                }
                if prediction_sigma.is_some() && predict_range.is_empty() {
                    return Err(
                        "--prediction-sigma requires at least one --predict-range".into(),
                    );
                }
                let predictions = predict_range
                    .iter()
                    .map(|range| TruingPredictionRequestV1 {
                        range_yd: match units {
                            UnitSystem::Imperial => *range,
                            UnitSystem::Metric => *range / 0.9144,
                        },
                        future_observation_sigma: prediction_sigma,
                    })
                    .collect();
                let initial_mv_fps = mv_prior.map_or(3000.0, |prior| prior.mean);
                let model = TruingModelInputsV1 {
                    muzzle_velocity_fps: initial_mv_fps,
                    ballistic_coefficient: bc,
                    drag_model,
                    mass_gr: weight_gr,
                    diameter_in: caliber_in,
                    zero_distance_yd: zero_yd,
                    sight_height_in: sight_in,
                    temperature_f: temp_f,
                    pressure_inhg: press_inhg,
                    humidity_pct: humidity,
                    altitude_ft: alt_ft,
                };
                model
                    .validate()
                    .map_err(|error| format!("invalid uncertainty-truing model: {error}"))?;
                let request = UncertaintyTruingRequestV1 {
                    model,
                    drop_unit,
                    observations: weighted_observations,
                    priors: TruingPriorsV1 {
                        muzzle_velocity_fps: mv_prior,
                        ballistic_coefficient: bc_prior,
                    },
                    predictions,
                };
                eprintln!(
                    "Fitting {} weighted observations (uncertainty-aware joint MV+BC MAP)...",
                    request.observations.len()
                );
                let report = run_uncertainty_truing_v1(&request)?;
                display_uncertainty_truing_result(&report, units, chrono_fps, output)?;
            } else if !observed.is_empty() {
                // MBA-1316: one or more --observed impacts -> joint MV+BC
                // calibration via the real forward model (always computed
                // locally; the online API only supports single-observation
                // velocity truing).
                run_multi_observation_truing(
                    measured_drop,
                    range_yd,
                    &observed,
                    drop_unit,
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
                    chrono_fps,
                    units,
                    output,
                )?;
            } else if use_local {
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
                        let velocity_adjustment =
                            chrono_fps.map(|c| result.effective_velocity_fps - c);
                        let adjustment_percent =
                            chrono_fps.map(|c| (result.effective_velocity_fps - c) / c * 100.0);

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
                    if !check_tos_accepted() && !prompt_tos_acceptance()? {
                        eprintln!(
                            "Cannot use online features without accepting Terms of Service."
                        );
                        std::process::exit(1);
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
                                        let velocity_adjustment =
                                            chrono_fps.map(|c| result.effective_velocity_fps - c);
                                        let adjustment_percent = chrono_fps.map(|c| {
                                            (result.effective_velocity_fps - c) / c * 100.0
                                        });

                                        let effective_vel = match units {
                                            UnitSystem::Imperial => result.effective_velocity_fps,
                                            UnitSystem::Metric => {
                                                result.effective_velocity_fps * 0.3048
                                            }
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
                                        eprintln!(
                                            "Error: Fallback calculation also failed: {}",
                                            fallback_err
                                        );
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

        Commands::PlanTruing {
            profile,
            velocity,
            bc,
            drag_model,
            mass,
            diameter,
            candidate_ranges,
            range_grid,
            observation_count,
            minimum_separation,
            measurement_resolution,
            drop_unit,
            zero_distance,
            sight_height,
            temperature,
            pressure,
            humidity,
            altitude,
            output,
        } => {
            let profile_data = profile
                .as_deref()
                .map(|name| load_profile_for_units(name, cli.units))
                .transpose()?;
            if let Some(profile) = &profile_data {
                let has_segments = profile
                    .bc_segments
                    .as_ref()
                    .is_some_and(|segments| !segments.is_empty())
                    || profile.use_bc_segments.unwrap_or(false);
                let has_custom_drag = profile
                    .drag_curve
                    .as_ref()
                    .is_some_and(|curve| !curve.is_empty())
                    || profile.drag_model.eq_ignore_ascii_case("custom");
                if has_segments || has_custom_drag {
                    return Err(format!(
                        "profile '{}' uses {} and cannot be reduced to the scalar G1/G7 BC parameter designed by plan-truing; use a scalar-BC profile or define a separate drag-deck scale model",
                        profile.name,
                        if has_custom_drag {
                            "a custom drag curve"
                        } else {
                            "velocity-banded BC values"
                        }
                    )
                    .into());
                }
            }

            let velocity = velocity
                .or_else(|| profile_data.as_ref().map(|profile| profile.velocity))
                .ok_or("--velocity is required unless --profile supplies it")?;
            let bc = bc
                .or_else(|| profile_data.as_ref().map(|profile| profile.bc))
                .ok_or("--bc is required unless --profile supplies it")?;
            let mass = mass
                .or_else(|| profile_data.as_ref().map(|profile| profile.mass))
                .ok_or("--mass is required unless --profile supplies it")?;
            let diameter = diameter
                .or_else(|| profile_data.as_ref().map(|profile| profile.diameter))
                .ok_or("--diameter is required unless --profile supplies it")?;
            let drag_model = drag_model
                .or_else(|| {
                    profile_data
                        .as_ref()
                        .map(|profile| parse_drag_model_arg(&profile.drag_model))
                })
                .unwrap_or(DragModelArg::G1);
            let zero_distance = zero_distance
                .or_else(|| {
                    profile_data
                        .as_ref()
                        .and_then(|profile| profile.zero_distance.or(profile.auto_zero))
                })
                .unwrap_or(100.0);
            let sight_height = sight_height
                .or_else(|| {
                    profile_data
                        .as_ref()
                        .and_then(|profile| profile.sight_height)
                })
                .unwrap_or(match cli.units {
                    UnitSystem::Imperial => 2.0,
                    UnitSystem::Metric => 50.0,
                });
            let temperature = UnitConverter::resolve_temperature(
                temperature.or_else(|| {
                    profile_data
                        .as_ref()
                        .map(|profile| profile.temperature)
                }),
                cli.units,
            )?;
            let pressure = UnitConverter::resolve_pressure(
                pressure.or_else(|| profile_data.as_ref().map(|profile| profile.pressure)),
                cli.units,
            )?;
            let humidity = humidity
                .or_else(|| profile_data.as_ref().map(|profile| profile.humidity))
                .unwrap_or(50.0);
            let altitude = altitude
                .or_else(|| profile_data.as_ref().map(|profile| profile.altitude))
                .unwrap_or(0.0);

            let candidate_ranges = match range_grid.as_deref() {
                Some(grid) => parse_truing_range_grid(grid)?,
                None => candidate_ranges,
            };
            let to_yards = |value: f64| match cli.units {
                UnitSystem::Imperial => value,
                UnitSystem::Metric => value / 0.9144,
            };
            let model = TruingModelInputsV1 {
                muzzle_velocity_fps: match cli.units {
                    UnitSystem::Imperial => velocity,
                    UnitSystem::Metric => velocity / 0.3048,
                },
                ballistic_coefficient: bc,
                drag_model,
                mass_gr: match cli.units {
                    UnitSystem::Imperial => mass,
                    UnitSystem::Metric => mass / GRAMS_PER_GRAIN,
                },
                diameter_in: match cli.units {
                    UnitSystem::Imperial => diameter,
                    UnitSystem::Metric => diameter / 25.4,
                },
                zero_distance_yd: to_yards(zero_distance),
                sight_height_in: match cli.units {
                    UnitSystem::Imperial => sight_height,
                    UnitSystem::Metric => sight_height / 25.4,
                },
                temperature_f: match cli.units {
                    UnitSystem::Imperial => temperature,
                    UnitSystem::Metric => temperature * 9.0 / 5.0 + 32.0,
                },
                pressure_inhg: match cli.units {
                    UnitSystem::Imperial => pressure,
                    UnitSystem::Metric => pressure / 33.8639,
                },
                humidity_pct: humidity,
                altitude_ft: match cli.units {
                    UnitSystem::Imperial => altitude,
                    UnitSystem::Metric => altitude / 0.3048,
                },
            };
            let request = TruingExperimentPlanRequestV1 {
                model,
                candidate_ranges_yd: candidate_ranges.into_iter().map(to_yards).collect(),
                observation_count,
                minimum_separation_yd: to_yards(minimum_separation),
                measurement_sigma_1sd: measurement_resolution,
                drop_unit,
            };
            eprintln!(
                "Evaluating {} candidate ranges for a {}-station truing plan...",
                request.candidate_ranges_yd.len(),
                request.observation_count
            );
            let plan = plan_truing_experiment_v1(&request)?;
            display_truing_experiment_plan(&plan, cli.units, output)?;
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
            let temperature = UnitConverter::resolve_temperature(temperature, cli.units)?;
            let pressure = UnitConverter::resolve_pressure(pressure, cli.units)?;
            // Load profile if specified
            let profile_data = profile.as_ref().map(|name| {
                load_profile_for_units(name, cli.units).unwrap_or_else(|e| {
                    eprintln!("Error: {}", e);
                    std::process::exit(1);
                })
            });

            let final_velocity = resolve_param(velocity, &profile_data, |p| p.velocity)
                .unwrap_or_else(|| {
                    eprintln!("Error: --velocity is required (or use --profile)");
                    std::process::exit(1);
                });
            let final_bc = resolve_param(bc, &profile_data, |p| p.bc).unwrap_or_else(|| {
                eprintln!("Error: --bc is required (or use --profile)");
                std::process::exit(1);
            });
            let final_mass = resolve_param(mass, &profile_data, |p| p.mass).unwrap_or_else(|| {
                eprintln!("Error: --mass is required (or use --profile)");
                std::process::exit(1);
            });
            let final_diameter = resolve_param(diameter, &profile_data, |p| p.diameter)
                .unwrap_or_else(|| {
                    eprintln!("Error: --diameter is required (or use --profile)");
                    std::process::exit(1);
                });
            let final_sight_height = sight_height
                .or_else(|| profile_data.as_ref().and_then(|p| p.sight_height))
                .unwrap_or(match cli.units {
                    UnitSystem::Imperial => 2.0,
                    UnitSystem::Metric => 50.0,
                });
            let final_drag_model = match (profile_data.as_ref(), velocity) {
                (Some(profile), None) => parse_drag_model_arg(&profile.drag_model),
                _ => drag_model,
            };

            handle_mpbr(
                final_velocity,
                final_bc,
                final_mass,
                final_diameter,
                final_drag_model,
                vital_zone,
                final_sight_height,
                temperature,
                pressure,
                humidity,
                altitude,
                cli.units,
                output,
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
            elevation_click_value,
            windage_click_value,
            sight_height,
            temperature,
            pressure,
            humidity,
            altitude,
            wind_speed,
            wind_direction,
            output,
        } => {
            let temperature = UnitConverter::resolve_temperature(temperature, cli.units)?;
            let pressure = UnitConverter::resolve_pressure(pressure, cli.units)?;
            let profile_data = profile.as_ref().map(|name| {
                load_profile_for_units(name, cli.units).unwrap_or_else(|e| {
                    eprintln!("Error: {}", e);
                    std::process::exit(1);
                })
            });

            // MBA-1355: resolve the elevation click graduation FIRST — before any of the
            // work below — so a missing graduation fails fast with a clear flag name.
            // come-ups has no windage column, so the resolved windage value (which would
            // otherwise just default to the elevation one) is intentionally discarded.
            let elevation_click: Option<ClickValue> = if matches!(adjustment_unit, AdjustmentUnit::Clicks) {
                match resolve_click_values(
                    elevation_click_value.as_deref(),
                    windage_click_value.as_deref(),
                    profile_data.as_ref(),
                ) {
                    Ok(Some((el, _wi))) => Some(el),
                    Ok(None) => None,
                    Err(e) => {
                        eprintln!("error: {e}");
                        std::process::exit(1);
                    }
                }
            } else {
                None
            };

            let final_velocity = resolve_param(velocity, &profile_data, |p| p.velocity)
                .unwrap_or_else(|| {
                    eprintln!("Error: --velocity is required (or use --profile)");
                    std::process::exit(1);
                });
            let final_bc = resolve_param(bc, &profile_data, |p| p.bc).unwrap_or_else(|| {
                eprintln!("Error: --bc is required (or use --profile)");
                std::process::exit(1);
            });
            let final_mass = resolve_param(mass, &profile_data, |p| p.mass).unwrap_or_else(|| {
                eprintln!("Error: --mass is required (or use --profile)");
                std::process::exit(1);
            });
            let final_diameter = resolve_param(diameter, &profile_data, |p| p.diameter)
                .unwrap_or_else(|| {
                    eprintln!("Error: --diameter is required (or use --profile)");
                    std::process::exit(1);
                });
            let final_sight_height = sight_height
                .or_else(|| profile_data.as_ref().and_then(|p| p.sight_height))
                .unwrap_or(match cli.units {
                    UnitSystem::Imperial => 2.0,
                    UnitSystem::Metric => 50.0,
                });
            let final_drag_model = match (profile_data.as_ref(), velocity) {
                (Some(profile), None) => parse_drag_model_arg(&profile.drag_model),
                _ => drag_model,
            };
            // MBA-1323 Phase 2: a `.a7p`-imported profile's velocity-BC segments / Mach-Cd
            // drag curve, resolved to engine shapes. come-ups has no --bc-segment/--drag-table
            // CLI flags of its own, so a saved profile is the only source for these.
            let bc_segments_data = profile_data
                .as_ref()
                .and_then(|p| p.bc_segments.as_ref())
                .map(|rows| bc_segments_from_profile(rows));
            let custom_drag_table = profile_data.as_ref().and_then(|p| p.drag_curve.as_ref()).map(|pts| {
                drag_table_from_profile(pts).unwrap_or_else(|e| {
                    eprintln!("Error: saved profile's drag curve is invalid: {e}");
                    std::process::exit(1);
                })
            });

            handle_come_ups(
                final_velocity,
                final_bc,
                final_mass,
                final_diameter,
                final_drag_model,
                zero_distance,
                start,
                end,
                step,
                adjustment_unit,
                elevation_click,
                final_sight_height,
                temperature,
                pressure,
                humidity,
                altitude,
                wind_speed,
                wind_direction,
                cli.units,
                output,
                bc_segments_data,
                custom_drag_table,
            )?;
        }

        Commands::Lead {
            profile,
            velocity,
            bc,
            mass,
            diameter,
            drag_model,
            sight_height,
            temperature,
            pressure,
            humidity,
            altitude,
            wind_speed,
            wind_direction,
            use_powder_sensitivity,
            powder_temp_sensitivity,
            powder_temp,
            powder_temp_curve,
            target_speed,
            target_angle,
            target_length,
            start,
            end,
            step,
            adjustment_unit,
            output,
        } => {
            let temperature = UnitConverter::resolve_temperature(temperature, cli.units)?;
            let pressure = UnitConverter::resolve_pressure(pressure, cli.units)?;
            let profile_data = profile.as_ref().map(|name| {
                load_profile_for_units(name, cli.units).unwrap_or_else(|e| {
                    eprintln!("Error: {}", e);
                    std::process::exit(1);
                })
            });

            let final_velocity = resolve_param(velocity, &profile_data, |p| p.velocity)
                .unwrap_or_else(|| {
                    eprintln!("Error: --velocity is required (or use --profile)");
                    std::process::exit(1);
                });
            let final_bc = resolve_param(bc, &profile_data, |p| p.bc).unwrap_or_else(|| {
                eprintln!("Error: --bc is required (or use --profile)");
                std::process::exit(1);
            });
            let final_mass = resolve_param(mass, &profile_data, |p| p.mass).unwrap_or_else(|| {
                eprintln!("Error: --mass is required (or use --profile)");
                std::process::exit(1);
            });
            let final_diameter = resolve_param(diameter, &profile_data, |p| p.diameter)
                .unwrap_or_else(|| {
                    eprintln!("Error: --diameter is required (or use --profile)");
                    std::process::exit(1);
                });
            let final_sight_height = sight_height
                .or_else(|| profile_data.as_ref().and_then(|p| p.sight_height))
                .unwrap_or(match cli.units {
                    UnitSystem::Imperial => 2.0,
                    UnitSystem::Metric => 50.0,
                });
            let final_drag_model = match (&profile_data, velocity) {
                (Some(p), None) => parse_drag_model_arg(&p.drag_model),
                _ => drag_model,
            };

            // Powder temperature plumbing (MBA-1325), identical resolution to the
            // `trajectory` command: --powder-temp serves as either the curve's lookup
            // temperature (defaults to --temperature) or the linear model's reference
            // temperature (defaults to the 70F/21C equivalent); powder_curve_temp_c
            // captures whether it was explicitly given (Some) so the curve can fall
            // back to ambient temperature (None) rather than a synthesized default.
            let powder_curve_temp_c: Option<f64> =
                powder_temp.map(|t| UnitConverter::temperature_to_metric(t, cli.units));
            let powder_temp = powder_temp.unwrap_or(match cli.units {
                UnitSystem::Imperial => DEFAULT_POWDER_REFERENCE_TEMP_F,
                UnitSystem::Metric => DEFAULT_POWDER_REFERENCE_TEMP_C,
            });
            let powder_temp_sensitivity = powder_temp_sensitivity.unwrap_or(match cli.units {
                UnitSystem::Imperial => 1.0,
                UnitSystem::Metric => 0.3048 / (5.0 / 9.0),
            });
            let powder_temp_curve_si: Option<Vec<(f64, f64)>> = match powder_temp_curve.as_deref()
            {
                Some(s) => Some(parse_powder_temp_curve(s, cli.units)?),
                None => None,
            };

            // MBA-1323 Phase 2: a `.a7p`-imported profile's velocity-BC segments / Mach-Cd
            // drag curve, resolved to engine shapes. lead has no --bc-segment/--drag-table
            // CLI flags of its own, so a saved profile is the only source for these.
            let bc_segments_data = profile_data
                .as_ref()
                .and_then(|p| p.bc_segments.as_ref())
                .map(|rows| bc_segments_from_profile(rows));
            let custom_drag_table = profile_data.as_ref().and_then(|p| p.drag_curve.as_ref()).map(|pts| {
                drag_table_from_profile(pts).unwrap_or_else(|e| {
                    eprintln!("Error: saved profile's drag curve is invalid: {e}");
                    std::process::exit(1);
                })
            });

            handle_lead(
                final_velocity,
                final_bc,
                final_mass,
                final_diameter,
                final_drag_model,
                final_sight_height,
                temperature,
                pressure,
                humidity,
                altitude,
                wind_speed,
                wind_direction,
                use_powder_sensitivity,
                powder_temp_sensitivity,
                powder_temp,
                powder_temp_curve_si,
                powder_curve_temp_c,
                target_speed,
                target_angle,
                target_length,
                start,
                end,
                step,
                adjustment_unit,
                cli.units,
                output,
                bc_segments_data,
                custom_drag_table,
            )?;
        }

        Commands::Powder {
            velocity,
            temperature,
            powder_temp_sensitivity,
            powder_temp,
            powder_temp_curve,
            sweep,
            mass,
            output,
        } => {
            let powder_temp_curve_si: Option<Vec<(f64, f64)>> = match powder_temp_curve.as_deref()
            {
                Some(s) => Some(parse_powder_temp_curve(s, cli.units)?),
                None => None,
            };
            handle_powder(
                velocity,
                temperature,
                powder_temp_sensitivity,
                powder_temp,
                powder_temp_curve_si,
                sweep,
                mass,
                cli.units,
                output,
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
            wind_angle,
            wind_angles,
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
            let temperature = UnitConverter::resolve_temperature(temperature, cli.units)?;
            let pressure = UnitConverter::resolve_pressure(pressure, cli.units)?;
            let profile_data = profile.as_ref().map(|name| {
                load_profile_for_units(name, cli.units).unwrap_or_else(|e| {
                    eprintln!("Error: {}", e);
                    std::process::exit(1);
                })
            });

            let final_velocity = resolve_param(velocity, &profile_data, |p| p.velocity)
                .unwrap_or_else(|| {
                    eprintln!("Error: --velocity is required (or use --profile)");
                    std::process::exit(1);
                });
            let final_bc = resolve_param(bc, &profile_data, |p| p.bc).unwrap_or_else(|| {
                eprintln!("Error: --bc is required (or use --profile)");
                std::process::exit(1);
            });
            let final_mass = resolve_param(mass, &profile_data, |p| p.mass).unwrap_or_else(|| {
                eprintln!("Error: --mass is required (or use --profile)");
                std::process::exit(1);
            });
            let final_diameter = resolve_param(diameter, &profile_data, |p| p.diameter)
                .unwrap_or_else(|| {
                    eprintln!("Error: --diameter is required (or use --profile)");
                    std::process::exit(1);
                });
            let final_sight_height = sight_height
                .or_else(|| profile_data.as_ref().and_then(|p| p.sight_height))
                .unwrap_or(match cli.units {
                    UnitSystem::Imperial => 2.0,
                    UnitSystem::Metric => 50.0,
                });
            let final_drag_model = match (profile_data.as_ref(), velocity) {
                (Some(profile), None) => parse_drag_model_arg(&profile.drag_model),
                _ => drag_model,
            };

            // Parse wind speeds
            let ws_vec: Vec<f64> = wind_speeds
                .split(',')
                .filter_map(|s| s.trim().parse::<f64>().ok())
                .collect();
            if ws_vec.is_empty() {
                eprintln!("Error: --wind-speeds must contain at least one valid number (e.g., '5,10,15,20')");
                std::process::exit(1);
            }

            // Resolve wind angle(s): neither flag → legacy full-value 90° card
            // (byte-identical to pre-MBA-727 output); --wind-angle → single
            // angle-aware card; --wind-angles → one card per angle.
            let (wind_angle_vec, legacy_labels): (Vec<f64>, bool) =
                if let Some(angle) = wind_angle {
                    (vec![angle], false)
                } else if let Some(csv) = wind_angles.as_ref() {
                    let mut angles = Vec::new();
                    for tok in csv.split(',') {
                        let tok = tok.trim();
                        let angle: f64 = tok.parse().map_err(|_| {
                            format!("--wind-angles contains an invalid number: '{}'", tok)
                        })?;
                        if angle.is_nan() || !(0.0..=360.0).contains(&angle) {
                            return Err(format!(
                                "--wind-angles value '{}' must be in [0, 360]",
                                tok
                            )
                            .into());
                        }
                        angles.push(angle);
                    }
                    // (no empty-check needed: split(',') yields >=1 token and every
                    // token either pushes or errors via `?`)
                    (angles, false)
                } else {
                    (vec![90.0], true)
                };

            handle_wind_card(
                final_velocity,
                final_bc,
                final_mass,
                final_diameter,
                final_drag_model,
                zero_distance,
                &ws_vec,
                &wind_angle_vec,
                legacy_labels,
                start,
                end,
                step,
                adjustment_unit,
                final_sight_height,
                temperature,
                pressure,
                humidity,
                altitude,
                cli.units,
                output,
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
            let temperature = UnitConverter::resolve_temperature(temperature, cli.units)?;
            let pressure = UnitConverter::resolve_pressure(pressure, cli.units)?;
            let profile_data = profile.as_ref().map(|name| {
                load_profile_for_units(name, cli.units).unwrap_or_else(|e| {
                    eprintln!("Error: {}", e);
                    std::process::exit(1);
                })
            });

            let final_mass = resolve_param(mass, &profile_data, |p| p.mass).unwrap_or_else(|| {
                eprintln!("Error: --mass is required (or use --profile)");
                std::process::exit(1);
            });
            let final_diameter = resolve_param(diameter, &profile_data, |p| p.diameter)
                .unwrap_or_else(|| {
                    eprintln!("Error: --diameter is required (or use --profile)");
                    std::process::exit(1);
                });
            let final_velocity = velocity.unwrap_or(match cli.units {
                UnitSystem::Imperial => 2700.0,
                UnitSystem::Metric => 823.0,
            });

            handle_stability(
                final_mass,
                final_diameter,
                length,
                twist_rate,
                final_velocity,
                temperature,
                pressure,
                altitude,
                cli.units,
                output,
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
            let temperature = UnitConverter::resolve_temperature(temperature, cli.units)?;
            let pressure = UnitConverter::resolve_pressure(pressure, cli.units)?;
            let profile_data = profile.as_ref().map(|name| {
                load_profile_for_units(name, cli.units).unwrap_or_else(|e| {
                    eprintln!("Error: {}", e);
                    std::process::exit(1);
                })
            });

            let final_velocity = resolve_param(velocity, &profile_data, |p| p.velocity)
                .unwrap_or_else(|| {
                    eprintln!("Error: --velocity is required (or use --profile)");
                    std::process::exit(1);
                });
            let final_bc = resolve_param(bc, &profile_data, |p| p.bc).unwrap_or_else(|| {
                eprintln!("Error: --bc is required (or use --profile)");
                std::process::exit(1);
            });
            let final_mass = resolve_param(mass, &profile_data, |p| p.mass).unwrap_or_else(|| {
                eprintln!("Error: --mass is required (or use --profile)");
                std::process::exit(1);
            });
            let final_diameter = resolve_param(diameter, &profile_data, |p| p.diameter)
                .unwrap_or_else(|| {
                    eprintln!("Error: --diameter is required (or use --profile)");
                    std::process::exit(1);
                });
            let final_sight_height = sight_height
                .or_else(|| profile_data.as_ref().and_then(|p| p.sight_height))
                .unwrap_or(match cli.units {
                    UnitSystem::Imperial => 2.0,
                    UnitSystem::Metric => 50.0,
                });
            let final_drag_model = match (profile_data.as_ref(), velocity) {
                (Some(profile), None) => parse_drag_model_arg(&profile.drag_model),
                _ => drag_model,
            };

            handle_range_table(
                final_velocity,
                final_bc,
                final_mass,
                final_diameter,
                final_drag_model,
                zero_distance,
                start,
                end,
                step,
                wind_speed,
                wind_direction,
                adjustment_unit,
                final_sight_height,
                temperature,
                pressure,
                humidity,
                altitude,
                cli.units,
                output,
            )?;
        }

        Commands::Compare {
            loads,
            profiles,
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
            let temperature = UnitConverter::resolve_temperature(temperature, cli.units)?;
            let pressure = UnitConverter::resolve_pressure(pressure, cli.units)?;

            // Assemble the load list: explicit --load specs first, then --profile loads,
            // preserving command-line order within each group. Load #1 is the delta baseline.
            let mut compare_loads: Vec<CompareLoad> = Vec::new();
            for spec in &loads {
                compare_loads.push(parse_compare_load_spec(spec, cli.units)?);
            }
            for name in &profiles {
                let p = load_profile_for_units(name, cli.units)?;
                // MBA-1323 Phase 2 shape: a profile's velocity-BC segments / Mach-Cd
                // drag curve, resolved to engine shapes and used by zeroing AND the
                // sampled runs (the scalar-BC-only caveat is lifted for compare).
                let bc_segments_data = p
                    .bc_segments
                    .as_ref()
                    .map(|rows| bc_segments_from_profile(rows));
                let custom_drag_table = p
                    .drag_curve
                    .as_ref()
                    .map(|pts| {
                        drag_table_from_profile(pts).map_err(|e| {
                            format!("profile '{name}': saved drag curve is invalid: {e}")
                        })
                    })
                    .transpose()?;
                compare_loads.push(CompareLoad {
                    name: name.clone(),
                    drag_model: parse_drag_model_arg(&p.drag_model),
                    bc: p.bc,
                    mass: p.mass,
                    velocity: p.velocity,
                    diameter: p.diameter,
                    bc_segments_data,
                    custom_drag_table,
                });
            }
            if compare_loads.len() < 2 {
                return Err("compare needs at least 2 loads (via --load and/or --profile)".into());
            }
            if compare_loads.len() > 8 {
                return Err("compare supports at most 8 loads".into());
            }
            if start >= end {
                return Err("--start must be less than --end".into());
            }

            handle_compare(
                compare_loads,
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
                cli.units,
                output,
            )?;
        }

        Commands::Profile { action } => {
            match action {
                ProfileAction::Save {
                    name,
                    velocity,
                    bc,
                    mass,
                    diameter,
                    drag_model,
                    twist_rate,
                    sight_height,
                    zero_distance,
                    temperature,
                    pressure,
                    humidity,
                    altitude,
                    bullet_name,
                    wind_speed,
                    wind_direction,
                    shooting_angle,
                    auto_zero,
                    twist_right,
                    use_bc_segments,
                    bullet_length,
                    elevation_click,
                    windage_click,
                } => {
                    let temperature = UnitConverter::resolve_temperature(temperature, cli.units)?;
                    let pressure = UnitConverter::resolve_pressure(pressure, cli.units)?;
                    let drag_str = match drag_model {
                        DragModelArg::G1 => "G1",
                        DragModelArg::G7 => "G7",
                    };
                    let units_str = match cli.units {
                        UnitSystem::Imperial => "imperial",
                        UnitSystem::Metric => "metric",
                    };

                    // MBA-1355: validate click graduations at save time so a saved profile
                    // can never store a value `resolve_click_values` would later reject.
                    if let Some(ref v) = elevation_click {
                        parse_click_value(v)
                            .map_err(|e| format!("--elevation-click '{v}' is invalid: {e}"))?;
                    }
                    if let Some(ref v) = windage_click {
                        parse_click_value(v)
                            .map_err(|e| format!("--windage-click '{v}' is invalid: {e}"))?;
                    }

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
                        wind_speed,
                        wind_direction,
                        shooting_angle,
                        auto_zero,
                        twist_right,
                        use_bc_segments,
                        bullet_length,
                        elevation_click,
                        windage_click,
                        // `profile save` only accepts a scalar BC/drag-model today; a full
                        // velocity-banded schedule or Mach/Cd curve can only be produced by
                        // `.a7p` import (see map_a7p_to_profile).
                        bc_segments: None,
                        drag_curve: None,
                    };

                    let path = save_profile(&profile)?;
                    eprintln!("Profile '{}' saved to {:?}", name, path);
                }

                #[cfg(feature = "profile-import")]
                ProfileAction::Import {
                    file,
                    name,
                    dry_run,
                    strict,
                } => {
                    use ballistics_engine::profile_import::{parse_a7p, EnvelopeStatus};
                    let bytes = fs::read(&file)
                        .map_err(|e| format!("cannot read {}: {e}", file.display()))?;
                    let doc = parse_a7p(&bytes)
                        .map_err(|e| format!("not a usable .a7p file: {e}"))?;
                    if strict {
                        if let EnvelopeStatus::Mismatch { expected, actual } = &doc.envelope {
                            return Err(format!(
                                "checksum mismatch (file says {expected}, payload hashes to {actual}) — refusing under --strict"
                            )
                            .into());
                        }
                    }
                    // Sanitize the --name override the same way a name derived from the
                    // file itself is sanitized (map_a7p_to_profile does that internally
                    // when there is no override) — an unsanitized override would otherwise
                    // be written verbatim into the profile store's file path.
                    let sanitized_name = name.as_deref().map(|raw| {
                        let sanitized = sanitize_profile_name(raw);
                        if sanitized != raw {
                            println!(
                                "note: profile name sanitized to '{sanitized}' for the profile store"
                            );
                        }
                        sanitized
                    });
                    let outcome = map_a7p_to_profile(&doc, sanitized_name.as_deref())?;
                    print!("{}", render_import_report(&outcome.report));
                    if dry_run {
                        println!("\nDry run — nothing saved.");
                    } else {
                        let path = save_profile(&outcome.profile)?;
                        println!(
                            "\nSaved profile '{}' -> {}",
                            outcome.profile.name,
                            path.display()
                        );
                    }
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
                            println!(
                                "│ {:<18} │{:>7.0} │{:>6.3} │{:>7.1} │{:>9.3} │ {:<8} │",
                                p.name, p.velocity, p.bc, p.mass, p.diameter, p.drag_model
                            );
                        }
                        println!("└────────────────────┴────────┴───────┴────────┴──────────┴──────────┘");
                    }
                }

                ProfileAction::Show { name } => {
                    let profile = load_profile(&name)?;
                    println!();
                    println!("Profile: {}", profile.name);
                    println!("╔════════════════════════════════════════╗");
                    println!(
                        "║  Velocity:      {:>10.1} {:<10}  ║",
                        profile.velocity,
                        if profile.units == "metric" {
                            "m/s"
                        } else {
                            "fps"
                        }
                    );
                    if profile.drag_curve.is_some() {
                        // No scalar BC applies to a full drag curve — see
                        // map_a7p_to_profile's CUSTOM handling. Printing the inert 0.0
                        // sentinel here would read as a real coefficient, so it is
                        // replaced by the "Drag curve:" summary line below instead.
                        println!("║  BC:            {:<23}║", "not applicable (CUSTOM)");
                    } else {
                        println!("║  BC:            {:>10.4}             ║", profile.bc);
                    }
                    println!(
                        "║  Mass:          {:>10.1} {:<10}  ║",
                        profile.mass,
                        if profile.units == "metric" { "g" } else { "gr" }
                    );
                    println!(
                        "║  Diameter:      {:>10.3} {:<10}  ║",
                        profile.diameter,
                        if profile.units == "metric" {
                            "mm"
                        } else {
                            "in"
                        }
                    );
                    println!("║  Drag model:    {:>10}             ║", profile.drag_model);
                    if let Some(tw) = profile.twist_rate {
                        println!("║  Twist rate:    {:>10.1}             ║", tw);
                    }
                    if let Some(sh) = profile.sight_height {
                        println!(
                            "║  Sight height:  {:>10.2} {:<10}  ║",
                            sh,
                            if profile.units == "metric" {
                                "mm"
                            } else {
                                "in"
                            }
                        );
                    }
                    if let Some(zd) = profile.zero_distance {
                        println!(
                            "║  Zero distance: {:>10.0} {:<10}  ║",
                            zd,
                            if profile.units == "metric" { "m" } else { "yd" }
                        );
                    }
                    if let Some(ref bn) = profile.bullet_name {
                        println!("║  Bullet:        {:<24}  ║", bn);
                    }
                    if let Some(ws) = profile.wind_speed {
                        println!(
                            "║  Wind speed:    {:>10.1} {:<10}  ║",
                            ws,
                            if profile.units == "metric" {
                                "m/s"
                            } else {
                                "mph"
                            }
                        );
                    }
                    if let Some(wd) = profile.wind_direction {
                        println!("║  Wind dir:      {:>10.1}°            ║", wd);
                    }
                    if let Some(sa) = profile.shooting_angle {
                        println!("║  Shoot angle:   {:>10.1}°            ║", sa);
                    }
                    if let Some(az) = profile.auto_zero {
                        println!(
                            "║  Auto-zero:     {:>10.0} {:<10}  ║",
                            az,
                            if profile.units == "metric" { "m" } else { "yd" }
                        );
                    }
                    if let Some(tr) = profile.twist_right {
                        println!(
                            "║  Twist:         {:>10}             ║",
                            if tr { "Right" } else { "Left" }
                        );
                    }
                    if let Some(ubc) = profile.use_bc_segments {
                        println!(
                            "║  BC segments:   {:>10}             ║",
                            if ubc { "Enabled" } else { "Disabled" }
                        );
                    }
                    if let Some(bl) = profile.bullet_length {
                        println!(
                            "║  Bullet length: {:>10.3} {:<10}  ║",
                            bl,
                            if profile.units == "metric" {
                                "mm"
                            } else {
                                "in"
                            }
                        );
                    }
                    if let Some(ref ec) = profile.elevation_click {
                        println!("║  Elev click:    {:<23}║", ec);
                    }
                    if let Some(ref wc) = profile.windage_click {
                        println!("║  Wind click:    {:<23}║", wc);
                    }
                    // MBA-1323 Phase 2: velocity-BC bands / Mach-Cd drag curve summaries
                    // (count + range). velocity_mps is always SI (see ProfileBcSegment).
                    if let Some(ref segs) = profile.bc_segments {
                        let vmin = segs
                            .iter()
                            .map(|s| s.velocity_mps)
                            .fold(f64::INFINITY, f64::min);
                        let vmax = segs
                            .iter()
                            .map(|s| s.velocity_mps)
                            .fold(f64::NEG_INFINITY, f64::max);
                        println!(
                            "║  BC bands:      {:<23}║",
                            format!("{} rows, {:.0}-{:.0} m/s", segs.len(), vmin, vmax)
                        );
                    }
                    if let Some(ref curve) = profile.drag_curve {
                        let mmin = curve.iter().map(|p| p.mach).fold(f64::INFINITY, f64::min);
                        let mmax = curve
                            .iter()
                            .map(|p| p.mach)
                            .fold(f64::NEG_INFINITY, f64::max);
                        println!(
                            "║  Drag curve:    {:<23}║",
                            format!("{} pts, Mach {:.2}-{:.2}", curve.len(), mmin, mmax)
                        );
                    }
                    println!(
                        "║  Temperature:   {:>10.1} {:<10}  ║",
                        profile.temperature,
                        if profile.units == "metric" {
                            "°C"
                        } else {
                            "°F"
                        }
                    );
                    println!(
                        "║  Pressure:      {:>10.2} {:<10}  ║",
                        profile.pressure,
                        if profile.units == "metric" {
                            "hPa"
                        } else {
                            "inHg"
                        }
                    );
                    println!("║  Humidity:      {:>10.1}%            ║", profile.humidity);
                    println!(
                        "║  Altitude:      {:>10.0} {:<10}  ║",
                        profile.altitude,
                        if profile.units == "metric" { "m" } else { "ft" }
                    );
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
    adjustment_unit: AdjustmentUnit,
}

// ============================================================================
// BC5D Segment Generation Helper (MBA-744)
// ============================================================================

/// Generate velocity-dependent BC segments from a Bc5dTableManager.
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments preserve the existing BC5D compatibility helper"
)]
fn generate_bc5d_segments(
    manager: &mut Bc5dTableManager,
    caliber: f64,
    base_bc: f64,
    drag_type: &str,
    weight_grains: f64,
    muzzle_velocity_fps: Option<f64>,
    bullet_length_in: Option<f64>,
    length_is_user: bool,
    print_ladder_units: Option<UnitSystem>,
) -> Option<Vec<BCSegmentData>> {
    let mut velocity_breakpoints: Vec<f64> = vec![
        4000.0, 3500.0, 3000.0, 2700.0, 2500.0, 2300.0, 2100.0, 2000.0, 1900.0, 1800.0, 1700.0,
        1600.0, 1500.0, 1400.0, 1350.0, 1300.0, 1250.0, 1200.0, 1150.0, 1100.0, 1050.0, 1000.0,
        950.0, 900.0, 850.0, 800.0, 700.0, 600.0, 500.0,
    ];

    if let Some(mv) = muzzle_velocity_fps {
        velocity_breakpoints.push(mv);
    }

    let mut valid_velocities: Vec<f64> = velocity_breakpoints
        .into_iter()
        .filter(|&v| v >= 500.0 && muzzle_velocity_fps.is_none_or(|mv| v <= mv))
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

        if let Ok(correction) = manager.lookup(
            caliber,
            weight_grains,
            base_bc,
            reference_muzzle_velocity,
            vel_mid,
            drag_type,
        ) {
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
        let length_display = bullet_length_in.unwrap_or_else(|| {
            // MBA-1135: mass-based length estimate (was a mass-blind caliber*3.5 heuristic).
            let est_m = ballistics_engine::stability::estimate_bullet_length_m(
                caliber * 0.0254,
                weight_grains * GRAINS_TO_KG,
            );
            if est_m > 0.0 {
                est_m / 0.0254
            } else {
                caliber * 3.5
            }
        });
        eprintln!(
            "BC5D Table: Generated {} velocity-dependent BC segments",
            segments.len()
        );
        // Length is informational only: the v2 5D table axes are
        // [drag_type][weight][bc][muzzle_vel][current_vel] — length is NOT a
        // lookup dimension. Label the value's source honestly (a previous
        // version printed "(est)" even for user-supplied values).
        eprintln!(
            "BC5D Table: {:.3}\" caliber, {:.1}gr, {:.3}\" length ({}; not a lookup axis)",
            caliber,
            weight_grains,
            length_display,
            if length_is_user { "user" } else { "est" }
        );
        let min_bc = segments
            .iter()
            .map(|s| s.bc_value)
            .fold(f64::INFINITY, f64::min);
        let max_bc = segments
            .iter()
            .map(|s| s.bc_value)
            .fold(f64::NEG_INFINITY, f64::max);
        eprintln!(
            "BC5D Table: BC range {:.5} - {:.5} across velocity envelope",
            min_bc, max_bc
        );
        // --print-bc-segments: dump the ladder as ready-to-paste --bc-segment
        // arguments (e.g. to run BC5D-equivalent corrections offline in the
        // WASM CLI, which accepts --bc-segment but cannot hold the tables).
        if let Some(units) = print_ladder_units {
            let (to_unit, unit_label) = match units {
                UnitSystem::Imperial => (1.0, "fps"),
                UnitSystem::Metric => (0.3048, "m/s"),
            };
            eprintln!(
                "BC5D Table: segment ladder (velocities in {unit_label}; paste as arguments):"
            );
            // Use shortest round-trip endpoint text: fixed decimal precision can collapse a
            // genuine narrow muzzle band into VMIN==VMAX, which the paste-back parser rejects.
            for seg in &segments {
                let display_min = seg.velocity_min * to_unit;
                let mut display_max = seg.velocity_max * to_unit;
                // Multiplication can itself map adjacent fps floats to the same metric float.
                // Preserve the strict parser contract with the smallest representable nudge.
                if display_max <= display_min {
                    display_max = display_min.next_up();
                }
                eprintln!(
                    "  --bc-segment {}:{}:{:.5}",
                    display_min,
                    display_max,
                    seg.bc_value
                );
            }
        }
        Some(segments)
    } else {
        eprintln!(
            "Warning: No BC5D table data found for caliber {:.3}\"",
            caliber
        );
        None
    }
}

// fallback_bullet_length_m moved to ballistics_engine::truing (MBA-1343).

/// How the mover Ring table column (MBA-1325) turns its raw mil angle into the display
/// value for the active `--adjustment-unit` (MBA-1355): every unit except Clicks is a
/// constant multiply; Clicks instead rounds to a whole click count via a resolved
/// [`ClickValue`].
#[derive(Debug, Clone, Copy)]
enum RingUnit {
    Factor(f64),
    Clicks(ClickValue),
}

fn run_trajectory(config: &TrajectoryConfig) -> Result<(), Box<dyn Error>> {
    // Destructure config for convenient access throughout the function
    let TrajectoryConfig {
        velocity,
        angle,
        bc,
        mass,
        diameter,
        bullet_length,
        drag_model,
        max_range,
        time_step,
        temperature,
        pressure,
        humidity,
        altitude,
        wind_speed,
        wind_direction,
        wind_vertical,
        ref wind_segments,
        output,
        full,
        plot,
        units,
        sight_height,
        bore_height,
        ignore_ground_impact,
        cant,
        use_bc_segments,
        use_cluster_bc,
        ref bc_segments_data,
        ref custom_drag_table,
        enable_magnus,
        enable_coriolis,
        enable_spin_drift,
        enable_aerodynamic_jump,
        enable_wind_shear,
        wind_shear_model,
        enable_pitch_damping,
        enable_precession,
        sample_trajectory,
        sample_interval,
        use_rk4,
        use_rk45,
        twist_rate,
        twist_right,
        latitude,
        shot_direction,
        shooting_angle,
        use_powder_sensitivity,
        powder_temp_sensitivity,
        powder_temp,
        ref powder_temp_curve,
        powder_curve_temp_c,
        target_speed,
        adjustment_unit,
        elevation_click,
        windage_click,
        ref pdf_metadata,
    } = *config;

    // Mover ring (MBA-1325): a per-point Ring column/field, additive across table/JSON/CSV,
    // enabled only when --target-speed > 0 (same units convention as `lead --target-speed`:
    // mph imperial / m/s metric). Resolved once here so every output branch below agrees.
    let target_speed_mps = UnitConverter::wind_to_metric(target_speed, units);
    let ring_enabled = target_speed_mps > 0.0;
    // Ring TABLE column angular unit honors --adjustment-unit (review fix M3). The MOA
    // conversion uses the crate's locked printed-table dial constant (MBA-724):
    // MOA_PER_UNIT_RATIO / MIL_PER_UNIT_RATIO == exactly 3.438 — deliberately NOT the
    // exact-angle 3437.7467/1000, so Ring(moa)/Ring(mil) keeps the same ratio as every
    // other MIL/MOA column pair this CLI prints (see moving_target.rs module docs).
    //
    // MBA-1355: Clicks has no fixed factor — it needs a resolved elevation graduation to
    // round the ring's mil angle to whole clicks (`RingUnit::Clicks`), rather than
    // scaling by a constant like every other unit (`RingUnit::Factor`). The Ring column
    // isn't cleanly "elevation" or "windage" (a moving-target lead ring is its own axis),
    // so it reuses the elevation click graduation, same as the PDF dope card's Drop
    // column below.
    let (ring_hdr, ring_unit) = match adjustment_unit {
        AdjustmentUnit::Mil => ("Ring(mil)", RingUnit::Factor(1.0)),
        AdjustmentUnit::Moa => (
            "Ring(moa)",
            RingUnit::Factor(
                ballistics_engine::moving_target::MOA_PER_UNIT_RATIO
                    / ballistics_engine::moving_target::MIL_PER_UNIT_RATIO,
            ),
        ),
        // MBA-1355: SMOA/IPHY are numerically identical (exact inches-per-hundred-yards);
        // only the header text differs, so they get their own labels but share the ratio.
        AdjustmentUnit::Smoa => ("Ring(smoa)", RingUnit::Factor(smoa_per_mil())),
        AdjustmentUnit::Iphy => ("Ring(iphy)", RingUnit::Factor(smoa_per_mil())),
        AdjustmentUnit::Clicks => (
            "Ring(clicks)",
            RingUnit::Clicks(elevation_click.expect(
                "adjustment_unit == Clicks implies elevation_click was resolved in main()'s \
                 Trajectory dispatch before run_trajectory was ever called",
            )),
        ),
    };

    // MBA-1135: track whether the twist is a synthesized default (shooter omitted --twist-rate)
    // so the stability summary can be honest about it rather than presenting an assumed-twist Sg
    // as a hard "STABLE" / "UNSTABLE" verdict.
    let twist_assumed = twist_rate.is_none();

    // Create ballistic inputs with all required fields
    let drag_model_enum = match drag_model {
        DragModelArg::G1 => DragModel::G1,
        DragModelArg::G7 => DragModel::G7,
    };
    let wind_direction_rad = wind_direction.to_radians();

    let inputs = BallisticInputs {
        // Core ballistics parameters
        bc_value: bc,
        bc_type: drag_model_enum,
        bullet_mass: mass,
        muzzle_velocity: velocity,
        bullet_diameter: diameter,
        bullet_length,

        // Targeting and positioning
        muzzle_angle: angle.to_radians(),
        target_distance: max_range,
        azimuth_angle: 0.0,
        shot_azimuth: shot_direction.map(|d| d.to_radians()).unwrap_or(0.0),
        shooting_angle: shooting_angle.to_radians(),
        cant_angle: cant.to_radians(),
        sight_height,
        muzzle_height: bore_height, // Bore height above ground from --bore-height CLI option
        target_height: 0.0,
        ground_threshold: if ignore_ground_impact {
            f64::NEG_INFINITY
        } else {
            0.0
        },

        // Environmental conditions
        altitude,
        temperature,
        pressure,
        humidity,
        latitude,

        // Wind conditions
        wind_speed,
        wind_angle: wind_direction_rad,

        // Bullet characteristics
        // MBA-1135: caliber/weight-aware default twist (Miller-inverse) instead of a fixed 1:12"
        // when the shooter omits --twist-rate. `twist_assumed` (above) records that this happened.
        twist_rate: twist_rate.unwrap_or_else(|| {
            ballistics_engine::stability::default_twist_inches(diameter, mass, velocity)
        }),
        is_twist_right: twist_right,
        caliber_inches: diameter / 0.0254,
        weight_grains: mass / GRAINS_TO_KG,
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
            // Convert (velocity / degree) to SI: velocity per a 1-degree DELTA.
            // The denominator must be a temperature delta, not the absolute-point
            // conversion of 1 F (which is -17.2 C) — that was MBA-963.
            UnitConverter::velocity_to_metric(powder_temp_sensitivity, units)
                / UnitConverter::temperature_delta_to_metric(1.0, units)
        } else {
            0.0
        },
        powder_temp: UnitConverter::temperature_to_metric(powder_temp, units),
        powder_temp_curve: powder_temp_curve.clone(),
        powder_curve_temp_c,
        tipoff_yaw: 0.0,
        tipoff_decay_distance: 50.0,
        // The schedule was resolved once before auto-zero and is shared by every local path.
        use_bc_segments,
        bc_segments: None,
        bc_segments_data: bc_segments_data.clone(),
        use_enhanced_spin_drift: enable_spin_drift,
        enable_aerodynamic_jump,
        use_form_factor: false,
        enable_wind_shear,
        // Forward the user's --wind-shear-model (MBA-965) instead of hardcoding
        // a single profile. "none" when shear is disabled so the engine skips it.
        wind_shear_model: if enable_wind_shear {
            wind_shear_model.as_engine_str().to_string()
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
        custom_drag_table: custom_drag_table.clone(),
    };

    // Set up wind conditions
    let wind = WindConditions {
        speed: wind_speed,
        direction: wind_direction_rad,
        vertical_speed: wind_vertical,
    };

    // Set up atmospheric conditions
    let atmosphere = AtmosphericConditions {
        temperature,
        pressure,
        humidity,
        altitude,
    };

    // Create solver
    let mut solver = TrajectorySolver::new(inputs.clone(), wind, atmosphere.clone());
    solver.set_max_range(max_range);
    solver.set_time_step(time_step);

    // Downrange-segmented wind (overrides the scalar wind when present).
    if !wind_segments.is_empty() {
        if enable_wind_shear {
            return Err(
                "--wind-segment cannot be combined with --enable-wind-shear \
                (downrange segments + altitude shear is not yet a defined model)"
                    .into(),
            );
        }
        // Note when a non-zero scalar wind is also set, since segments take precedence.
        if wind_speed != 0.0 || wind_vertical != 0.0 {
            eprintln!(
                "note: --wind-segment overrides --wind-speed/--wind-direction/--wind-vertical \
                 (scalar wind ignored)"
            );
        }
        // Warn if the segments don't cover the whole trajectory (wind is zero beyond
        // the last segment's until-distance). A 1 m epsilon avoids float-noise warnings
        // when a segment is set right at max_range.
        let coverage_m = wind_segments
            .iter()
            .map(|seg| seg.until_m)
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

    // MBA-1285: when a custom drag table is in play, warn (best-effort, stderr) if the
    // shot's Mach range runs outside the table's measured domain. The engine already
    // holds the nearest endpoint Cd out there; this just flags that as approximate.
    if let Some(table) = custom_drag_table.as_ref() {
        let (_, sos) = ballistics_engine::atmosphere::calculate_atmosphere(
            altitude,
            Some(temperature),
            Some(pressure),
            humidity,
        );
        if sos.is_finite() && sos > 0.0 {
            let muzzle_mach = velocity / sos;
            let impact_mach = result.impact_velocity / sos;
            if let (Some(&lo), Some(&hi)) =
                (table.mach_values.first(), table.mach_values.last())
            {
                if muzzle_mach > hi || impact_mach < lo {
                    eprintln!(
                        "Warning: shot Mach range [{impact_mach:.2}, {muzzle_mach:.2}] extends beyond \
                         the drag table domain [{lo:.2}, {hi:.2}]; the nearest tabulated Cd is held \
                         outside that range (approximate)."
                    );
                }
            }
        }
    }

    // Report the summary using the SAME twist the integration path actually
    // used (MBA-964). The path builds BallisticInputs with
    // `twist_rate.unwrap_or(12.0)`, so it always flies with a real twist; the
    // summary must therefore compute stability/spin-drift from that effective
    // twist rather than only when `--twist-rate` was explicitly supplied,
    // otherwise the JSON/table reports null while drift is silently applied.
    let effective_twist_rate = inputs.twist_rate;

    // Calculate stability coefficient if twist rate is provided
    let stability = if effective_twist_rate > 0.0 {
        let (resolved_temp_c, resolved_pressure_hpa) =
            ballistics_engine::atmosphere::resolve_station_conditions(
                temperature,
                pressure,
                altitude,
            );
        ballistics_engine::stability::compute_stability_coefficient(
            &inputs,
            (altitude, resolved_temp_c, resolved_pressure_hpa, 1.0),
        )
    } else {
        0.0
    };

    // Report the same canonical Litz spin drift applied to the trajectory points. The
    // empirical t^1.83 fit already reflects real spin history, so layering a retained-spin
    // multiplier on top would double-count decay and diverge from every solver path.
    let spin_drift = if enable_spin_drift && effective_twist_rate > 0.0 && stability > 0.0 {
        ballistics_engine::spin_drift::litz_drift_meters(
            stability,
            result.time_of_flight,
            twist_right,
        )
    } else {
        0.0
    };

    // Detect runs that ended because the integrator hit its 100 s time cap
    // rather than at a real impact or the requested range (MBA-969). This
    // happens for steep --ignore-ground-impact shots: ground termination is
    // disabled, so the bullet keeps "falling" past the ground and the only
    // terminator left is the time cap. Without this the summary prints the
    // capped state as a bogus ground impact (e.g. y = -12390 yd at t = 99.99 s).
    const INTEGRATION_TIME_CAP_S: f64 = 100.0;
    let reached_max_range = result.max_range >= max_range - (max_range.abs() * 1e-3).max(1e-6);
    let reached_time_cap = result
        .points
        .last()
        .map(|p| p.time >= INTEGRATION_TIME_CAP_S - 0.05)
        .unwrap_or(false)
        && !reached_max_range;

    // Format output
    match output {
        OutputFormat::Json => {
            // Honor --units like the table/csv branches do (MBA-962): convert
            // all distance/velocity/energy fields to the user's unit system and
            // record which system the numbers are in.
            let units_label = match units {
                UnitSystem::Metric => "metric",
                UnitSystem::Imperial => "imperial",
            };
            let trajectory_result = TrajectoryResult {
                units: units_label.to_string(),
                max_range: UnitConverter::distance_from_metric(result.max_range, units),
                max_height: UnitConverter::distance_from_metric(result.max_height, units),
                time_of_flight: result.time_of_flight,
                impact_velocity: UnitConverter::velocity_from_metric(result.impact_velocity, units),
                impact_energy: UnitConverter::energy_from_metric(result.impact_energy, units),
                stability_coefficient: if stability > 0.0 {
                    Some(stability)
                } else {
                    None
                },
                spin_drift: if spin_drift.abs() > 0.0001 {
                    Some(UnitConverter::distance_from_metric(spin_drift, units))
                } else {
                    None
                },
                // Pitch-damping diagnostics (MBA-966), only present when
                // --enable-pitch-damping populated them on the engine result.
                min_pitch_damping: result.min_pitch_damping,
                transonic_mach: result.transonic_mach,
                // Precession/nutation diagnostics (MBA-966), only present when
                // --enable-precession populated them. Angles stay in radians,
                // matching the table output.
                max_yaw_angle: result.max_yaw_angle,
                max_precession_angle: result.max_precession_angle,
                final_pitch_angle: result.angular_state.as_ref().map(|s| s.pitch_angle),
                final_yaw_angle: result.angular_state.as_ref().map(|s| s.yaw_angle),
                trajectory: if full {
                    result
                        .points
                        .into_iter()
                        .map(|p| {
                            // Mover ring (MBA-1325), additive: only populated when
                            // --target-speed enabled it; None/None serialize to "absent"
                            // (skip_serializing_if) so JSON without the flag is unchanged.
                            let (mover_ring_m, mover_ring_mil) = if ring_enabled {
                                let (ring_m, ring_mil) =
                                    ballistics_engine::mover_ring(target_speed_mps, p.time, p.position.x);
                                (Some(ring_m), ring_mil)
                            } else {
                                (None, None)
                            };
                            TrajectoryPoint {
                                time: p.time,
                                // Output contract is unchanged: the `x` field is lateral
                                // (drift), `z` is downrange. With McCoy internals these map
                                // to position.z (lateral) and position.x (downrange).
                                x: UnitConverter::distance_from_metric(p.position.z, units),
                                y: UnitConverter::distance_from_metric(p.position.y, units),
                                z: UnitConverter::distance_from_metric(p.position.x, units),
                                velocity: UnitConverter::velocity_from_metric(
                                    p.velocity_magnitude,
                                    units,
                                ),
                                energy: UnitConverter::energy_from_metric(p.kinetic_energy, units),
                                mover_ring_m,
                                mover_ring_mil,
                            }
                        })
                        .collect()
                } else {
                    vec![]
                },
                legend: TrajectoryLegend::for_units(units),
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
                    let header = format!(
                        "distance_{},drop_{},drift_{},velocity_{},energy_{},time_s",
                        dist_unit, defl_unit, defl_unit, vel_unit, energy_unit
                    );
                    // Mover ring (MBA-1325): extra column, header carries the unit; only
                    // emitted when --target-speed enabled it (additive, matches JSON).
                    if ring_enabled {
                        println!("{},ring_mil", header);
                    } else {
                        println!("{}", header);
                    }
                    for s in sampled {
                        let distance = UnitConverter::distance_from_metric(s.distance_m, units);
                        // Imperial drop/drift in inches (meters * 39.3701); Metric in meters.
                        let (drop, drift) = match units {
                            UnitSystem::Imperial => (s.drop_m * 39.3701, s.wind_drift_m * 39.3701),
                            UnitSystem::Metric => (s.drop_m, s.wind_drift_m),
                        };
                        let vel = UnitConverter::velocity_from_metric(s.velocity_mps, units);
                        let energy = UnitConverter::energy_from_metric(s.energy_j, units);
                        let row = format!(
                            "{:.2},{:.2},{:.2},{:.2},{:.2},{:.4}",
                            distance, drop, drift, vel, energy, s.time_s
                        );
                        if ring_enabled {
                            let (_, ring_mil) =
                                ballistics_engine::mover_ring(target_speed_mps, s.time_s, s.distance_m);
                            match ring_mil {
                                Some(mil) => println!("{},{:.3}", row, mil),
                                None => println!("{},", row),
                            }
                        } else {
                            println!("{}", row);
                        }
                    }
                } else {
                    // Output raw trajectory points (all integration steps)
                    let header = format!(
                        "time_s,x_{},y_{},z_{},velocity_{},energy_{}",
                        dist_unit, dist_unit, dist_unit, vel_unit, energy_unit
                    );
                    if ring_enabled {
                        println!("{},ring_mil", header);
                    } else {
                        println!("{}", header);
                    }
                    for p in result.points {
                        // Output contract unchanged: x column is lateral, z is downrange.
                        // McCoy internals: lateral=position.z, downrange=position.x.
                        let x = UnitConverter::distance_from_metric(p.position.z, units);
                        let y = UnitConverter::distance_from_metric(p.position.y, units);
                        let z = UnitConverter::distance_from_metric(p.position.x, units);
                        let vel = UnitConverter::velocity_from_metric(p.velocity_magnitude, units);
                        let energy = UnitConverter::energy_from_metric(p.kinetic_energy, units);
                        let row = format!(
                            "{:.4},{:.2},{:.2},{:.2},{:.2},{:.2}",
                            p.time, x, y, z, vel, energy
                        );
                        if ring_enabled {
                            // Downrange is position.x (McCoy frame) regardless of the CSV's
                            // lateral/downrange column swap above.
                            let (_, ring_mil) = ballistics_engine::mover_ring(target_speed_mps, p.time, p.position.x);
                            match ring_mil {
                                Some(mil) => println!("{},{:.3}", row, mil),
                                None => println!("{},", row),
                            }
                        } else {
                            println!("{}", row);
                        }
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
            if reached_time_cap {
                // Integration stopped at the time cap, not at an impact: the
                // "impact" velocity/energy/ToF are the capped state, not a real
                // terminal-ballistics result, so don't label them as impact.
                println!(
                    "║ Final Velocity:    {:>8.2} {:3}       ║",
                    velocity_display, velocity_unit
                );
                println!(
                    "║ Time at Cap:       {:>8.3} s          ║",
                    result.time_of_flight
                );
                println!("║ Impact:            No impact (cap)     ║");
            } else {
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
            }
            if stability > 0.0 {
                println!("╠════════════════════════════════════════╣");
                println!("║ Stability (SG):    {:>8.2}            ║", stability);
                if twist_assumed {
                    // MBA-1135: the twist was synthesized (shooter omitted --twist-rate), so this
                    // Sg is essentially the design target of the assumed twist. Do NOT present a
                    // hard STABLE/UNSTABLE verdict; a note after the box explains the assumption.
                    println!("║ Status:            {:>8}            ║", "assumed");
                } else {
                    let stability_status = if stability < 1.0 {
                        "UNSTABLE"
                    } else if stability < 1.5 {
                        "MARGINAL"
                    } else {
                        "STABLE  "
                    };
                    println!("║ Status:            {:>8}            ║", stability_status);
                }
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
                } else if min_damping > -2.0 {
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

            // MBA-1135: be explicit when the twist was synthesized rather than supplied. The Sg
            // (and any spin-drift) above are estimates computed from an assumed twist, not values
            // derived from a real barrel; do not read the Sg as an evaluated stability verdict.
            if twist_assumed && stability > 0.0 {
                let twist_note = match units {
                    UnitSystem::Imperial => format!("1:{:.1}\"", effective_twist_rate),
                    UnitSystem::Metric => format!("1:{:.1}mm", effective_twist_rate * 25.4),
                };
                println!(
                    "note: twist rate not supplied; assumed {} for the stability/spin-drift \
                     estimate — Sg shown is not an evaluated stability verdict.",
                    twist_note
                );
            }

            // Report termination. A run that ran out the integration time cap
            // (MBA-969) is NOT a ground impact even though its capped y is far
            // below zero, so flag it explicitly instead of claiming an impact.
            if reached_time_cap {
                println!();
                println!(
                    "No impact: integration reached the {:.0} s time cap (use a shorter --max-range \
                     or remove --ignore-ground-impact).",
                    INTEGRATION_TIME_CAP_S
                );
            } else if let Some(last_point) = result.points.last() {
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

                let step = if result.points.len() > 20 {
                    result.points.len() / 20
                } else {
                    1
                };

                // Mover ring (MBA-1325) is a whole extra column, so the table is rendered
                // as two distinct layouts rather than splicing a conditional cell into a
                // fixed format string — this keeps the ring_enabled == false path (the
                // default) byte-identical to the pre-MBA-1325 table.
                if ring_enabled {
                    println!("┌──────────┬──────────┬──────────┬──────────┬──────────┬──────────┐");
                    println!(
                        "│ Time (s) │  X {:5} │  Y {:5} │ Vel{:5} │Energy{:5}│ {}│",
                        dist_hdr, dist_hdr, vel_hdr, energy_hdr, ring_hdr
                    );
                    println!("├──────────┼──────────┼──────────┼──────────┼──────────┼──────────┤");

                    for (i, p) in result.points.iter().enumerate() {
                        if i % step == 0 || i == result.points.len() - 1 {
                            let x_display = UnitConverter::distance_from_metric(p.position.x, units); // X column = downrange (position.x; McCoy frame)
                            let y_display = UnitConverter::distance_from_metric(p.position.y, units);
                            let vel_display =
                                UnitConverter::velocity_from_metric(p.velocity_magnitude, units);
                            let energy_display =
                                UnitConverter::energy_from_metric(p.kinetic_energy, units);
                            let (_, ring_mil) = ballistics_engine::mover_ring(target_speed_mps, p.time, p.position.x);
                            let ring_cell = match ring_mil {
                                Some(mil) => match ring_unit {
                                    RingUnit::Factor(f) => format!("{:>8.2}", mil * f),
                                    // clicks_for(drop_yd, range_yd, click) only needs the
                                    // drop_yd/range_yd RATIO — passing (mil, 1000.0) reuses
                                    // it directly on an already-computed mil angle (ring_mil
                                    // is `ring_m / downrange_m * 1000`), same trick as
                                    // smoa_per_mil()'s factor-ratio reuse just above.
                                    RingUnit::Clicks(click) => {
                                        format!("{:>8}", clicks_for(mil, 1000.0, &click))
                                    }
                                },
                                None => format!("{:>8}", "-"),
                            };

                            println!(
                                "│ {:>8.3} │ {:>8.2} │ {:>8.2} │ {:>8.2} │ {:>8.2} │ {} │",
                                p.time, x_display, y_display, vel_display, energy_display, ring_cell
                            );
                        }
                    }
                    println!("└──────────┴──────────┴──────────┴──────────┴──────────┴──────────┘");
                } else {
                    println!("┌──────────┬──────────┬──────────┬──────────┬──────────┐");
                    println!(
                        "│ Time (s) │  X {:5} │  Y {:5} │ Vel{:5} │Energy{:5}│",
                        dist_hdr, dist_hdr, vel_hdr, energy_hdr
                    );
                    println!("├──────────┼──────────┼──────────┼──────────┼──────────┤");

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

            // Inline terminal chart (MBA-1320): two stacked charts, drop vs. range then
            // lateral drift vs. range, from the SAME per-point data the "Trajectory
            // Points:" table above prints — result.points, McCoy frame (x=downrange,
            // y=vertical, z=lateral) — converted through the same UnitConverter calls and
            // the same range_unit as that table's X/Y columns. This intentionally reuses
            // the raw, non-LOS-adjusted vertical/lateral values (the table's "Y" column),
            // NOT the sight-line-relative "drop"/"wind_drift" computed by
            // --sample-trajectory (TrajectorySample::drop_m/wind_drift_m, which needs
            // sight-height/target-height geometry this code path doesn't have and displays
            // in inches/mm rather than the range unit) — the two are different, deliberate
            // vertical conventions; don't conflate them.
            //
            // Independent of --full: result.points already holds every raw integration
            // point regardless of whether the table above printed them, so --plot works
            // standalone without needing --full too.
            if let Some(plot_style) = plot {
                if !result.points.is_empty() {
                    let drop_label = format!("drop ({})", range_unit);
                    let drop_points: Vec<(f64, f64)> = result
                        .points
                        .iter()
                        .map(|p| {
                            (
                                UnitConverter::distance_from_metric(p.position.x, units),
                                UnitConverter::distance_from_metric(p.position.y, units),
                            )
                        })
                        .collect();
                    println!();
                    println!("Drop vs Range:");
                    println!(
                        "{}",
                        terminal_plot::render_chart(
                            &[(drop_label.as_str(), drop_points.as_slice())],
                            72,
                            12,
                            plot_style.as_canvas_style(),
                        )
                    );

                    let drift_label = format!("drift ({})", range_unit);
                    let drift_points: Vec<(f64, f64)> = result
                        .points
                        .iter()
                        .map(|p| {
                            (
                                UnitConverter::distance_from_metric(p.position.x, units),
                                UnitConverter::distance_from_metric(p.position.z, units),
                            )
                        })
                        .collect();
                    println!();
                    println!("Lateral Drift vs Range:");
                    println!(
                        "{}",
                        terminal_plot::render_chart(
                            &[(drift_label.as_str(), drift_points.as_slice())],
                            72,
                            12,
                            plot_style.as_canvas_style(),
                        )
                    );
                }
            }
        }

        #[cfg(feature = "pdf")]
        OutputFormat::Pdf => {
            // PDF output requires metadata and output file
            let pdf_meta = pdf_metadata.as_ref().ok_or(
                "PDF output requires metadata (use --target-speed, --powder, --bullet-name, etc.)",
            )?;
            let output_path = pdf_meta
                .output_file
                .as_ref()
                .ok_or("PDF output requires --output-file")?;

            // Get sampled trajectory points (required for dope card)
            let sampled = result
                .sampled_points
                .as_ref()
                .ok_or("PDF output requires --sample-trajectory flag for trajectory data")?;

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
                solver_mode: if cfg!(feature = "online") {
                    "online".to_string()
                } else {
                    "offline".to_string()
                },
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
                unit_label: match pdf_meta.adjustment_unit {
                    AdjustmentUnit::Mil => "MIL".to_string(),
                    AdjustmentUnit::Moa => "MOA".to_string(),
                    AdjustmentUnit::Smoa => "SMOA".to_string(),
                    AdjustmentUnit::Iphy => "IPHY".to_string(),
                    AdjustmentUnit::Clicks => "CLICKS".to_string(),
                },
            };

            // Convert sampled trajectory to dope card rows
            // Filter from zero range onwards (typically 100 yards)
            let rows: Vec<DopeCardRow> = sampled
                .iter()
                .filter(|s| {
                    s.distance_m >= UnitConverter::distance_to_metric(100.0, UnitSystem::Imperial)
                })
                .map(|s| {
                    // Convert distance to yards for range
                    let range_yd =
                        UnitConverter::distance_from_metric(s.distance_m, UnitSystem::Imperial);
                    // Convert drop to yards (s.drop_m is already in meters, positive = below line of sight)
                    let drop_yd = s.drop_m / 0.9144; // meters to yards
                                                     // Convert drift to yards
                    let drift_yd = s.wind_drift_m / 0.9144;
                    // Lead: how far the target moves during the bullet's flight, as an
                    // angular hold. Movement (yd) = speed (yd/s) * time; 1760/3600 = mph->yd/s.
                    // Shared moving-target math (MBA-1287): mph → m/s, perpendicular
                    // (90°) crossing — lead_m = speed·TOF exactly, as before; the yd
                    // conversion chain differs from the old 1760/3600 factor by ≤1 ulp.
                    let lead_m = ballistics_engine::lead_from_tof(
                        pdf_meta.target_speed_mph * 0.44704,
                        90.0,
                        s.time_s,
                        s.distance_m,
                    )
                    .lead_m;
                    let lead_yd = lead_m / 0.9144;
                    let unit = pdf_meta.adjustment_unit;

                    DopeCardRow {
                        range_yd: range_yd.round() as u32,
                        // Drop: positive = dial up (bullet below LOS). drop_m is already
                        // positive-below-LOS (sample_trajectory: los_y - y_interp), matching the
                        // come-up / range tables, so do NOT negate it (this column was sign-flipped).
                        // Rendered in MIL/MOA/SMOA/IPHY (MBA-724) or whole clicks (MBA-1355) per
                        // --adjustment-unit; drop uses the elevation graduation, wind/lead the
                        // windage one.
                        drop_adj: adjustment_display(drop_yd, range_yd, unit, elevation_click),
                        // Wind: positive = dial right for wind from right
                        wind_adj: adjustment_display(drift_yd, range_yd, unit, windage_click),
                        // Lead for a moving target
                        lead_adj: adjustment_display(lead_yd, range_yd, unit, windage_click),
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
            eprintln!(
                "  {} ranges from {} to {} yards",
                rows.len(),
                rows.first().map(|r| r.range_yd).unwrap_or(0),
                rows.last().map(|r| r.range_yd).unwrap_or(0)
            );
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
                println!("zero_angle,{:.4},degrees", response.zero_angle.to_degrees());
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
                println!(
                    "║ BC Confidence:     {:>8.1}%           ║",
                    confidence * 100.0
                );
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
            eprintln!(
                "Hint: Use local calculation (without --online) for PDF dope card generation."
            );
        }
    }

    Ok(())
}

// ============================================================================
// WEZ (Weapon Employment Zone) sweep mode -- MBA-1317
//
// The compute path (target-size parsing, the seeded sweep, and variance attribution) lives in
// ballistics_engine::wez (MBA-1343 Phase B) so non-CLI front ends (e.g. the WASM terminal) can
// reuse it; only the three renderers (summary table / statistics CSV / full JSON) remain here.
// ============================================================================

/// Dominant-bucket label shared by the summary-table and statistics-CSV renderers.
fn wez_dominant_label(row: &WezRow) -> String {
    if row.attribution_unavailable {
        "n/a".to_string()
    } else {
        row.dominant_error_source
            .map(|b| b.to_string())
            .unwrap_or_else(|| "none".to_string())
    }
}

/// `-o full`: the entire [`WezResult`] as pretty JSON (the 0.25.0 JSON contract).
fn print_wez_full(result: &WezResult) -> Result<(), Box<dyn Error>> {
    println!("{}", serde_json::to_string_pretty(result)?);
    Ok(())
}

/// `-o statistics`: one CSV row per range step.
fn print_wez_statistics(result: &WezResult, units: UnitSystem) {
    println!("range,p_hit,dominant_error_source,wind_call_share,mv_sd_share,other_share");
    for row in &result.rows {
        println!(
            "{:.2},{:.4},{},{:.4},{:.4},{:.4}",
            UnitConverter::distance_from_metric(row.range_m, units),
            row.p_hit,
            wez_dominant_label(row),
            row.wind_call_share,
            row.mv_sd_share,
            row.other_share
        );
    }
}

/// Default `-o summary`: the human-readable sweep table.
fn print_wez_summary(result: &WezResult, units: UnitSystem) {
    let dist_unit = match units {
        UnitSystem::Imperial => "yd",
        UnitSystem::Metric => "m",
    };
    let wind_unit = match units {
        UnitSystem::Imperial => "mph",
        UnitSystem::Metric => "m/s",
    };
    println!(
        "WEZ sweep: {} sims/step, wind call {:.2} {wind_unit} + wind std {:.2} {wind_unit} (quadrature) = {:.2} {wind_unit} effective",
        result.num_sims_per_step,
        UnitConverter::wind_from_metric(result.wind_call_error_mps, units),
        UnitConverter::wind_from_metric(result.wind_speed_std_mps, units),
        UnitConverter::wind_from_metric(result.combined_wind_speed_std_mps, units),
    );
    println!(
        "┌────────────┬──────────┬───────────────┬───────────┬───────────┬───────────┐"
    );
    println!(
        "│ Range ({dist_unit:>3}) │  P(hit)  │ Dominant      │ Wind call │  MV SD    │ Other/grp │"
    );
    println!(
        "├────────────┼──────────┼───────────────┼───────────┼───────────┼───────────┤"
    );
    for row in &result.rows {
        println!(
            "│ {:>10.1} │ {:>7.1}% │ {:<13} │ {:>8.1}% │ {:>8.1}% │ {:>8.1}% │",
            UnitConverter::distance_from_metric(row.range_m, units),
            row.p_hit * 100.0,
            wez_dominant_label(row),
            row.wind_call_share * 100.0,
            row.mv_sd_share * 100.0,
            row.other_share * 100.0,
        );
    }
    println!(
        "└────────────┴──────────┴───────────────┴───────────┴───────────┴───────────┘"
    );
}

#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments mirror the stable Monte Carlo CLI command shape (MBA-1317)"
)]
fn run_monte_carlo_wez(
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
    wind_direction_std: f64,
    wind_speed: f64,
    wind_direction: f64,
    wind_vertical: f64,
    wind_call_error: f64,
    target_size: TargetSizeMetric,
    wez_start: f64,
    wez_end: f64,
    wez_step: f64,
    custom_drag_table: Option<ballistics_engine::drag::DragTable>,
    cant: f64,
    output: MonteCarloOutput,
    units: UnitSystem,
) -> Result<(), Box<dyn Error>> {
    let result = compute_wez(
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
        wind_direction_std,
        wind_speed,
        wind_direction,
        wind_vertical,
        wind_call_error,
        target_size,
        wez_start,
        wez_end,
        wez_step,
        // The native `monte-carlo` subcommand has no --drag-model flag, so the sweep
        // always runs the G1 default `BallisticInputs` has always used (the WASM
        // terminal's monte-carlo, which does expose --drag-model, passes it through).
        DragModel::G1,
        custom_drag_table,
        cant,
    )?;

    match output {
        MonteCarloOutput::Full => print_wez_full(&result)?,
        MonteCarloOutput::Statistics => print_wez_statistics(&result, units),
        MonteCarloOutput::Summary => print_wez_summary(&result, units),
    }

    Ok(())
}

#[cfg(test)]
mod wez_tests {
    use super::*;

    // The compute-path tests (parse_target_size, step function, monotonicity, variance shares)
    // moved to ballistics_engine::wez alongside the code they exercise (MBA-1343 Phase B). This
    // one stays: it pins the CLI wrapper `run_monte_carlo_wez` (compute + print dispatch)
    // end-to-end.

    // ---- Solver-invalid input is a clean error, not a panic (MBA-1317) ---------------------

    /// `-v 0` passes clap's own `0.0..=6000` velocity range check (clap permits the lower bound)
    /// but fails `TrajectorySolver::solve()`'s stricter "> 0" gate. Before the fix,
    /// `wez_solve_target_plane`'s `.expect(...)` on that failure turned this into a panic
    /// (exit 101 with a backtrace) instead of the clean `Err` the base (non-WEZ) `monte-carlo`
    /// command already returns for the identical input. This must return `Err`, never panic.
    #[test]
    fn zero_muzzle_velocity_is_a_clean_error_not_a_panic() {
        let result = run_monte_carlo_wez(
            /* velocity */ 0.0,
            /* angle */ 0.0,
            /* bc */ 0.475,
            /* mass */ 0.010_886,
            /* diameter */ 0.007_82,
            /* num_sims */ 10,
            /* velocity_std */ 0.0,
            /* angle_std */ 0.0,
            /* bc_std */ 0.0,
            /* wind_std */ 0.0,
            /* wind_direction_std */ 0.0,
            /* wind_speed */ 0.0,
            /* wind_direction */ 0.0,
            /* wind_vertical */ 0.0,
            /* wind_call_error */ 0.0,
            TargetSizeMetric::Rect {
                width_m: 0.4572,
                height_m: 0.762,
            },
            /* wez_start */ 200.0,
            /* wez_end */ 300.0,
            /* wez_step */ 100.0,
            /* custom_drag_table */ None,
            /* cant */ 0.0,
            MonteCarloOutput::Summary,
            UnitSystem::Metric,
        );
        let err = result.expect_err("a zero muzzle velocity must not solve successfully");
        assert!(
            err.to_string().contains("muzzle_velocity"),
            "expected a muzzle_velocity validation error, got: {err}"
        );
    }
}

#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments mirror the stable Monte Carlo CLI command shape"
)]
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
    wind_direction_std: f64,
    wind_speed: f64,
    wind_direction: f64,
    wind_vertical: f64,
    target_distance: Option<f64>,
    target_radius: f64,
    custom_drag_table: Option<ballistics_engine::drag::DragTable>,
    cant: f64,
    output: MonteCarloOutput,
) -> Result<(), Box<dyn Error>> {
    // Create base inputs. MBA-967: use the same bore-height/ground convention as the
    // `trajectory` subcommand (standard 1.5 m bore height; this helper works in metric) so each
    // simulation stops at a realistic ground impact instead of flying to the integrator's range
    // cap. Without this, "Mean Range" reports the ~1000 m cap rather than the ground-impact range.
    let bore_height_metric = 1.5_f64;
    // The base inputs are what each Monte Carlo sample perturbs, so the drag table must be set
    // here for the deck to apply to every simulated shot (MBA-1285).
    let base_inputs = BallisticInputs {
        muzzle_velocity: velocity,
        muzzle_angle: angle.to_radians(),
        bc_value: bc,
        bullet_mass: mass,
        bullet_diameter: diameter,
        muzzle_height: bore_height_metric,
        ground_threshold: 0.0,
        custom_drag_table,
        cant_angle: cant.to_radians(),
        ..Default::default()
    };

    // Create base wind conditions
    let base_wind = WindConditions {
        speed: wind_speed,
        direction: wind_direction.to_radians(),
        vertical_speed: wind_vertical,
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
    let results = ballistics_engine::run_monte_carlo_with_wind_and_direction_std_dev(
        base_inputs,
        base_wind,
        mc_params,
        wind_direction_std.to_radians(),
    )?;

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

    // CEP is conditional on actually reaching the target plane. Samples that fall short stay in
    // all result vectors (and in the hit-probability denominator), but their finite -1e9 marker is
    // not a physical target-plane impact and must not become the median (MBA-1159).
    let cep = results.target_plane_cep();
    let target_shortfall_fraction = results.target_shortfall_fraction();

    // Calculate hit probability if target distance specified
    // MBA-971: hit probability is the fraction of samples landing within the hit radius of the
    // point of aim at the target plane (shared MonteCarloResults::hit_probability), matching the
    // FFI. The old range-precision notion (ground-impact range within 1 m of the target) reported
    // 0% for any target short of the impact range.
    let hit_probability = target_distance.map(|_target| results.hit_probability(target_radius));

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
            if let Some(cep) = cep {
                println!("║ CEP (arrivals):    {:>8.2} m          ║", cep);
            } else {
                println!("║ CEP (arrivals):         N/A            ║");
            }
            println!(
                "║ Target Shortfall:  {:>8.1} %          ║",
                target_shortfall_fraction * 100.0
            );
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
                target_shortfall_fraction,
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

#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments mirror the stable zero-calculation CLI command shape"
)]
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
    custom_drag_table: Option<ballistics_engine::drag::DragTable>,
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
        custom_drag_table,
        ..Default::default()
    };

    // Set up atmospheric conditions
    let atmosphere = AtmosphericConditions {
        temperature,
        pressure,
        humidity,
        altitude,
    };

    // The bullet starts at bore height. A same-elevation zero crosses the horizontal line of
    // sight at y = sight_height; target_height is an offset at the target.
    let los_target_height = sight_height + target_height;

    // Calculate zero angle with atmospheric conditions
    let zero_angle = ballistics_engine::calculate_zero_angle_with_conditions(
        inputs.clone(),
        target_distance,
        los_target_height,
        WindConditions::default(),
        atmosphere.clone(),
    )?;

    // Calculate trajectory at zero angle to get additional info
    let mut zeroed_inputs = inputs;
    zeroed_inputs.muzzle_angle = zero_angle;

    let mut solver = TrajectorySolver::new(zeroed_inputs, Default::default(), atmosphere);
    solver.set_max_range(target_distance * 3.0);
    let trajectory = solver.solve()?;

    // `zero_angle` is the bore angle above horizontal. The sight line runs from
    // (0, sight_height) to (target_distance, sight_height + target_height), so its slope is
    // target_height / target_distance: the sight-height translation cancels between endpoints.
    let sight_adjustment_moa = zero_angle.to_degrees() * 60.0
        - (target_height / target_distance * 3437.75);
    // The bullet starts below the sight line by `sight_height`, which can already exceed the
    // 5 cm lower bound at the muzzle. Ignore that initial offset and report only a lower-bound
    // exit after the trajectory has first entered the point-blank band.
    let mut entered_point_blank_band = false;
    let point_blank_range = trajectory
        .points
        .iter()
        .find_map(|p| {
            let los_y = sight_height + target_height * (p.position.x / target_distance);
            let height_from_los = p.position.y - los_y;

            if height_from_los >= -0.05 {
                entered_point_blank_band = true;
                None
            } else if entered_point_blank_band {
                Some(p.position.x)
            } else {
                None
            }
        })
        .unwrap_or(trajectory.max_range);

    match output {
        OutputFormat::Json => {
            let units_label = match units {
                UnitSystem::Metric => "metric",
                UnitSystem::Imperial => "imperial",
            };
            let result = serde_json::json!({
                "units": units_label,
                "zero_angle_degrees": zero_angle.to_degrees(),
                "zero_angle_moa": zero_angle.to_degrees() * 60.0,
                "zero_angle_mrad": zero_angle * 1000.0,
                "sight_adjustment_moa": sight_adjustment_moa,
                "max_ordinate": UnitConverter::distance_from_metric(trajectory.max_height, units),
                "point_blank_range": UnitConverter::distance_from_metric(point_blank_range, units),
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
    let weight_grains = mass / GRAINS_TO_KG;
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

/// Parse a `"d,v;d,v;..."` data string into `(f64, f64)` pairs, tolerating surrounding
/// quotes and whitespace. Errors on any malformed pair.
fn parse_data_pairs(s: &str) -> Result<Vec<(f64, f64)>, Box<dyn Error>> {
    let cleaned = s.trim().trim_matches('\'').trim_matches('"');
    let mut out = Vec::new();
    for pair in cleaned.split(';') {
        let pair = pair.trim();
        if pair.is_empty() {
            continue;
        }
        let parts: Vec<&str> = pair.split(',').collect();
        if parts.len() != 2 {
            return Err(
                format!("Malformed data pair '{}': expected \"distance,value\".", pair).into(),
            );
        }
        let d: f64 = parts[0]
            .trim()
            .parse()
            .map_err(|_| format!("Invalid distance '{}'.", parts[0].trim()))?;
        let v: f64 = parts[1]
            .trim()
            .parse()
            .map_err(|_| format!("Invalid value '{}'.", parts[1].trim()))?;
        out.push((d, v));
    }
    if out.is_empty() {
        return Err("No valid data pairs found.".into());
    }
    Ok(out)
}

/// Parse the `--drag-model` selector into the list of drag models to estimate.
fn parse_drag_models(s: &str) -> Result<Vec<DragModel>, Box<dyn Error>> {
    match s.trim().to_lowercase().as_str() {
        "g1" => Ok(vec![DragModel::G1]),
        "g7" => Ok(vec![DragModel::G7]),
        "both" | "all" | "g1,g7" | "g1g7" => Ok(vec![DragModel::G1, DragModel::G7]),
        other => Err(format!("Unknown --drag-model '{}'; use g1, g7, or both.", other).into()),
    }
}

/// Estimate BC for every (drag model × available data basis) combination and print the
/// results. `drop_points` are `(m, m)`, `vel_points` are `(m, m/s)`.
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments preserve the existing BC-estimation compatibility helper"
)]
fn run_bc_estimation_multi(
    velocity: f64,
    mass: f64,
    diameter: f64,
    drop_points: &[(f64, f64)],
    vel_points: &[(f64, f64)],
    models: &[DragModel],
    atmosphere: AtmosphericConditions,
    zero_range: Option<f64>,
    sight_height: f64,
    units: UnitSystem,
    output: OutputFormat,
) -> Result<(), Box<dyn Error>> {
    use ballistics_engine::{estimate_bc_fit, BcFitMode};

    // The data bases we can fit against, in a stable order (drop first, then velocity).
    struct Basis<'a> {
        mode: BcFitMode,
        points: &'a [(f64, f64)],
    }
    let mut bases: Vec<Basis> = Vec::new();
    if !drop_points.is_empty() {
        bases.push(Basis {
            mode: BcFitMode::Drop,
            points: drop_points,
        });
    }
    if !vel_points.is_empty() {
        bases.push(Basis {
            mode: BcFitMode::Velocity,
            points: vel_points,
        });
    }

    struct Variant {
        model: DragModel,
        mode: BcFitMode,
        bc: f64,
        rms_user: f64,
        rms_unit: &'static str,
        n: usize,
        at_bound: bool,
    }
    let mut variants: Vec<Variant> = Vec::new();
    for &model in models {
        for basis in &bases {
            // Zero range only shapes a drop fit; a velocity fit is frame-independent.
            let zr = match basis.mode {
                BcFitMode::Drop => zero_range,
                BcFitMode::Velocity => None,
            };
            let est = estimate_bc_fit(
                velocity,
                mass,
                diameter,
                basis.points,
                model,
                basis.mode,
                atmosphere.clone(),
                zr,
                sight_height,
            )?;
            // Convert the RMS residual back to user units for a readable fit-quality column.
            let (rms_user, rms_unit) = match basis.mode {
                BcFitMode::Drop => match units {
                    UnitSystem::Imperial => (est.rms_error / 0.0254, "in"),
                    UnitSystem::Metric => (est.rms_error, "m"),
                },
                BcFitMode::Velocity => match units {
                    UnitSystem::Imperial => (est.rms_error / 0.3048, "fps"),
                    UnitSystem::Metric => (est.rms_error, "m/s"),
                },
            };
            variants.push(Variant {
                model,
                mode: basis.mode,
                bc: est.bc,
                rms_user,
                rms_unit,
                n: basis.points.len(),
                at_bound: est.at_bound,
            });
        }
    }

    let model_name = |m: DragModel| match m {
        DragModel::G7 => "G7",
        DragModel::G1 => "G1",
        _ => "G?",
    };
    let basis_name = |mode: BcFitMode| match mode {
        BcFitMode::Drop => "drop",
        BcFitMode::Velocity => "velocity",
    };

    match output {
        OutputFormat::Json => {
            let vs: Vec<_> = variants
                .iter()
                .map(|v| {
                    serde_json::json!({
                        "drag_model": model_name(v.model),
                        "fit_basis": basis_name(v.mode),
                        "estimated_bc": (v.bc * 1000.0).round() / 1000.0,
                        "fit_rms": (v.rms_user * 1000.0).round() / 1000.0,
                        "fit_rms_unit": v.rms_unit,
                        "n_points": v.n,
                        "reliable": !v.at_bound,
                    })
                })
                .collect();
            let result = serde_json::json!({
                "velocity": velocity,
                "mass": mass,
                "diameter": diameter,
                "variants": vs,
            });
            println!("{}", serde_json::to_string_pretty(&result)?);
        }

        OutputFormat::Csv => {
            println!("drag_model,fit_basis,estimated_bc,fit_rms,fit_rms_unit,n_points,reliable");
            for v in &variants {
                println!(
                    "{},{},{:.4},{:.3},{},{},{}",
                    model_name(v.model),
                    basis_name(v.mode),
                    v.bc,
                    v.rms_user,
                    v.rms_unit,
                    v.n,
                    !v.at_bound
                );
            }
        }

        OutputFormat::Table => {
            println!("BC Estimation");
            println!(
                "  {:<6} {:<20} {:>12}   {:<10} ",
                "Model", "Fit basis", "Estimated BC", "Fit RMS"
            );
            println!("  {:-<6} {:-<20} {:->12}   {:-<10}", "", "", "", "");
            for v in &variants {
                println!(
                    "  {:<6} {:<20} {:>12.3}   {:>6.2} {:<4}{}",
                    model_name(v.model),
                    format!("{} ({} pts)", basis_name(v.mode), v.n),
                    v.bc,
                    v.rms_user,
                    v.rms_unit,
                    if v.at_bound { " ⚠ UNRELIABLE (hit BC limit)" } else { "" }
                );
            }
            if variants.iter().any(|v| v.at_bound) {
                println!();
                println!("  ⚠ One or more fits ran to the BC search limit — the data did not");
                println!("    determine a real value. Add more (and longer-range) points, and");
                println!("    check --zero-range / --temperature / --pressure match the data.");
            }
        }

        OutputFormat::Pdf => {
            eprintln!("Error: PDF output is not supported for BC estimation.");
            eprintln!("Hint: Use --output json, --output csv, or --output table instead.");
        }
    }

    Ok(())
}

/// Display velocity truing results in the specified format
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments preserve the existing velocity-truing display helper"
)]
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
            println!(
                "{:.1},{},{},{},{},{},{:.4},{:.2},{}",
                effective_vel,
                vel_unit,
                velocity_adjustment
                    .map(|v| {
                        let adj = match units {
                            UnitSystem::Imperial => v,
                            UnitSystem::Metric => v * 0.3048,
                        };
                        format!("{:.1}", adj)
                    })
                    .unwrap_or_default(),
                adjustment_percent
                    .map(|v| format!("{:.2}", v))
                    .unwrap_or_default(),
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
            println!(
                "║  Effective Muzzle Velocity: {:>8.1} {:>4}                 ║",
                effective_vel, vel_unit
            );
            if let Some(adj) = velocity_adjustment {
                let adj_display = match units {
                    UnitSystem::Imperial => adj,
                    UnitSystem::Metric => adj * 0.3048,
                };
                println!(
                    "║  Adjustment from Chrono:    {:>+8.1} {:>4}                 ║",
                    adj_display, vel_unit
                );
                if let Some(pct) = adjustment_percent {
                    println!(
                        "║  Adjustment Percentage:     {:>+8.2}%                      ║",
                        pct
                    );
                }
            }
            println!("╠════════════════════════════════════════════════════════════╣");
            println!(
                "║  Confidence:                {:>8}                        ║",
                confidence
            );
            println!(
                "║  Iterations:                {:>8}                        ║",
                iterations
            );
            println!(
                "║  Final Error:               {:>8.4} MIL                  ║",
                final_error_mil
            );
            println!(
                "║  Calculated Drop:           {:>8.2} MIL                  ║",
                calculated_drop_mil
            );
            println!(
                "║  Measured Drop:             {:>8.2} MIL                  ║",
                measured_drop
            );
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

// calculate_drop_at_velocity, TrueVelocityLocalResult, and calculate_true_velocity_local
// moved to ballistics_engine::truing (MBA-1343 Phase C), completing the truing lib surface
// (single-observation binary search + multi-observation joint fit) so the WASM terminal can
// reuse the exact compute paths. solve_trajectory_drop moved earlier (Phase A).

// ============================================================================
// MBA-1316: multi-observation joint MV + BC calibration (truing v2)
// ============================================================================
//
// The compute core (constants, TruingObservation, TruingForwardModel, the
// fitters, identifiability gates and MultiTruingReport) moved to
// ballistics_engine::truing (MBA-1343) so the WASM terminal can reuse it.
// Only the thin CLI wrapper and the rendering stay here.

/// Parse an uncertainty-aware additional observation.  The optional third
/// component is a per-observation 1-sigma override; otherwise `default_sigma`
/// is used.  This parser is reached only when the opt-in uncertainty surface is
/// active, leaving the legacy RANGE:DROP parser and its error contract intact.
fn parse_uncertain_truing_observation(
    token: &str,
    units: UnitSystem,
    default_sigma: f64,
) -> Result<(TruingObservation, f64), String> {
    let parts: Vec<&str> = token.split(':').collect();
    if !(2..=3).contains(&parts.len()) {
        return Err(format!(
            "invalid --observed '{token}': expected RANGE:DROP or RANGE:DROP:SIGMA"
        ));
    }
    let base = format!("{}:{}", parts[0], parts[1]);
    let observation = parse_truing_observation(&base, units)?;
    let sigma = if parts.len() == 3 {
        parts[2]
            .trim()
            .parse::<f64>()
            .map_err(|_| format!("invalid --observed sigma '{}' in '{token}'", parts[2]))?
    } else {
        default_sigma
    };
    if !sigma.is_finite() || sigma <= 0.0 {
        return Err(format!(
            "invalid --observed sigma in '{token}': sigma must be positive and finite"
        ));
    }
    Ok((observation, sigma))
}

/// Parse an explicit normal prior token as `(mean, one_sigma)`.
fn parse_normal_prior(token: &str, flag: &str) -> Result<(f64, f64), String> {
    let (mean, sigma) = token
        .split_once(':')
        .ok_or_else(|| format!("{flag} expects MEAN:SIGMA"))?;
    let mean = mean
        .trim()
        .parse::<f64>()
        .map_err(|_| format!("{flag} has an invalid mean in '{token}'"))?;
    let sigma = sigma
        .trim()
        .parse::<f64>()
        .map_err(|_| format!("{flag} has an invalid sigma in '{token}'"))?;
    if !mean.is_finite() {
        return Err(format!("{flag} mean must be finite"));
    }
    if !sigma.is_finite() || sigma <= 0.0 {
        return Err(format!("{flag} sigma must be positive and finite"));
    }
    Ok((mean, sigma))
}

/// Expand the explicit `START:END:STEP` planner grid.  The grid is never
/// inferred: the end is included only when it lies on a generated step (within
/// floating-point tolerance), and the hard cap keeps accidental CLI input from
/// launching an unbounded number of trajectory solves.
fn parse_truing_range_grid(token: &str) -> Result<Vec<f64>, String> {
    let parts: Vec<&str> = token.split(':').collect();
    if parts.len() != 3 {
        return Err(format!(
            "--range-grid expects START:END:STEP (got '{token}')"
        ));
    }
    let mut values = [0.0_f64; 3];
    for (i, part) in parts.iter().enumerate() {
        values[i] = part
            .trim()
            .parse::<f64>()
            .map_err(|_| format!("--range-grid contains an invalid number in '{token}'"))?;
    }
    let [start, end, step] = values;
    if !start.is_finite() || !end.is_finite() || !step.is_finite() {
        return Err("--range-grid values must be finite".to_string());
    }
    if start <= 0.0 || end < start || step <= 0.0 {
        return Err("--range-grid requires 0 < START <= END and STEP > 0".to_string());
    }
    let count = ((end - start) / step + 1e-12).floor() as usize + 1;
    if count > 1_000 {
        return Err(format!(
            "--range-grid expands to {count} candidates; maximum is 1000"
        ));
    }
    Ok((0..count).map(|i| start + i as f64 * step).collect())
}


/// Orchestrate the multi-observation joint MV+BC calibration and print the
/// result. `measured_drop`/`range_yd` are the primary observation (drop already
/// in `drop_unit`); `observed` holds the additional `RANGE:DROP` tokens.
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments mirror the stable true-velocity CLI command shape"
)]
fn run_multi_observation_truing(
    measured_drop: f64,
    range_yd: f64,
    observed: &[String],
    drop_unit: DropUnit,
    bc_input: f64,
    drag_model: DragModelArg,
    mass_gr: f64,
    diameter_in: f64,
    zero_yd: f64,
    sight_in: f64,
    temp_f: f64,
    press_inhg: f64,
    humidity: f64,
    alt_ft: f64,
    bc_segments: &Option<Vec<BCSegmentData>>,
    chrono_fps: Option<f64>,
    units: UnitSystem,
    output: OutputFormat,
) -> Result<(), Box<dyn Error>> {
    // Assemble the observation set: primary (--range/--measured-drop) first,
    // then every --observed token (the typed core takes parsed observations;
    // MBA-1343 review).
    let mut observations: Vec<TruingObservation> = Vec::with_capacity(observed.len() + 1);
    observations.push(TruingObservation {
        range_yd,
        drop: measured_drop,
    });
    for token in observed {
        observations.push(parse_truing_observation(token, units)?);
    }

    // Pre-validate before announcing the fit: the core validates again (cheaply),
    // but running the same check here first keeps the historical stderr contract —
    // an invalid set (e.g. duplicate ranges) errors WITHOUT a "Fitting ..."
    // progress line, exactly as it always has.
    validate_truing_observations(&observations)?;

    eprintln!(
        "Fitting {} observations (joint MV+BC calibration)...",
        observations.len()
    );

    let report = run_multi_observation_truing_core(
        &observations,
        drop_unit,
        bc_input,
        drag_model,
        mass_gr,
        diameter_in,
        zero_yd,
        sight_in,
        temp_f,
        press_inhg,
        humidity,
        alt_ft,
        bc_segments,
    )?;
    display_multi_truing_result(&report, drop_unit, units, chrono_fps, output);
    Ok(())
}

// rms_at and truing_quality_line moved to ballistics_engine::truing (MBA-1343).

/// Render a multi-observation truing report as table / JSON / CSV.
fn display_multi_truing_result(
    report: &MultiTruingReport,
    drop_unit: DropUnit,
    units: UnitSystem,
    chrono_fps: Option<f64>,
    output: OutputFormat,
) {
    let vel_unit = match units {
        UnitSystem::Imperial => "fps",
        UnitSystem::Metric => "m/s",
    };
    let range_unit = match units {
        UnitSystem::Imperial => "yd",
        UnitSystem::Metric => "m",
    };
    let drop_label = drop_unit.label();
    let mv_display = match units {
        UnitSystem::Imperial => report.fitted_mv_fps,
        UnitSystem::Metric => report.fitted_mv_fps * 0.3048,
    };
    let range_display = |range_yd: f64| match units {
        UnitSystem::Imperial => range_yd,
        UnitSystem::Metric => range_yd * 0.9144,
    };
    // Chrono comparison (chrono_fps is already in fps).
    let (adj_display, adj_pct) = match chrono_fps {
        Some(c) => {
            let adj_fps = report.fitted_mv_fps - c;
            let adj_disp = match units {
                UnitSystem::Imperial => adj_fps,
                UnitSystem::Metric => adj_fps * 0.3048,
            };
            let pct = if c != 0.0 { adj_fps / c * 100.0 } else { 0.0 };
            (Some(adj_disp), Some(pct))
        }
        None => (None, None),
    };

    match output {
        OutputFormat::Json => {
            let obs_json: Vec<serde_json::Value> = report
                .observations
                .iter()
                .enumerate()
                .map(|(i, o)| {
                    serde_json::json!({
                        format!("range_{range_unit}"): range_display(o.range_yd),
                        format!("observed_drop_{drop_label}"): o.drop,
                        format!("predicted_drop_{drop_label}"): report.predicted[i],
                        format!("residual_{drop_label}"): report.residuals[i],
                    })
                })
                .collect();

            let mut json_output = serde_json::json!({
                "mode": if report.bc_fitted { "joint_mv_bc" } else { "mv_only" },
                "fitted_muzzle_velocity": mv_display,
                "velocity_unit": vel_unit,
                "bc_fitted": report.bc_fitted,
                "fitted_bc": report.fitted_bc,
                "input_bc": report.bc_input,
                "observations": obs_json,
                format!("rms_residual_{drop_label}"): report.rms,
                "iterations": report.iterations,
                "converged": report.converged,
                "bc_sensitivity_ratio": report.sensitivity_ratio,
                "condition_number": if report.condition_number.is_finite() {
                    serde_json::json!(report.condition_number)
                } else {
                    serde_json::Value::Null
                },
                "quality": report.quality,
                "legend": {
                    "units": {
                        "range": range_unit,
                        "drop": drop_label,
                        "velocity": vel_unit,
                    },
                },
            });
            if !report.reason.is_empty() {
                json_output["bc_hold_reason"] = serde_json::json!(report.reason);
            }
            if let Some(adj) = adj_display {
                json_output["velocity_adjustment"] = serde_json::json!(adj);
            }
            if let Some(pct) = adj_pct {
                json_output["adjustment_percent"] = serde_json::json!(pct);
            }
            match serde_json::to_string_pretty(&json_output) {
                Ok(s) => println!("{s}"),
                Err(e) => eprintln!("Error serializing JSON: {e}"),
            }
        }
        OutputFormat::Csv => {
            println!("range_{range_unit},observed_drop_{drop_label},predicted_drop_{drop_label},residual_{drop_label}");
            for (i, o) in report.observations.iter().enumerate() {
                println!(
                    "{:.1},{:.4},{:.4},{:+.4}",
                    range_display(o.range_yd),
                    o.drop,
                    report.predicted[i],
                    report.residuals[i]
                );
            }
            println!();
            println!("fitted_muzzle_velocity_{vel_unit},bc_fitted,fitted_bc,input_bc,rms_residual_{drop_label},iterations,converged,bc_sensitivity_ratio,condition_number");
            println!(
                "{:.1},{},{:.4},{:.4},{:.4},{},{},{:.4},{}",
                mv_display,
                report.bc_fitted,
                report.fitted_bc,
                report.bc_input,
                report.rms,
                report.iterations,
                report.converged,
                report.sensitivity_ratio,
                if report.condition_number.is_finite() {
                    format!("{:.1}", report.condition_number)
                } else {
                    "inf".to_string()
                },
            );
        }
        OutputFormat::Table => {
            println!();
            println!("=== VELOCITY + BC TRUING (multi-observation) ===");
            println!();
            println!(
                "  Fitted muzzle velocity: {:>9.1} {}",
                mv_display, vel_unit
            );
            if report.bc_fitted {
                println!(
                    "  Fitted BC:              {:>9.4}  (input {:.4})",
                    report.fitted_bc, report.bc_input
                );
            } else {
                println!(
                    "  BC:                     {:>9.4}  (held; not fitted)",
                    report.fitted_bc
                );
            }
            if let Some(adj) = adj_display {
                println!("  Adjustment from chrono: {:>+9.1} {}", adj, vel_unit);
                if let Some(pct) = adj_pct {
                    println!("  Adjustment percentage:  {:>+9.2}%", pct);
                }
            }
            println!();
            println!(
                "  {:>10}  {:>14}  {:>14}  {:>12}",
                format!("Range ({range_unit})"),
                format!("Observed ({drop_label})"),
                format!("Predicted ({drop_label})"),
                format!("Resid ({drop_label})"),
            );
            println!("  {}", "-".repeat(56));
            for (i, o) in report.observations.iter().enumerate() {
                println!(
                    "  {:>10.1}  {:>14.3}  {:>14.3}  {:>+12.3}",
                    range_display(o.range_yd),
                    o.drop,
                    report.predicted[i],
                    report.residuals[i]
                );
            }
            println!("  {}", "-".repeat(56));
            println!(
                "  RMS residual: {:.3} {}   |   iterations: {}{}",
                report.rms,
                drop_label,
                report.iterations,
                if report.converged {
                    ""
                } else {
                    " (not fully converged)"
                }
            );
            println!();
            println!("  {}", report.quality);
            if !report.reason.is_empty() {
                println!("  Note: {}", report.reason);
            }
            println!(
                "  Diagnostics: BC sensitivity ratio {:.4}, conditioning {}",
                report.sensitivity_ratio,
                if report.condition_number.is_finite() {
                    format!("{:.0}", report.condition_number)
                } else {
                    "inf".to_string()
                }
            );
            println!();
        }
        OutputFormat::Pdf => {
            eprintln!("Error: PDF output is not supported for velocity truing results.");
            eprintln!("Hint: Use --output json, --output csv, or --output table instead.");
        }
    }
}

/// Convert a Rust enum's stable CamelCase variant name into the snake_case code
/// used by the uncertainty JSON/CSV contracts.  Keeping this generic makes a
/// newly added warning/failure code visible instead of silently dropping it.
fn enum_variant_code(value: &impl std::fmt::Debug) -> String {
    let name = format!("{value:?}");
    let mut code = String::with_capacity(name.len() + 4);
    for (index, ch) in name.chars().enumerate() {
        if ch.is_ascii_uppercase() && index > 0 {
            code.push('_');
        }
        code.push(ch.to_ascii_lowercase());
    }
    code
}

/// Render MBA-1346's deterministic pre-shoot observation-range design.
fn display_truing_experiment_plan(
    plan: &TruingExperimentPlanV1,
    units: UnitSystem,
    output: OutputFormat,
) -> Result<(), Box<dyn Error>> {
    let range_scale = match units {
        UnitSystem::Imperial => 1.0,
        UnitSystem::Metric => 0.9144,
    };
    let range_unit = match units {
        UnitSystem::Imperial => "yd",
        UnitSystem::Metric => "m",
    };
    let range = |yards: f64| yards * range_scale;
    let mode = match plan.mode {
        TruingPlanModeV1::JointMvBc => "joint_mv_bc",
        TruingPlanModeV1::MvOnly => "mv_only",
    };

    match output {
        OutputFormat::Json => {
            let selected: Vec<serde_json::Value> = plan
                .selected_stations
                .iter()
                .map(|station| {
                    serde_json::json!({
                        "input_index": station.input_index,
                        format!("range_{range_unit}"): range(station.range_yd),
                        format!("predicted_drop_{}", plan.measurement_drop_unit): station.predicted_drop,
                        "scaled_mv_sensitivity": station.scaled_mv_sensitivity,
                        "scaled_bc_sensitivity": station.scaled_bc_sensitivity,
                        "leave_one_out_information_gain_nats": station.leave_one_out_information_gain_nats,
                    })
                })
                .collect();
            let rejected: Vec<serde_json::Value> = plan
                .rejected_candidates
                .iter()
                .map(|candidate| {
                    serde_json::json!({
                        "input_index": candidate.input_index,
                        format!("range_{range_unit}"): if candidate.range_yd.is_finite() {
                            serde_json::json!(range(candidate.range_yd))
                        } else {
                            serde_json::Value::Null
                        },
                        "reason": enum_variant_code(&candidate.reason),
                        "detail": candidate.detail,
                    })
                })
                .collect();
            let info = &plan.information;
            let value = serde_json::json!({
                "schema_version": "truing-experiment-plan-v1",
                "mode": mode,
                "requested_observation_count": plan.requested_observation_count,
                format!("minimum_separation_{range_unit}"): range(plan.minimum_separation_yd),
                "measurement_model": {
                    "sigma_1sd": plan.measurement_sigma_1sd,
                    "drop_unit": plan.measurement_drop_unit,
                    "assumption": "independent Gaussian impact-reading error at every station",
                },
                "selected_stations": selected,
                format!("unselected_candidate_ranges_{range_unit}"): plan
                    .unselected_candidate_ranges_yd
                    .iter()
                    .map(|value| range(*value))
                    .collect::<Vec<_>>(),
                "rejected_candidates": rejected,
                "information": {
                    "bc_sensitivity_ratio": info.sensitivity_ratio,
                    "condition_number": info.condition_number,
                    "minimum_singular_value": info.minimum_singular_value,
                    "maximum_singular_value": info.maximum_singular_value,
                    "weak_axis_fractional_sigma_1sd": info.weak_axis_fractional_sigma_1sd,
                    "log_determinant": info.log_determinant,
                    "expected_information_gain_nats": info.expected_information_gain_nats,
                },
                "optimizer": {
                    "search_strategy": enum_variant_code(&plan.search_strategy),
                    "eligible_candidate_count": plan.eligible_candidate_count,
                    "raw_combination_count": plan.raw_combination_count,
                    "evaluated_design_count": plan.evaluated_design_count,
                },
                "warnings": plan.warnings,
                "legend": {
                    "range_unit": range_unit,
                    "drop_unit": plan.measurement_drop_unit,
                    "scaled_sensitivity": "predicted drop change in declared observation sigmas for a unit fractional parameter change",
                    "information_gain": "0.5 log det(I + F), where I is the disclosed identity reference information in fractional MV/BC coordinates; design score only",
                },
            });
            println!("{}", serde_json::to_string_pretty(&value)?);
        }
        OutputFormat::Csv => {
            println!("plan_mode,measurement_sigma_1sd,measurement_drop_unit,requested_observation_count,minimum_separation_{range_unit},bc_sensitivity_ratio,condition_number,minimum_singular_value,maximum_singular_value,weak_axis_fractional_sigma_1sd,log_determinant,expected_information_gain_nats");
            println!(
                "{mode},{:.8},{},{},{:.8},{:.10},{},{:.10},{:.10},{},{},{:.10}",
                plan.measurement_sigma_1sd,
                plan.measurement_drop_unit,
                plan.requested_observation_count,
                range(plan.minimum_separation_yd),
                plan.information.sensitivity_ratio,
                plan.information
                    .condition_number
                    .map(|value| format!("{value:.10}"))
                    .unwrap_or_default(),
                plan.information.minimum_singular_value,
                plan.information.maximum_singular_value,
                plan.information
                    .weak_axis_fractional_sigma_1sd
                    .map(|value| format!("{value:.10}"))
                    .unwrap_or_default(),
                plan.information
                    .log_determinant
                    .map(|value| format!("{value:.10}"))
                    .unwrap_or_default(),
                plan.information.expected_information_gain_nats,
            );
            println!();
            println!("range_{range_unit},status,predicted_drop_{},scaled_mv_sensitivity,scaled_bc_sensitivity,information_gain_nats,detail", plan.measurement_drop_unit);
            for station in &plan.selected_stations {
                println!(
                    "{:.6},selected,{:.8},{:.8},{:.8},{:.8},",
                    range(station.range_yd),
                    station.predicted_drop,
                    station.scaled_mv_sensitivity,
                    station.scaled_bc_sensitivity,
                    station.leave_one_out_information_gain_nats,
                );
            }
            for candidate in &plan.unselected_candidate_ranges_yd {
                println!("{:.6},eligible_not_selected,,,,,", range(*candidate));
            }
            for candidate in &plan.rejected_candidates {
                let display_range = if candidate.range_yd.is_finite() {
                    range(candidate.range_yd).to_string()
                } else {
                    String::new()
                };
                println!(
                    "{},rejected:{},,,,,{}",
                    display_range,
                    enum_variant_code(&candidate.reason),
                    candidate.detail.replace(',', ";")
                );
            }
            if !plan.warnings.is_empty() {
                println!();
                println!("warning");
                for warning in &plan.warnings {
                    println!("{}", warning.replace(',', ";").replace('\n', " "));
                }
            }
        }
        OutputFormat::Table => {
            println!();
            println!("=== TRUING EXPERIMENT PLAN ===");
            println!(
                "  Recommendation: {}",
                if plan.mode == TruingPlanModeV1::JointMvBc {
                    "joint MV + BC"
                } else {
                    "MV only (available ranges do not identify BC)"
                }
            );
            println!(
                "  Measurement assumption: independent 1-sigma reading error = {:.4} {}",
                plan.measurement_sigma_1sd, plan.measurement_drop_unit
            );
            println!(
                "  Selected ranges: {} {range_unit}",
                plan.selected_stations
                    .iter()
                    .map(|station| format!("{:.1}", range(station.range_yd)))
                    .collect::<Vec<_>>()
                    .join(", ")
            );
            println!(
                "  Constraints: {} stations, minimum separation {:.1} {range_unit}",
                plan.requested_observation_count,
                range(plan.minimum_separation_yd)
            );
            println!();
            println!(
                "  BC sensitivity ratio: {:.5}",
                plan.information.sensitivity_ratio
            );
            println!(
                "  Condition number: {}",
                plan.information
                    .condition_number
                    .map(|value| format!("{value:.2}"))
                    .unwrap_or_else(|| "unbounded".to_string())
            );
            println!(
                "  Singular values: min {:.5}, max {:.5}; weak-axis fractional 1-sigma: {}",
                plan.information.minimum_singular_value,
                plan.information.maximum_singular_value,
                plan.information
                    .weak_axis_fractional_sigma_1sd
                    .map(|value| format!("{value:.4}"))
                    .unwrap_or_else(|| "unbounded".to_string())
            );
            println!();
            println!(
                "  {:>10}  {:>12}  {:>13}  {:>13}  {:>12}",
                format!("Range ({range_unit})"),
                format!("Drop ({})", plan.measurement_drop_unit),
                "MV sensitivity",
                "BC sensitivity",
                "Info (nats)",
            );
            for station in &plan.selected_stations {
                println!(
                    "  {:>10.1}  {:>12.4}  {:>13.5}  {:>13.5}  {:>12.5}",
                    range(station.range_yd),
                    station.predicted_drop,
                    station.scaled_mv_sensitivity,
                    station.scaled_bc_sensitivity,
                    station.leave_one_out_information_gain_nats,
                );
            }
            println!();
            println!(
                "  Search: {} ({} of {} candidate sets evaluated)",
                enum_variant_code(&plan.search_strategy),
                plan.evaluated_design_count,
                plan.raw_combination_count,
            );
            if !plan.rejected_candidates.is_empty() {
                println!("  Rejected candidates:");
                for candidate in &plan.rejected_candidates {
                    let display_range = if candidate.range_yd.is_finite() {
                        format!("{:.1} {range_unit}", range(candidate.range_yd))
                    } else {
                        "non-finite range".to_string()
                    };
                    println!(
                        "    - {display_range}: {} ({})",
                        enum_variant_code(&candidate.reason),
                        candidate.detail
                    );
                }
            }
            if !plan.warnings.is_empty() {
                println!("  Notes:");
                for warning in &plan.warnings {
                    println!("    - {warning}");
                }
            }
            println!();
        }
        OutputFormat::Pdf => {
            return Err("PDF output is not supported for truing experiment plans".into());
        }
    }
    Ok(())
}

/// Render the opt-in MBA-1353 weighted MAP and local Gaussian approximation.
/// The legacy renderer above is deliberately not reused or modified: callers
/// reach this function only after explicitly declaring observation sigma.
fn display_uncertainty_truing_result(
    report: &UncertaintyTruingReportV1,
    units: UnitSystem,
    chrono_fps: Option<f64>,
    output: OutputFormat,
) -> Result<(), Box<dyn Error>> {
    let velocity_scale = match units {
        UnitSystem::Imperial => 1.0,
        UnitSystem::Metric => 0.3048,
    };
    let velocity_unit = match units {
        UnitSystem::Imperial => "fps",
        UnitSystem::Metric => "m/s",
    };
    let range_scale = match units {
        UnitSystem::Imperial => 1.0,
        UnitSystem::Metric => 0.9144,
    };
    let range_unit = match units {
        UnitSystem::Imperial => "yd",
        UnitSystem::Metric => "m",
    };
    let drop_unit = report.drop_unit.label();
    let velocity = |fps: f64| fps * velocity_scale;
    let range = |yards: f64| yards * range_scale;

    match output {
        OutputFormat::Json => {
            let observations: Vec<serde_json::Value> = report
                .observations
                .iter()
                .map(|observation| {
                    serde_json::json!({
                        format!("range_{range_unit}"): range(observation.range_yd),
                        format!("observed_drop_{drop_unit}"): observation.observed_drop,
                        format!("sigma_1sd_{drop_unit}"): observation.sigma,
                        format!("predicted_drop_{drop_unit}"): observation.predicted_drop,
                        format!("residual_{drop_unit}"): observation.residual,
                        "standardized_residual": observation.standardized_residual,
                    })
                })
                .collect();
            let priors = serde_json::json!({
                "muzzle_velocity": report.priors.muzzle_velocity_fps.map(|prior| serde_json::json!({
                    "mean": velocity(prior.mean),
                    "sigma_1sd": velocity(prior.sigma),
                    "unit": velocity_unit,
                })),
                "ballistic_coefficient": report.priors.ballistic_coefficient.map(|prior| serde_json::json!({
                    "mean": prior.mean,
                    "sigma_1sd": prior.sigma,
                })),
            });
            let approximation = match &report.approximation {
                TruingApproximationV1::Available(approximation) => {
                    let mv_interval = approximation.muzzle_velocity_interval_95;
                    let bc_interval = approximation.ballistic_coefficient_interval_95;
                    serde_json::json!({
                        "available": true,
                        "method": "laplace_gauss_newton_known_sigma",
                        "interval_probability": mv_interval.probability,
                        "muzzle_velocity": {
                            "estimate": velocity(mv_interval.estimate),
                            "standard_deviation": velocity(mv_interval.standard_deviation),
                            "lower": velocity(mv_interval.lower),
                            "upper": velocity(mv_interval.upper),
                            "unit": velocity_unit,
                        },
                        "ballistic_coefficient": {
                            "estimate": bc_interval.estimate,
                            "standard_deviation": bc_interval.standard_deviation,
                            "lower": bc_interval.lower,
                            "upper": bc_interval.upper,
                        },
                        "covariance": {
                            format!("mv_variance_{velocity_unit}2"): approximation.covariance.mv_variance_fps2 * velocity_scale * velocity_scale,
                            format!("mv_bc_covariance_{velocity_unit}"): approximation.covariance.mv_bc_covariance_fps * velocity_scale,
                            "bc_variance": approximation.covariance.bc_variance,
                        },
                        "mv_bc_correlation": approximation.mv_bc_correlation,
                        "scaled_information_condition_number": approximation.scaled_information_condition_number,
                    })
                }
                TruingApproximationV1::Unavailable(failure) => serde_json::json!({
                    "available": false,
                    "method": "laplace_gauss_newton_known_sigma",
                    "failure": {
                        "code": enum_variant_code(&failure.code),
                        "message": failure.message,
                    },
                }),
            };
            let predictive_bands: Vec<serde_json::Value> = report
                .predictive_bands
                .iter()
                .map(|band| {
                    let interval = |value: ballistics_engine::truing_uncertainty::GaussianIntervalV1| {
                        serde_json::json!({
                            "standard_deviation": value.standard_deviation,
                            "lower": value.lower,
                            "upper": value.upper,
                            "probability": value.probability,
                        })
                    };
                    serde_json::json!({
                        format!("range_{range_unit}"): range(band.range_yd),
                        format!("predicted_drop_{drop_unit}"): band.predicted_drop,
                        "latent_model_drop_interval": band.latent_interval_95.map(interval),
                        "future_observation_interval": band.future_observation_interval_95.map(interval),
                    })
                })
                .collect();
            let warnings: Vec<serde_json::Value> = report
                .warnings
                .iter()
                .map(|warning| {
                    serde_json::json!({
                        "code": enum_variant_code(&warning.code),
                        "message": warning.message,
                    })
                })
                .collect();
            let diagnostics = &report.diagnostics;
            let mut value = serde_json::json!({
                "schema_version": "truing-uncertainty-v1",
                "mode": "joint_mv_bc_map",
                "fitted_muzzle_velocity": velocity(report.map_muzzle_velocity_fps),
                "velocity_unit": velocity_unit,
                "fitted_bc": report.map_ballistic_coefficient,
                "iterations": report.iterations,
                "converged": report.converged,
                "observations": observations,
                "priors": priors,
                "diagnostics": {
                    "chi_square": diagnostics.chi_square,
                    "prior_penalty": diagnostics.prior_penalty,
                    "penalized_chi_square": diagnostics.penalized_chi_square,
                    "effective_parameter_count": diagnostics.effective_parameter_count,
                    "effective_degrees_of_freedom": diagnostics.effective_degrees_of_freedom,
                    "reduced_chi_square": diagnostics.reduced_chi_square,
                    "bc_sensitivity_ratio": diagnostics.bc_sensitivity_ratio,
                    "data_condition_number": if diagnostics.data_condition_number.is_finite() {
                        serde_json::json!(diagnostics.data_condition_number)
                    } else {
                        serde_json::Value::Null
                    },
                    "map_scaled_gradient_inf_norm": diagnostics.map_scaled_gradient_inf_norm,
                    "map_convergence_criterion": diagnostics
                        .map_convergence_criterion
                        .map(|criterion| enum_variant_code(&criterion)),
                    "map_objective_poll_radius": diagnostics.map_objective_poll_radius,
                    "map_max_objective_poll_improvement": diagnostics
                        .map_max_objective_poll_improvement,
                    "map_objective_poll_evaluations": diagnostics
                        .map_objective_poll_evaluations,
                },
                "approximation": approximation,
                "predictive_bands": predictive_bands,
                "warnings": warnings,
                "legend": {
                    "units": {
                        "range": range_unit,
                        "drop": drop_unit,
                        "velocity": velocity_unit,
                    },
                    "observation_sigma": "absolute known one-standard-deviation reading error",
                    "latent_model_drop_interval": "parameter uncertainty only",
                    "future_observation_interval": "parameter uncertainty plus explicitly declared future reading error",
                },
            });
            if let Some(chrono) = chrono_fps {
                value["velocity_adjustment"] =
                    serde_json::json!(velocity(report.map_muzzle_velocity_fps - chrono));
                value["adjustment_percent"] = if chrono != 0.0 {
                    serde_json::json!(
                        (report.map_muzzle_velocity_fps - chrono) / chrono * 100.0
                    )
                } else {
                    serde_json::Value::Null
                };
            }
            println!("{}", serde_json::to_string_pretty(&value)?);
        }
        OutputFormat::Csv => {
            let convergence_criterion = report
                .diagnostics
                .map_convergence_criterion
                .map(|criterion| enum_variant_code(&criterion))
                .unwrap_or_default();
            println!("range_{range_unit},observed_drop_{drop_unit},sigma_1sd_{drop_unit},predicted_drop_{drop_unit},residual_{drop_unit},standardized_residual");
            for observation in &report.observations {
                println!(
                    "{:.6},{:.8},{:.8},{:.8},{:.8},{:.8}",
                    range(observation.range_yd),
                    observation.observed_drop,
                    observation.sigma,
                    observation.predicted_drop,
                    observation.residual,
                    observation.standardized_residual,
                );
            }
            println!();
            println!("fitted_muzzle_velocity_{velocity_unit},fitted_bc,chi_square,prior_penalty,penalized_chi_square,effective_parameter_count,effective_degrees_of_freedom,reduced_chi_square,bc_sensitivity_ratio,data_condition_number,converged,map_convergence_criterion,map_scaled_gradient_inf_norm,map_objective_poll_radius,map_max_objective_poll_improvement,map_objective_poll_evaluations");
            println!(
                "{:.8},{:.10},{:.8},{:.8},{:.8},{},{},{},{:.8},{},{},{},{:.10},{},{},{}",
                velocity(report.map_muzzle_velocity_fps),
                report.map_ballistic_coefficient,
                report.diagnostics.chi_square,
                report.diagnostics.prior_penalty,
                report.diagnostics.penalized_chi_square,
                report
                    .diagnostics
                    .effective_parameter_count
                    .map(|value| format!("{value:.8}"))
                    .unwrap_or_default(),
                report
                    .diagnostics
                    .effective_degrees_of_freedom
                    .map(|value| format!("{value:.8}"))
                    .unwrap_or_default(),
                report
                    .diagnostics
                    .reduced_chi_square
                    .map(|value| format!("{value:.8}"))
                    .unwrap_or_default(),
                report.diagnostics.bc_sensitivity_ratio,
                if report.diagnostics.data_condition_number.is_finite() {
                    format!("{:.8}", report.diagnostics.data_condition_number)
                } else {
                    "inf".to_string()
                },
                report.converged,
                convergence_criterion,
                report.diagnostics.map_scaled_gradient_inf_norm,
                report
                    .diagnostics
                    .map_objective_poll_radius
                    .map(|value| format!("{value:.10}"))
                    .unwrap_or_default(),
                report
                    .diagnostics
                    .map_max_objective_poll_improvement
                    .map(|value| format!("{value:.10}"))
                    .unwrap_or_default(),
                report.diagnostics.map_objective_poll_evaluations,
            );
            println!();
            println!("prior_parameter,mean,sigma_1sd,unit");
            if let Some(prior) = report.priors.muzzle_velocity_fps {
                println!(
                    "muzzle_velocity,{:.8},{:.8},{velocity_unit}",
                    velocity(prior.mean),
                    velocity(prior.sigma)
                );
            }
            if let Some(prior) = report.priors.ballistic_coefficient {
                println!(
                    "ballistic_coefficient,{:.10},{:.10},dimensionless",
                    prior.mean, prior.sigma
                );
            }
            println!();
            match &report.approximation {
                TruingApproximationV1::Available(approximation) => {
                    println!("parameter,estimate,standard_deviation,interval_95_lower,interval_95_upper,interval_probability,unit");
                    let mv = approximation.muzzle_velocity_interval_95;
                    let bc = approximation.ballistic_coefficient_interval_95;
                    println!(
                        "muzzle_velocity,{:.8},{:.8},{:.8},{:.8},{:.6},{velocity_unit}",
                        velocity(mv.estimate),
                        velocity(mv.standard_deviation),
                        velocity(mv.lower),
                        velocity(mv.upper),
                        mv.probability,
                    );
                    println!(
                        "ballistic_coefficient,{:.10},{:.10},{:.10},{:.10},{:.6},dimensionless",
                        bc.estimate,
                        bc.standard_deviation,
                        bc.lower,
                        bc.upper,
                        bc.probability,
                    );
                    println!();
                    println!("covariance_component,value");
                    println!(
                        "mv_variance_{velocity_unit}2,{:.12}",
                        approximation.covariance.mv_variance_fps2
                            * velocity_scale
                            * velocity_scale
                    );
                    println!(
                        "mv_bc_covariance_{velocity_unit},{:.12}",
                        approximation.covariance.mv_bc_covariance_fps * velocity_scale
                    );
                    println!(
                        "bc_variance,{:.12}",
                        approximation.covariance.bc_variance
                    );
                    println!(
                        "mv_bc_correlation,{:.12}",
                        approximation.mv_bc_correlation
                    );
                    println!(
                        "scaled_information_condition_number,{:.12}",
                        approximation.scaled_information_condition_number
                    );
                }
                TruingApproximationV1::Unavailable(failure) => {
                    println!("approximation_failure_code,message");
                    println!(
                        "{},{}",
                        enum_variant_code(&failure.code),
                        failure.message.replace(',', ";").replace('\n', " ")
                    );
                }
            }
            if !report.predictive_bands.is_empty() {
                println!();
                println!("prediction_range_{range_unit},predicted_drop_{drop_unit},latent_lower_{drop_unit},latent_upper_{drop_unit},future_lower_{drop_unit},future_upper_{drop_unit}");
                for band in &report.predictive_bands {
                    println!(
                        "{:.6},{:.8},{},{},{},{}",
                        range(band.range_yd),
                        band.predicted_drop,
                        band.latent_interval_95
                            .map(|interval| format!("{:.8}", interval.lower))
                            .unwrap_or_default(),
                        band.latent_interval_95
                            .map(|interval| format!("{:.8}", interval.upper))
                            .unwrap_or_default(),
                        band.future_observation_interval_95
                            .map(|interval| format!("{:.8}", interval.lower))
                            .unwrap_or_default(),
                        band.future_observation_interval_95
                            .map(|interval| format!("{:.8}", interval.upper))
                            .unwrap_or_default(),
                    );
                }
            }
            if !report.warnings.is_empty() {
                println!();
                println!("warning_code,message");
                for warning in &report.warnings {
                    println!(
                        "{},{}",
                        enum_variant_code(&warning.code),
                        warning.message.replace(',', ";").replace('\n', " ")
                    );
                }
            }
        }
        OutputFormat::Table => {
            println!();
            println!("=== UNCERTAINTY-AWARE MV + BC TRUING ===");
            println!("  Method: weighted joint MAP + local Gaussian approximation");
            println!("  Observation errors: absolute known 1-sigma ({drop_unit})");
            println!();
            println!(
                "  MAP muzzle velocity: {:>10.2} {velocity_unit}",
                velocity(report.map_muzzle_velocity_fps)
            );
            println!(
                "  MAP ballistic coefficient: {:>8.5}",
                report.map_ballistic_coefficient
            );
            match &report.approximation {
                TruingApproximationV1::Available(approximation) => {
                    let mv = approximation.muzzle_velocity_interval_95;
                    let bc = approximation.ballistic_coefficient_interval_95;
                    println!(
                        "  MV 95% interval: [{:.2}, {:.2}] {velocity_unit} (SD {:.2})",
                        velocity(mv.lower),
                        velocity(mv.upper),
                        velocity(mv.standard_deviation)
                    );
                    println!(
                        "  BC 95% interval: [{:.5}, {:.5}] (SD {:.5})",
                        bc.lower, bc.upper, bc.standard_deviation
                    );
                    println!(
                        "  MV/BC correlation: {:+.5}",
                        approximation.mv_bc_correlation
                    );
                }
                TruingApproximationV1::Unavailable(failure) => {
                    println!(
                        "  Gaussian approximation: unavailable ({})",
                        enum_variant_code(&failure.code)
                    );
                    println!("    {}", failure.message);
                }
            }
            println!();
            println!(
                "  {:>10}  {:>12}  {:>9}  {:>12}  {:>10}",
                format!("Range ({range_unit})"),
                format!("Observed ({drop_unit})"),
                "Sigma",
                format!("Predicted ({drop_unit})"),
                "Std resid",
            );
            for observation in &report.observations {
                println!(
                    "  {:>10.1}  {:>12.4}  {:>9.4}  {:>12.4}  {:>+10.3}",
                    range(observation.range_yd),
                    observation.observed_drop,
                    observation.sigma,
                    observation.predicted_drop,
                    observation.standardized_residual,
                );
            }
            println!();
            println!(
                "  Chi-square: {:.3}; effective DOF: {}; reduced chi-square: {}",
                report.diagnostics.chi_square,
                report
                    .diagnostics
                    .effective_degrees_of_freedom
                    .map(|value| format!("{value:.3}"))
                    .unwrap_or_else(|| "unavailable".to_string()),
                report
                    .diagnostics
                    .reduced_chi_square
                    .map(|value| format!("{value:.3}"))
                    .unwrap_or_else(|| "unavailable".to_string()),
            );
            if !report.predictive_bands.is_empty() {
                println!();
                println!("  Predictive drop bands (95%):");
                for band in &report.predictive_bands {
                    let latent = band
                        .latent_interval_95
                        .map(|interval| format!("[{:.4}, {:.4}]", interval.lower, interval.upper))
                        .unwrap_or_else(|| "unavailable".to_string());
                    let future = band
                        .future_observation_interval_95
                        .map(|interval| format!("; future [{:.4}, {:.4}]", interval.lower, interval.upper))
                        .unwrap_or_default();
                    println!(
                        "    {:>8.1} {range_unit}: mean {:.4} {drop_unit}; model {latent}{future}",
                        range(band.range_yd),
                        band.predicted_drop,
                    );
                }
            }
            if !report.warnings.is_empty() {
                println!();
                println!("  Warnings:");
                for warning in &report.warnings {
                    println!(
                        "    - {}: {}",
                        enum_variant_code(&warning.code),
                        warning.message
                    );
                }
            }
            println!();
        }
        OutputFormat::Pdf => {
            return Err("PDF output is not supported for uncertainty-aware truing".into());
        }
    }
    Ok(())
}

// ============================================================================
// MPBR, Come-Ups, Wind Card Handler Functions
// ============================================================================

/// Shared: build BallisticInputs + atmosphere + wind from common parameters (all in metric)
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments preserve the shared CLI trajectory compatibility helper"
)]
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
    // MBA-1323 Phase 2: saved-profile velocity-BC segments / Mach-Cd drag curve, already
    // resolved to engine shapes by the caller (bc_segments_from_profile /
    // drag_table_from_profile). `None` for every caller that does not (yet) source these from
    // a saved profile — see handle_come_ups/handle_lead for the callers that do.
    bc_segments_data: Option<Vec<BCSegmentData>>,
    custom_drag_table: Option<ballistics_engine::drag::DragTable>,
) -> (BallisticInputs, WindConditions, AtmosphericConditions) {
    let drag_model_enum = match drag_model {
        DragModelArg::G1 => DragModel::G1,
        DragModelArg::G7 => DragModel::G7,
    };
    let wind_direction_rad = wind_direction.to_radians();
    let use_bc_segments = bc_segments_data.is_some();

    let inputs = BallisticInputs {
        bc_value: bc,
        bc_type: drag_model_enum,
        bullet_mass: mass,
        muzzle_velocity: velocity,
        bullet_diameter: diameter,
        bullet_length: fallback_bullet_length_m(diameter, mass),
        muzzle_angle: 0.0,
        target_distance: max_range,
        sight_height,
        altitude,
        temperature,
        pressure,
        humidity,
        wind_speed,
        wind_angle: wind_direction_rad,
        use_rk4: true,           // Required for non-Euler solver
        use_adaptive_rk45: true, // Use RK45 adaptive (default solver)
        enable_trajectory_sampling: true,
        sample_interval,
        caliber_inches: diameter / 0.0254,
        weight_grains: mass / GRAINS_TO_KG,
        twist_rate: 12.0,
        is_twist_right: true,
        use_bc_segments,
        bc_segments_data,
        custom_drag_table,
        ..Default::default()
    };

    // wind_direction enters from the CLI in degrees; both engine structures use radians.
    let wind = WindConditions {
        speed: wind_speed,
        direction: wind_direction_rad,
        ..Default::default()
    };

    let atmosphere = AtmosphericConditions {
        temperature,
        pressure,
        humidity,
        altitude,
    };

    (inputs, wind, atmosphere)
}

/// Run a trajectory and return sampled points at the given zero angle
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments preserve the shared sampled-trajectory compatibility helper"
)]
fn run_sampled_trajectory(
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
    zero_angle_rad: f64,
    // MBA-1323 Phase 2: see build_trajectory_components's doc comment on these two.
    bc_segments_data: Option<Vec<BCSegmentData>>,
    custom_drag_table: Option<ballistics_engine::drag::DragTable>,
) -> Result<Vec<trajectory_sampling::TrajectorySample>, Box<dyn Error>> {
    let (mut inputs, wind, atmosphere) = build_trajectory_components(
        velocity,
        bc,
        mass,
        diameter,
        drag_model,
        sight_height,
        temperature,
        pressure,
        humidity,
        altitude,
        wind_speed,
        wind_direction,
        max_range,
        sample_interval,
        bc_segments_data,
        custom_drag_table,
    );
    inputs.muzzle_angle = zero_angle_rad;

    let mut solver = TrajectorySolver::new(inputs, wind, atmosphere);
    solver.set_max_range(max_range);
    solver.set_time_step(0.001);
    let result = solver.solve()?;

    Ok(result.sampled_points.unwrap_or_default())
}

/// Resolve bullet parameters: CLI arg overrides profile value
fn resolve_param(
    cli_val: Option<f64>,
    profile: &Option<ProfileData>,
    getter: fn(&ProfileData) -> f64,
) -> Option<f64> {
    cli_val.or_else(|| profile.as_ref().map(getter))
}

/// MPBR handler
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments mirror the stable MPBR CLI command shape"
)]
fn handle_mpbr(
    velocity: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    drag_model: DragModelArg,
    vital_zone: f64,   // in user units (inches or cm)
    sight_height: f64, // in user units
    temperature: f64,
    pressure: f64,
    humidity: f64,
    altitude: f64,
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
            bullet_length: fallback_bullet_length_m(diameter_m, mass_kg),
            sight_height: sight_height_m,
            use_rk4: true,
            ..Default::default()
        };

        let atmosphere = AtmosphericConditions {
            temperature: temperature_c,
            pressure: pressure_hpa,
            humidity,
            altitude: altitude_m,
        };

        let zero_angle = match ballistics_engine::calculate_zero_angle_with_conditions(
            zero_inputs,
            test_zero_m,
            sight_height_m,
            WindConditions::default(),
            atmosphere,
        ) {
            Ok(a) => a,
            Err(_) => {
                zero_high_m = test_zero_m;
                continue;
            }
        };

        // Run trajectory with this zero, sampling every ~1 yard
        let samples = match run_sampled_trajectory(
            velocity_m,
            bc,
            mass_kg,
            diameter_m,
            drag_model,
            sight_height_m,
            temperature_c,
            pressure_hpa,
            humidity,
            altitude_m,
            0.0,
            0.0,
            test_zero_m * 1.5, // max range past zero
            UnitConverter::distance_to_metric(1.0, UnitSystem::Imperial), // ~1 yd sample interval
            zero_angle,
            // MPBR does not yet consume saved-profile bc_segments/drag_curve (MBA-1323
            // Phase 2 follow-up) — see CLI_USAGE.md's a7p import section.
            None,
            None,
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
        bullet_length: fallback_bullet_length_m(diameter_m, mass_kg),
        sight_height: sight_height_m,
        use_rk4: true,
        ..Default::default()
    };

    let atmosphere = AtmosphericConditions {
        temperature: temperature_c,
        pressure: pressure_hpa,
        humidity,
        altitude: altitude_m,
    };

    let final_zero_angle = ballistics_engine::calculate_zero_angle_with_conditions(
        final_inputs,
        best_zero_m,
        sight_height_m,
        WindConditions::default(),
        atmosphere,
    )?;

    let final_samples = run_sampled_trajectory(
        velocity_m,
        bc,
        mass_kg,
        diameter_m,
        drag_model,
        sight_height_m,
        temperature_c,
        pressure_hpa,
        humidity,
        altitude_m,
        0.0,
        0.0,
        best_zero_m * 2.0,
        UnitConverter::distance_to_metric(1.0, UnitSystem::Imperial),
        final_zero_angle,
        None,
        None,
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
        let prev_drop = if sign_flip {
            -final_samples[i - 1].drop_m
        } else {
            final_samples[i - 1].drop_m
        };
        let curr_drop = if sign_flip {
            -final_samples[i].drop_m
        } else {
            final_samples[i].drop_m
        };

        // Near zero: trajectory goes from negative/zero to positive (crossing LOS upward)
        if !found_near && prev_drop <= 0.0 && curr_drop > 0.0 && final_samples[i].distance_m > 5.0 {
            // Interpolate
            let frac = (-prev_drop) / (curr_drop - prev_drop);
            near_zero_m = final_samples[i - 1].distance_m
                + frac * (final_samples[i].distance_m - final_samples[i - 1].distance_m);
            found_near = true;
        }

        // Far zero: trajectory goes from positive to negative (crossing LOS downward)
        if found_near && !found_far && prev_drop > 0.0 && curr_drop <= 0.0 {
            let frac = prev_drop / (prev_drop - curr_drop);
            far_zero_m = final_samples[i - 1].distance_m
                + frac * (final_samples[i].distance_m - final_samples[i - 1].distance_m);
            found_far = true;
        }

        // MPBR: where drop goes below -half_vital_m
        if found_far && curr_drop < -half_vital_m && mpbr_m == 0.0 {
            let frac = (-half_vital_m - prev_drop) / (curr_drop - prev_drop);
            mpbr_m = final_samples[i - 1].distance_m
                + frac * (final_samples[i].distance_m - final_samples[i - 1].distance_m);
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
            println!(
                "max_ordinate_distance,{:.1},{}",
                max_ord_dist_display, dist_unit
            );
            println!("impact_velocity,{:.0},{}", impact_vel_display, vel_unit);
            println!("impact_energy,{:.0},{}", impact_energy_display, energy_unit);
        }
        OutputFormat::Table | OutputFormat::Pdf => {
            println!();
            println!(
                "MPBR Analysis (vital zone: {:.1} {})",
                vital_zone_display, size_unit
            );
            println!("╔════════════════════════════════════════╗");
            println!(
                "║  MPBR:            {:>6.0} {:<14}║",
                mpbr_display, dist_unit
            );
            println!(
                "║  Optimal zero:    {:>6.0} {:<14}║",
                optimal_zero_display, dist_unit
            );
            println!(
                "║  Near zero:       {:>6.0} {:<14}║",
                near_zero_display, dist_unit
            );
            println!(
                "║  Far zero:        {:>6.0} {:<14}║",
                far_zero_display, dist_unit
            );
            println!(
                "║  Max ordinate:    {:>6.1} {} at {:.0} {}  ║",
                max_ord_display, size_unit, max_ord_dist_display, dist_unit
            );
            println!(
                "║  Impact velocity: {:>6.0} {:<14}║",
                impact_vel_display, vel_unit
            );
            println!(
                "║  Impact energy:   {:>6.0} {:<14}║",
                impact_energy_display, energy_unit
            );
            println!("╚════════════════════════════════════════╝");
        }
    }

    Ok(())
}

/// Come-ups handler
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments mirror the stable come-ups CLI command shape"
)]
fn handle_come_ups(
    velocity: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    drag_model: DragModelArg,
    zero_distance: f64,
    start: f64,
    end: f64,
    step: f64,
    adjustment_unit: AdjustmentUnit,
    // Resolved turret elevation click graduation for `--adjustment-unit clicks` (MBA-1355):
    // Some(...) iff adjustment_unit == Clicks (resolved once, eagerly, in the ComeUps
    // command's dispatch — see resolve_click_values); None otherwise.
    elevation_click: Option<ClickValue>,
    sight_height: f64,
    temperature: f64,
    pressure: f64,
    humidity: f64,
    altitude: f64,
    wind_speed: f64,
    wind_direction: f64,
    units: UnitSystem,
    output: OutputFormat,
    // MBA-1323 Phase 2: saved-profile velocity-BC segments / Mach-Cd drag curve, already
    // resolved to engine shapes by the caller. See build_trajectory_components's doc comment.
    bc_segments_data: Option<Vec<BCSegmentData>>,
    custom_drag_table: Option<ballistics_engine::drag::DragTable>,
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

    // The zero-angle solve must use the SAME velocity-keyed BC / drag curve as the flight
    // below it (otherwise a segment or curve that changes early-flight drag would mis-zero
    // the shot) — same reasoning as the `trajectory` command's shared resolution.
    let zero_inputs = BallisticInputs {
        bc_value: bc,
        bc_type: drag_model_enum,
        bullet_mass: mass_kg,
        muzzle_velocity: velocity_m,
        bullet_diameter: diameter_m,
        bullet_length: fallback_bullet_length_m(diameter_m, mass_kg),
        sight_height: sight_height_m,
        use_rk4: true,
        use_bc_segments: bc_segments_data.is_some(),
        bc_segments_data: bc_segments_data.clone(),
        custom_drag_table: custom_drag_table.clone(),
        ..Default::default()
    };

    let atmosphere = AtmosphericConditions {
        temperature: temperature_c,
        pressure: pressure_hpa,
        humidity,
        altitude: altitude_m,
    };

    let zero_angle = ballistics_engine::calculate_zero_angle_with_conditions(
        zero_inputs,
        zero_distance_m,
        sight_height_m,
        WindConditions::default(),
        atmosphere,
    )?;

    // Run trajectory with sampling
    let samples = run_sampled_trajectory(
        velocity_m,
        bc,
        mass_kg,
        diameter_m,
        drag_model,
        sight_height_m,
        temperature_c,
        pressure_hpa,
        humidity,
        altitude_m,
        wind_speed_m,
        wind_direction, // degrees, converted internally
        end_m * 1.1,
        sample_m,
        zero_angle,
        bc_segments_data,
        custom_drag_table,
    )?;

    // Build output rows at the requested range intervals
    let adj_label = match adjustment_unit {
        AdjustmentUnit::Mil => "MIL",
        AdjustmentUnit::Moa => "MOA",
        AdjustmentUnit::Smoa => "SMOA",
        AdjustmentUnit::Iphy => "IPHY",
        AdjustmentUnit::Clicks => "CLICKS",
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
            (a.distance_m - range_m)
                .abs()
                .partial_cmp(&(b.distance_m - range_m).abs())
                .unwrap()
        });

        if let Some(sample) = closest {
            if (sample.distance_m - range_m).abs() < sample_m * 1.5 {
                let drop_yd = UnitConverter::distance_from_metric(sample.drop_m, units);
                let range_display = UnitConverter::distance_from_metric(sample.distance_m, units);
                let drop_adj =
                    adjustment_display(drop_yd, range_display, adjustment_unit, elevation_click);
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
            let json_rows: Vec<serde_json::Value> = rows
                .iter()
                .map(|r| {
                    serde_json::json!({
                        "range": r.range,
                        "drop": r.drop_adj,
                        "come_up": r.come_up,
                        "velocity": r.velocity,
                        "energy": r.energy,
                        "time": r.time,
                    })
                })
                .collect();
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
            println!(
                "range_{},drop_{},come_up_{},velocity_{},energy_{},time_s",
                dist_unit,
                adj_label.to_lowercase(),
                adj_label.to_lowercase(),
                vel_unit,
                energy_unit
            );
            // MBA-1355: clicks are whole numbers — drop the 3-decimal formatting used by
            // every angular unit so a clicks CSV doesn't print "10.000".
            let is_clicks = matches!(adjustment_unit, AdjustmentUnit::Clicks);
            for r in &rows {
                if is_clicks {
                    println!(
                        "{:.0},{:.0},{:.0},{:.0},{:.0},{:.3}",
                        r.range, r.drop_adj, r.come_up, r.velocity, r.energy, r.time
                    );
                } else {
                    println!(
                        "{:.0},{:.3},{:.3},{:.0},{:.0},{:.3}",
                        r.range, r.drop_adj, r.come_up, r.velocity, r.energy, r.time
                    );
                }
            }
        }
        OutputFormat::Table | OutputFormat::Pdf => {
            println!();
            println!(
                "Come-Up Table (zero: {:.0} {}, {})",
                zero_distance, dist_unit, adj_label
            );
            // MBA-1355 Task 2 fix pass: the Drop column's border/data width grows with
            // the adjustment label so "Drop (SMOA)"/"Drop (IPHY)" (4 chars) and
            // "Drop (CLICKS)" (6 chars) no longer overflow the original 3-char-label
            // (MIL/MOA) column. MIL/MOA keep the original 10-wide column exactly —
            // `come_up_drop_label_width` floors at 3 — so `default_trajectory_header_is_stable`
            // stays byte-identical.
            let drop_label_w = come_up_drop_label_width(adj_label);
            let drop_dashes = "─".repeat(drop_label_w + 7);
            let ten = "─".repeat(10);
            println!("┌{ten}┬{drop_dashes}┬{ten}┬{ten}┬{ten}┬{ten}┐");
            println!("{}", come_up_header_line(dist_unit, adj_label, vel_unit));
            println!("├{ten}┼{drop_dashes}┼{ten}┼{ten}┼{ten}┼{ten}┤");
            // MBA-1355: clicks are whole numbers — the Drop/Come-Up columns drop the
            // 3-decimal formatting used by every angular unit so the table shows clean
            // integer click counts instead of e.g. "10.000".
            let is_clicks = matches!(adjustment_unit, AdjustmentUnit::Clicks);
            let drop_field_w = drop_label_w + 6;
            for (i, r) in rows.iter().enumerate() {
                let drop_str = if is_clicks {
                    format!("{:>drop_field_w$.0}", r.drop_adj)
                } else {
                    format!("{:>drop_field_w$.3}", r.drop_adj)
                };
                let come_up_str = if i == 0 {
                    "    —     ".to_string()
                } else if is_clicks {
                    format!("{:>9.0} ", r.come_up)
                } else {
                    format!("{:>9.3} ", r.come_up)
                };
                println!(
                    "│{:>9.0} │{} │{}│{:>9.0} │{:>9.0} │{:>9.3} │",
                    r.range, drop_str, come_up_str, r.velocity, r.energy, r.time
                );
            }
            println!("└{ten}┴{drop_dashes}┴{ten}┴{ten}┴{ten}┴{ten}┘");
        }
    }

    Ok(())
}

/// Width of the Drop column's `(LABEL)` slot: the original 3-char assumption (MIL/MOA)
/// floors it, so those two keep the exact pre-existing 10-wide column; SMOA/IPHY (4
/// chars) and CLICKS (6 chars) widen it so the label doesn't overflow the border
/// (MBA-1355 Task 2 fix pass).
fn come_up_drop_label_width(adj_label: &str) -> usize {
    adj_label.len().max(3)
}

/// Come-ups Table header line, extracted so its exact text for the default unit can be
/// pinned by a test without spinning up the whole CLI (MBA-1355 regression guard: a
/// future Clicks-formatting change must not silently reformat the existing MIL/MOA/
/// SMOA/IPHY header).
fn come_up_header_line(dist_unit: &str, adj_label: &str, vel_unit: &str) -> String {
    let w = come_up_drop_label_width(adj_label);
    format!(
        "│Range ({:>2})|Drop ({:>w$})|Come-Up   │ Vel ({:>3})│Energy    │ Time (s) │",
        dist_unit, adj_label, vel_unit, w = w
    )
}

/// Standalone powder-temperature velocity resolution (MBA-737) — no trajectory solve.
/// Delegates to the same `resolve_powder_adjusted_velocity` the trajectory and lead
/// solvers call. NOTE: this command always applies the linear model, while
/// trajectory/lead gate it behind --use-powder-sensitivity (a curve applies there
/// unconditionally).
#[allow(clippy::too_many_arguments)]
fn handle_powder(
    velocity: Option<f64>,
    temperature: Option<f64>,
    powder_temp_sensitivity: Option<f64>,
    powder_temp: Option<f64>,
    powder_temp_curve_si: Option<Vec<(f64, f64)>>,
    sweep: Option<String>,
    mass: Option<f64>,
    units: UnitSystem,
    output: OutputFormat,
) -> Result<(), String> {
    use ballistics_engine::cli_api::parse_powder_sweep;
    use ballistics_engine::resolve_powder_adjusted_velocity;

    let has_curve = powder_temp_curve_si.as_ref().is_some_and(|c| !c.is_empty());
    let curve_ref = powder_temp_curve_si.as_deref();

    // Display-unit defaults identical to the trajectory/lead commands.
    let sens_display = powder_temp_sensitivity.unwrap_or(match units {
        UnitSystem::Imperial => 1.0,
        UnitSystem::Metric => 0.3048 / (5.0 / 9.0),
    });
    let sens_si = UnitConverter::velocity_to_metric(sens_display, units)
        / UnitConverter::temperature_delta_to_metric(1.0, units);
    // Linear model: --powder-temp is the reference temperature the stated velocity
    // was measured at (70F equivalent when omitted). Curve: it is the lookup temp.
    let ref_temp_c = powder_temp
        .map(|t| UnitConverter::temperature_to_metric(t, units))
        .unwrap_or(DEFAULT_POWDER_REFERENCE_TEMP_C);
    let powder_curve_temp_c: Option<f64> =
        powder_temp.map(|t| UnitConverter::temperature_to_metric(t, units));
    let ambient_c: Option<f64> =
        temperature.map(|t| UnitConverter::temperature_to_metric(t, units));
    let velocity_si: Option<f64> = velocity.map(|v| UnitConverter::velocity_to_metric(v, units));
    let mass_kg: Option<f64> = mass.map(|m| UnitConverter::mass_to_metric(m, units));

    // Validation: the linear model needs a nominal velocity; a single (non-sweep)
    // resolution additionally needs a temperature to resolve at.
    if !has_curve && velocity.is_none() {
        return Err(
            "--velocity is required for the linear model (or supply --powder-temp-curve)"
                .to_string(),
        );
    }
    let sweep_temps: Option<Vec<f64>> = match sweep.as_deref() {
        Some(s) => Some(parse_powder_sweep(s)?),
        None => None,
    };
    if sweep_temps.is_none() {
        if !has_curve && temperature.is_none() {
            return Err("--temperature is required (the shot-day air temperature), or use --sweep"
                .to_string());
        }
        if has_curve && powder_temp.is_none() && temperature.is_none() {
            return Err(
                "--powder-temp or --temperature is required with --powder-temp-curve, or use --sweep"
                    .to_string(),
            );
        }
    }

    // Resolve one velocity (m/s) at a display-unit temperature. For the linear model
    // the temperature is the shot-day ambient; for a curve it is the powder
    // temperature the curve is interpolated at.
    let resolve_at_c = |temp_c: f64| -> f64 {
        resolve_powder_adjusted_velocity(
            velocity_si.unwrap_or(0.0),
            temp_c,
            !has_curve, // linear model implied by running this command without a curve
            sens_si,
            ref_temp_c,
            curve_ref,
            None, // per-row lookup IS temp_c; never pin sweep rows to one powder temp
        )
    };

    let energy_at = |v_mps: f64| -> Option<f64> {
        mass_kg.map(|m| UnitConverter::energy_from_metric(0.5 * m * v_mps * v_mps, units))
    };
    let (vel_unit, temp_unit, energy_unit) = match units {
        UnitSystem::Imperial => ("fps", "°F", "ft·lb"),
        UnitSystem::Metric => ("m/s", "°C", "J"),
    };

    struct PowderRow {
        temp_display: f64,
        velocity_display: f64,
        shift_display: Option<f64>,
        energy_display: Option<f64>,
    }
    let row_at = |temp_display: f64, temp_c: f64| -> PowderRow {
        let v_mps = resolve_at_c(temp_c);
        let v_display = UnitConverter::velocity_from_metric(v_mps, units);
        PowderRow {
            temp_display,
            velocity_display: v_display,
            shift_display: velocity.map(|nominal| v_display - nominal),
            energy_display: energy_at(v_mps),
        }
    };

    let rows: Vec<PowderRow> = match &sweep_temps {
        Some(temps) => temps
            .iter()
            .map(|&t| row_at(t, UnitConverter::temperature_to_metric(t, units)))
            .collect(),
        None => {
            // Single resolution. Linear: at the ambient temperature. Curve: at the
            // powder temperature (falling back to ambient), matching the solver.
            let (t_display, t_c) = if has_curve {
                match (powder_temp, powder_curve_temp_c, temperature, ambient_c) {
                    (Some(td), Some(tc), _, _) => (td, tc),
                    (_, _, Some(td), Some(tc)) => (td, tc),
                    _ => unreachable!("validated above"),
                }
            } else {
                (temperature.expect("validated above"), ambient_c.expect("validated above"))
            };
            vec![row_at(t_display, t_c)]
        }
    };

    let model_desc = if has_curve {
        let curve = curve_ref.expect("has_curve");
        let lo = UnitConverter::temperature_from_metric(curve[0].0, units);
        let hi = UnitConverter::temperature_from_metric(curve[curve.len() - 1].0, units);
        format!(
            "measured curve, {} points ({:.1}\u{2013}{:.1} {})",
            curve.len(),
            lo,
            hi,
            temp_unit
        )
    } else {
        format!("linear, {:.2} {}/{}", sens_display, vel_unit, temp_unit)
    };

    match output {
        OutputFormat::Table => {
            println!("Powder Temperature Velocity");
            println!("===========================");
            println!("  Model:              {}", model_desc);
            if !has_curve {
                println!(
                    "  Reference temp:     {:.1} {}",
                    UnitConverter::temperature_from_metric(ref_temp_c, units),
                    temp_unit
                );
            }
            if let Some(v) = velocity {
                println!("  Nominal velocity:   {:.1} {}", v, vel_unit);
            }
            println!();
            if rows.len() == 1 {
                let r = &rows[0];
                let temp_label = if has_curve { "Powder temp:" } else { "Shot temp:" };
                println!("  {:<20}{:.1} {}", temp_label, r.temp_display, temp_unit);
                match r.shift_display {
                    Some(s) => println!(
                        "  Resolved velocity:  {:.1} {}  ({:+.1})",
                        r.velocity_display, vel_unit, s
                    ),
                    None => println!(
                        "  Resolved velocity:  {:.1} {}",
                        r.velocity_display, vel_unit
                    ),
                }
                if let Some(e) = r.energy_display {
                    println!("  Muzzle energy:      {:.0} {}", e, energy_unit);
                }
            } else {
                let has_shift = rows.first().is_some_and(|r| r.shift_display.is_some());
                let has_energy = rows.first().is_some_and(|r| r.energy_display.is_some());
                print!("  {:>10}  {:>14}", format!("Temp ({})", temp_unit), format!("Velocity ({})", vel_unit));
                if has_shift {
                    print!("  {:>12}", format!("Shift ({})", vel_unit));
                }
                if has_energy {
                    print!("  {:>15}", format!("Energy ({})", energy_unit));
                }
                println!();
                for r in &rows {
                    print!("  {:>10.1}  {:>14.1}", r.temp_display, r.velocity_display);
                    if let Some(s) = r.shift_display {
                        print!("  {:>12.1}", s);
                    }
                    if let Some(e) = r.energy_display {
                        print!("  {:>15.0}", e);
                    }
                    println!();
                }
            }
        }
        OutputFormat::Json => {
            let round1 = |x: f64| (x * 10.0).round() / 10.0;
            let rows_json: Vec<serde_json::Value> = rows
                .iter()
                .map(|r| {
                    let mut o = serde_json::json!({
                        "temperature": round1(r.temp_display),
                        "velocity": round1(r.velocity_display),
                    });
                    if let Some(s) = r.shift_display {
                        o["shift"] = serde_json::json!(round1(s));
                    }
                    if let Some(e) = r.energy_display {
                        o["energy"] = serde_json::json!(e.round());
                    }
                    o
                })
                .collect();
            let mut result = serde_json::json!({
                "command": "powder",
                "units": match units {
                    UnitSystem::Imperial => "imperial",
                    UnitSystem::Metric => "metric",
                },
                "model": if has_curve { "curve" } else { "linear" },
            });
            if !has_curve {
                result["sensitivity"] = serde_json::json!(sens_display);
                result["reference_temp"] =
                    serde_json::json!(round1(UnitConverter::temperature_from_metric(ref_temp_c, units)));
            } else {
                result["curve_points"] = serde_json::json!(curve_ref.expect("has_curve").len());
            }
            if let Some(v) = velocity {
                result["nominal_velocity"] = serde_json::json!(v);
            }
            if let Some(m) = mass {
                result["mass"] = serde_json::json!(m);
            }
            if sweep_temps.is_some() {
                result["sweep"] = serde_json::json!(rows_json);
            } else {
                result["resolved"] = rows_json.into_iter().next().unwrap_or(serde_json::json!(null));
            }
            println!("{}", serde_json::to_string_pretty(&result).map_err(|e| e.to_string())?);
        }
        OutputFormat::Csv => {
            let has_shift = rows.first().is_some_and(|r| r.shift_display.is_some());
            let has_energy = rows.first().is_some_and(|r| r.energy_display.is_some());
            let mut header = format!("temperature_{},velocity_{}", temp_unit_ascii(units), vel_unit_ascii(units));
            if has_shift {
                header.push_str(&format!(",shift_{}", vel_unit_ascii(units)));
            }
            if has_energy {
                header.push_str(&format!(",energy_{}", energy_unit_ascii(units)));
            }
            println!("{}", header);
            for r in &rows {
                let mut line = format!("{:.1},{:.1}", r.temp_display, r.velocity_display);
                if let Some(s) = r.shift_display {
                    line.push_str(&format!(",{:.1}", s));
                }
                if let Some(e) = r.energy_display {
                    line.push_str(&format!(",{:.0}", e));
                }
                println!("{}", line);
            }
        }
        OutputFormat::Pdf => {
            return Err("PDF output is not supported for the powder command".to_string());
        }
    }
    Ok(())
}

/// ASCII (CSV-header-safe) unit suffixes for the powder command's CSV output.
fn temp_unit_ascii(units: UnitSystem) -> &'static str {
    match units {
        UnitSystem::Imperial => "f",
        UnitSystem::Metric => "c",
    }
}
fn vel_unit_ascii(units: UnitSystem) -> &'static str {
    match units {
        UnitSystem::Imperial => "fps",
        UnitSystem::Metric => "mps",
    }
}
fn energy_unit_ascii(units: UnitSystem) -> &'static str {
    match units {
        UnitSystem::Imperial => "ftlb",
        UnitSystem::Metric => "j",
    }
}

/// Moving-target lead table handler (MBA-1287).
///
/// Runs `ballistics_engine::calculate_lead` once per range in the sweep. A single
/// out-of-domain range (e.g. an inbound target that overtakes the shooter) must not
/// abort the whole table: the error is captured per-row and rendered as a dashed row
/// (table) or an `"error"` field (JSON), and the sweep continues.
#[allow(clippy::too_many_arguments)]
fn handle_lead(
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
    use_powder_sensitivity: bool,
    powder_temp_sensitivity: f64,
    powder_temp: f64,
    powder_temp_curve: Option<Vec<(f64, f64)>>,
    powder_curve_temp_c: Option<f64>,
    target_speed: f64,
    target_angle: f64,
    target_length: Option<f64>,
    start: f64,
    end: f64,
    step: f64,
    adjustment_unit: AdjustmentUnit,
    units: UnitSystem,
    output: OutputFormat,
    // MBA-1323 Phase 2: saved-profile velocity-BC segments / Mach-Cd drag curve, already
    // resolved to engine shapes by the caller. See build_trajectory_components's doc comment.
    bc_segments_data: Option<Vec<BCSegmentData>>,
    custom_drag_table: Option<ballistics_engine::drag::DragTable>,
) -> Result<(), Box<dyn Error>> {
    // MBA-1355: `lead` (moving-target) is out of scope for `--adjustment-unit clicks` — reject
    // before doing any work rather than silently falling back to MIL below.
    reject_clicks_out_of_scope(adjustment_unit);
    // Convert to metric
    let velocity_m = UnitConverter::velocity_to_metric(velocity, units);
    let mass_kg = UnitConverter::mass_to_metric(mass, units);
    let diameter_m = UnitConverter::diameter_to_metric(diameter, units);
    let sight_height_m = UnitConverter::sight_height_to_metric(sight_height, units);
    let temperature_c = UnitConverter::temperature_to_metric(temperature, units);
    let pressure_hpa = UnitConverter::pressure_to_metric(pressure, units);
    let altitude_m = UnitConverter::altitude_to_metric(altitude, units);
    let wind_speed_m = UnitConverter::wind_to_metric(wind_speed, units);
    let target_speed_mps = UnitConverter::wind_to_metric(target_speed, units);
    let target_length_m = target_length.map(|l| match units {
        UnitSystem::Imperial => l * 0.0254, // inches to meters
        UnitSystem::Metric => l / 1000.0,   // mm to meters
    });
    let end_m = UnitConverter::distance_to_metric(end, units);
    let step_m = UnitConverter::distance_to_metric(step, units);

    let (mut inputs, wind, atmosphere) = build_trajectory_components(
        velocity_m,
        bc,
        mass_kg,
        diameter_m,
        drag_model,
        sight_height_m,
        temperature_c,
        pressure_hpa,
        humidity,
        altitude_m,
        wind_speed_m,
        wind_direction,
        end_m,
        step_m,
        bc_segments_data,
        custom_drag_table,
    );

    // Powder temperature (MBA-1325): identical resolution to `run_trajectory` so a
    // powder-temp-curve/-sensitivity muzzle-velocity correction reproduces exactly —
    // TrajectorySolver::new (invoked inside calculate_lead) resolves the correction from
    // these BallisticInputs fields, so setting them here is the entire plumbing needed.
    inputs.use_powder_sensitivity = use_powder_sensitivity;
    inputs.powder_temp_sensitivity = if use_powder_sensitivity {
        UnitConverter::velocity_to_metric(powder_temp_sensitivity, units)
            / UnitConverter::temperature_delta_to_metric(1.0, units)
    } else {
        0.0
    };
    inputs.powder_temp = UnitConverter::temperature_to_metric(powder_temp, units);
    inputs.powder_temp_curve = powder_temp_curve;
    inputs.powder_curve_temp_c = powder_curve_temp_c;

    let adj_label = match adjustment_unit {
        AdjustmentUnit::Mil => "MIL",
        AdjustmentUnit::Moa => "MOA",
        AdjustmentUnit::Smoa => "SMOA",
        AdjustmentUnit::Iphy => "IPHY",
        AdjustmentUnit::Clicks => {
            debug_assert!(false, "Clicks must be resolved via clicks_for() before display");
            "MIL"
        }
    };
    let (dist_unit, speed_unit) = match units {
        UnitSystem::Imperial => ("yd", "mph"),
        UnitSystem::Metric => ("m", "m/s"),
    };
    let has_bodies = target_length_m.is_some();

    struct LeadRow {
        range: f64,
        result: Result<ballistics_engine::LeadSolution, ballistics_engine::LeadError>,
        bodies: Option<f64>,
    }

    let mut rows: Vec<LeadRow> = Vec::new();
    let mut current_range = start;
    while current_range <= end + 0.1 {
        let range_m = UnitConverter::distance_to_metric(current_range, units);

        let result = ballistics_engine::calculate_lead(
            inputs.clone(),
            wind.clone(),
            atmosphere.clone(),
            target_speed_mps,
            target_angle,
            range_m,
        );

        let bodies = match (&result, target_length_m) {
            (Ok(sol), Some(len_m)) if len_m > 0.0 => Some(sol.lead_m / len_m),
            _ => None,
        };

        rows.push(LeadRow {
            range: current_range,
            result,
            bodies,
        });

        current_range += step;
    }

    match output {
        OutputFormat::Json => {
            let json_rows: Vec<serde_json::Value> = rows
                .iter()
                .map(|r| match &r.result {
                    Ok(sol) => {
                        let mut row = serde_json::json!({
                            "range": r.range,
                            "tof_s": sol.time_of_flight_s,
                            "lead": UnitConverter::distance_from_metric(sol.lead_m, units),
                            "lead_mil": sol.lead_mil,
                            "lead_moa": sol.lead_moa,
                            "intercept_range": UnitConverter::distance_from_metric(
                                sol.corrected_range_m,
                                units
                            ),
                            "iterations": sol.iterations,
                        });
                        if let Some(bodies) = r.bodies {
                            row["bodies"] = serde_json::json!(bodies);
                        }
                        row
                    }
                    Err(e) => serde_json::json!({
                        "range": r.range,
                        "error": e.to_string(),
                    }),
                })
                .collect();

            let json = serde_json::json!({
                "target_speed": target_speed,
                "target_speed_unit": speed_unit,
                "target_angle": target_angle,
                "units": match units {
                    UnitSystem::Imperial => "imperial",
                    UnitSystem::Metric => "metric",
                },
                "distance_unit": dist_unit,
                "adjustment_unit": adj_label,
                "rows": json_rows,
            });
            println!("{}", serde_json::to_string_pretty(&json)?);
        }
        OutputFormat::Csv => {
            if has_bodies {
                println!(
                    "range_{dist_unit},tof_s,lead_{dist_unit},lead_{unit},intercept_{dist_unit},iterations,bodies,error",
                    dist_unit = dist_unit,
                    unit = adj_label.to_lowercase()
                );
            } else {
                println!(
                    "range_{dist_unit},tof_s,lead_{dist_unit},lead_{unit},intercept_{dist_unit},iterations,error",
                    dist_unit = dist_unit,
                    unit = adj_label.to_lowercase()
                );
            }
            for r in &rows {
                match &r.result {
                    Ok(sol) => {
                        let lead_disp = UnitConverter::distance_from_metric(sol.lead_m, units);
                        let lead_adj = match adjustment_unit {
                            AdjustmentUnit::Mil => sol.lead_mil,
                            AdjustmentUnit::Moa => sol.lead_moa,
                            AdjustmentUnit::Smoa | AdjustmentUnit::Iphy => sol.lead_mil * smoa_per_mil(),
                            AdjustmentUnit::Clicks => {
                                debug_assert!(false, "Clicks must be resolved via clicks_for() before display");
                                sol.lead_mil
                            }
                        };
                        let intercept_disp =
                            UnitConverter::distance_from_metric(sol.corrected_range_m, units);
                        if has_bodies {
                            println!(
                                "{:.0},{:.3},{:.2},{:.3},{:.1},{},{},",
                                r.range,
                                sol.time_of_flight_s,
                                lead_disp,
                                lead_adj,
                                intercept_disp,
                                sol.iterations,
                                r.bodies.map(|b| format!("{:.2}", b)).unwrap_or_default()
                            );
                        } else {
                            println!(
                                "{:.0},{:.3},{:.2},{:.3},{:.1},{},",
                                r.range,
                                sol.time_of_flight_s,
                                lead_disp,
                                lead_adj,
                                intercept_disp,
                                sol.iterations
                            );
                        }
                    }
                    Err(e) => {
                        if has_bodies {
                            println!("{:.0},,,,,,,{}", r.range, e);
                        } else {
                            println!("{:.0},,,,,,{}", r.range, e);
                        }
                    }
                }
            }
        }
        OutputFormat::Table | OutputFormat::Pdf => {
            println!();
            println!(
                "Moving-Target Lead Table (target speed: {:.1} {}, angle: {:.0}\u{b0}, {})",
                target_speed, speed_unit, target_angle, adj_label
            );
            println!("Positive lead = hold in the direction of target travel.");
            println!();

            if has_bodies {
                println!("┌──────────┬──────────┬──────────┬──────────┬──────────┬──────────┐");
                println!(
                    "│Range ({:>2})│TOF (s)   │Lead ({:>3})│Lead ({:>3})│Intercept │Bodies    │",
                    dist_unit, dist_unit, adj_label
                );
                println!("├──────────┼──────────┼──────────┼──────────┼──────────┼──────────┤");
                for r in &rows {
                    match &r.result {
                        Ok(sol) => {
                            let lead_disp = UnitConverter::distance_from_metric(sol.lead_m, units);
                            let lead_adj = match adjustment_unit {
                                AdjustmentUnit::Mil => sol.lead_mil,
                                AdjustmentUnit::Moa => sol.lead_moa,
                                AdjustmentUnit::Smoa | AdjustmentUnit::Iphy => {
                                    sol.lead_mil * smoa_per_mil()
                                }
                                AdjustmentUnit::Clicks => {
                                    debug_assert!(
                                        false,
                                        "Clicks must be resolved via clicks_for() before display"
                                    );
                                    sol.lead_mil
                                }
                            };
                            let intercept_disp =
                                UnitConverter::distance_from_metric(sol.corrected_range_m, units);
                            println!(
                                "│{:>9.0} │{:>9.3} │{:>9.2} │{:>9.3} │{:>9.1} │{:>9.2} │",
                                r.range,
                                sol.time_of_flight_s,
                                lead_disp,
                                lead_adj,
                                intercept_disp,
                                r.bodies.unwrap_or(0.0)
                            );
                        }
                        Err(e) => {
                            println!(
                                "│{:>9.0} │{:>9} │{:>9} │{:>9} │{:>9} │{:>9} │  {}",
                                r.range, "--", "--", "--", "--", "--", e
                            );
                        }
                    }
                }
                println!("└──────────┴──────────┴──────────┴──────────┴──────────┴──────────┘");
            } else {
                println!("┌──────────┬──────────┬──────────┬──────────┬──────────┐");
                println!(
                    "│Range ({:>2})│TOF (s)   │Lead ({:>3})│Lead ({:>3})│Intercept │",
                    dist_unit, dist_unit, adj_label
                );
                println!("├──────────┼──────────┼──────────┼──────────┼──────────┤");
                for r in &rows {
                    match &r.result {
                        Ok(sol) => {
                            let lead_disp = UnitConverter::distance_from_metric(sol.lead_m, units);
                            let lead_adj = match adjustment_unit {
                                AdjustmentUnit::Mil => sol.lead_mil,
                                AdjustmentUnit::Moa => sol.lead_moa,
                                AdjustmentUnit::Smoa | AdjustmentUnit::Iphy => {
                                    sol.lead_mil * smoa_per_mil()
                                }
                                AdjustmentUnit::Clicks => {
                                    debug_assert!(
                                        false,
                                        "Clicks must be resolved via clicks_for() before display"
                                    );
                                    sol.lead_mil
                                }
                            };
                            let intercept_disp =
                                UnitConverter::distance_from_metric(sol.corrected_range_m, units);
                            println!(
                                "│{:>9.0} │{:>9.3} │{:>9.2} │{:>9.3} │{:>9.1} │",
                                r.range, sol.time_of_flight_s, lead_disp, lead_adj, intercept_disp
                            );
                        }
                        Err(e) => {
                            println!(
                                "│{:>9.0} │{:>9} │{:>9} │{:>9} │{:>9} │  {}",
                                r.range, "--", "--", "--", "--", e
                            );
                        }
                    }
                }
                println!("└──────────┴──────────┴──────────┴──────────┴──────────┘");
            }
        }
    }

    Ok(())
}

/// Wind card handler
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments mirror the stable wind-card CLI command shape"
)]
fn handle_wind_card(
    velocity: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    drag_model: DragModelArg,
    zero_distance: f64,
    wind_speeds: &[f64],
    wind_angles: &[f64],
    legacy_labels: bool,
    start: f64,
    end: f64,
    step: f64,
    adjustment_unit: AdjustmentUnit,
    sight_height: f64,
    temperature: f64,
    pressure: f64,
    humidity: f64,
    altitude: f64,
    units: UnitSystem,
    output: OutputFormat,
) -> Result<(), Box<dyn Error>> {
    // MBA-1355: wind card is out of scope for `--adjustment-unit clicks` — reject before doing
    // any work rather than silently falling back to MIL below.
    reject_clicks_out_of_scope(adjustment_unit);
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
        bullet_length: fallback_bullet_length_m(diameter_m, mass_kg),
        sight_height: sight_height_m,
        use_rk4: true,
        ..Default::default()
    };

    let atmosphere = AtmosphericConditions {
        temperature: temperature_c,
        pressure: pressure_hpa,
        humidity,
        altitude: altitude_m,
    };

    let zero_angle = ballistics_engine::calculate_zero_angle_with_conditions(
        zero_inputs,
        zero_distance_m,
        sight_height_m,
        WindConditions::default(),
        atmosphere,
    )?;

    let adj_label = match adjustment_unit {
        AdjustmentUnit::Mil => "MIL",
        AdjustmentUnit::Moa => "MOA",
        AdjustmentUnit::Smoa => "SMOA",
        AdjustmentUnit::Iphy => "IPHY",
        AdjustmentUnit::Clicks => {
            debug_assert!(false, "Clicks must be resolved via clicks_for() before display");
            "MIL"
        }
    };

    let (dist_unit, wind_unit) = match units {
        UnitSystem::Imperial => ("yd", "mph"),
        UnitSystem::Metric => ("m", "m/s"),
    };

    // For each wind angle (wind-FROM degrees, matches --wind-direction), run
    // trajectory at each wind speed and collect drift. Legacy default is a
    // single 90° (full-value crosswind from the right) angle.
    //
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

    let mut json_cards: Vec<serde_json::Value> = Vec::new();

    for (angle_idx, &angle_deg) in wind_angles.iter().enumerate() {
        let mut all_drifts: Vec<Vec<f64>> = vec![Vec::new(); ranges.len()];

        for &ws in wind_speeds {
            let ws_m = UnitConverter::wind_to_metric(ws, units);

            let samples = run_sampled_trajectory(
                velocity_m,
                bc,
                mass_kg,
                diameter_m,
                drag_model,
                sight_height_m,
                temperature_c,
                pressure_hpa,
                humidity,
                altitude_m,
                ws_m,
                angle_deg,
                end_m * 1.1,
                sample_m,
                zero_angle,
                // wind-card does not yet consume saved-profile bc_segments/drag_curve
                // (MBA-1323 Phase 2 follow-up) — see CLI_USAGE.md's a7p import section.
                None,
                None,
            )?;

            for (ri, &range_display) in ranges.iter().enumerate() {
                let range_m = UnitConverter::distance_to_metric(range_display, units);

                let closest = samples.iter().min_by(|a, b| {
                    (a.distance_m - range_m)
                        .abs()
                        .partial_cmp(&(b.distance_m - range_m).abs())
                        .unwrap()
                });

                let drift_adj = if let Some(sample) = closest {
                    if (sample.distance_m - range_m).abs() < sample_m * 1.5 {
                        let drift_yd =
                            UnitConverter::distance_from_metric(sample.wind_drift_m, units);
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

        let wind_rows: Vec<WindRow> = ranges
            .iter()
            .enumerate()
            .map(|(i, &range)| WindRow {
                range,
                drifts: all_drifts[i].clone(),
            })
            .collect();

        match output {
            OutputFormat::Json => {
                let json_rows: Vec<serde_json::Value> = wind_rows
                    .iter()
                    .map(|r| {
                        let mut row = serde_json::json!({ "range": r.range });
                        for (j, &ws) in wind_speeds.iter().enumerate() {
                            row[format!("wind_{}", ws)] =
                                serde_json::json!(r.drifts.get(j).unwrap_or(&0.0));
                        }
                        row
                    })
                    .collect();
                let mut card = serde_json::json!({
                    "zero_distance": zero_distance,
                    "adjustment_unit": adj_label,
                    "distance_unit": dist_unit,
                    "wind_unit": wind_unit,
                    "wind_speeds": wind_speeds,
                    "data": json_rows,
                });
                if legacy_labels {
                    card["crosswind"] = serde_json::json!("full-value (90°)");
                } else {
                    card["wind_angle"] = serde_json::json!(angle_deg);
                }
                json_cards.push(card);
            }
            OutputFormat::Csv => {
                if !legacy_labels {
                    if angle_idx > 0 {
                        println!();
                    }
                    println!("# wind_angle={}", angle_deg);
                }
                let ws_headers: Vec<String> = wind_speeds
                    .iter()
                    .map(|ws| format!("wind_{}_{}", ws, wind_unit))
                    .collect();
                println!("range_{},{}", dist_unit, ws_headers.join(","));
                for r in &wind_rows {
                    let drift_strs: Vec<String> = r
                        .drifts
                        .iter()
                        .map(|d| {
                            let d = if d.abs() < 1e-9 { 0.0 } else { *d };
                            format!("{:.1}", d)
                        })
                        .collect();
                    println!("{:.0},{}", r.range, drift_strs.join(","));
                }
            }
            OutputFormat::Table | OutputFormat::Pdf => {
                println!();
                if legacy_labels {
                    println!(
                        "Wind Card (zero: {:.0} {}, {}, full-value crosswind)",
                        zero_distance, dist_unit, adj_label
                    );
                } else {
                    println!(
                        "Wind Card (zero: {:.0} {}, {}) — wind angle {}° (wind-FROM: 0=head, 90=right, 180=tail, 270=left)",
                        zero_distance, dist_unit, adj_label, angle_deg
                    );
                }

                // Header
                let col_width = 10;
                let range_header = format!("Range ({:>2})", dist_unit);
                let mut header = format!("┌{:─>w$}", "", w = col_width);
                for _ in wind_speeds {
                    header += &format!("┬{:─>w$}", "", w = col_width);
                }
                header += "┐";
                println!("{}", header);

                let mut label_row = format!("│{:<w$}", range_header, w = col_width);
                for ws in wind_speeds {
                    label_row += &format!("│{:>8} {} ", ws, wind_unit);
                }
                label_row += "│";
                println!("{}", label_row);

                let mut sep = format!("├{:─>w$}", "", w = col_width);
                for _ in wind_speeds {
                    sep += &format!("┼{:─>w$}", "", w = col_width);
                }
                sep += "┤";
                println!("{}", sep);

                for r in &wind_rows {
                    let mut row_str = format!("│{:>9.0} ", r.range);
                    for d in &r.drifts {
                        let d = if d.abs() < 1e-9 { 0.0 } else { *d };
                        row_str += &format!("│{:>9.1} ", d);
                    }
                    row_str += "│";
                    println!("{}", row_str);
                }

                let mut footer = format!("└{:─>w$}", "", w = col_width);
                for _ in wind_speeds {
                    footer += &format!("┴{:─>w$}", "", w = col_width);
                }
                footer += "┘";
                println!("{}", footer);
            }
        }
    }

    if let OutputFormat::Json = output {
        if json_cards.len() == 1 {
            println!("{}", serde_json::to_string_pretty(&json_cards[0])?);
        } else {
            let array = serde_json::Value::Array(json_cards);
            println!("{}", serde_json::to_string_pretty(&array)?);
        }
    }

    Ok(())
}

// ============================================================================
// Stability Analysis Handler (MBA-734)
// ============================================================================

/// Calculate the Miller stability factor and minimum twist rates
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments mirror the stable stability CLI command shape"
)]
fn handle_stability(
    mass: f64,
    diameter: f64,
    length: f64,
    twist_rate: f64,
    velocity: f64,
    temperature: f64,
    pressure: f64,
    altitude: f64,
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
        let mut low = 1.0_f64; // 1 inch/turn (extremely fast)
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
                &test_inputs,
                atmo_params,
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
            println!(
                "{:.2},{},{:.1},{:.1},{:.1},{:.3},{:.0}",
                sg, status, twist_rate, min15_out, min10_out, length, velocity
            );
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
            println!(
                "║ Muzzle Velocity:     {:>10.0} {:<3}  ║",
                velocity, vel_unit
            );
            println!("╚════════════════════════════════════════╝");
        }
    }

    Ok(())
}

// ============================================================================
// Range Table Handler (MBA-733)
// ============================================================================

/// Generate a comprehensive range table combining drop, wind, velocity, energy, and ToF
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments mirror the stable range-table CLI command shape"
)]
fn handle_range_table(
    velocity: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    drag_model: DragModelArg,
    zero_distance: f64,
    start: f64,
    end: f64,
    step: f64,
    wind_speed: f64,
    wind_direction: f64,
    adjustment_unit: AdjustmentUnit,
    sight_height: f64,
    temperature: f64,
    pressure: f64,
    humidity: f64,
    altitude: f64,
    units: UnitSystem,
    output: OutputFormat,
) -> Result<(), Box<dyn Error>> {
    // MBA-1355: range-table is out of scope for `--adjustment-unit clicks` — reject before doing
    // any work rather than silently falling back to MIL below.
    reject_clicks_out_of_scope(adjustment_unit);
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
        bullet_length: fallback_bullet_length_m(diameter_m, mass_kg),
        sight_height: sight_height_m,
        use_rk4: true,
        ..Default::default()
    };

    let atmosphere = AtmosphericConditions {
        temperature: temperature_c,
        pressure: pressure_hpa,
        humidity,
        altitude: altitude_m,
    };

    let zero_angle = ballistics_engine::calculate_zero_angle_with_conditions(
        zero_inputs,
        zero_distance_m,
        sight_height_m,
        WindConditions::default(),
        atmosphere,
    )?;

    // Run trajectory WITH wind (for wind drift values)
    // range-table does not yet consume saved-profile bc_segments/drag_curve (MBA-1323
    // Phase 2 follow-up) — see CLI_USAGE.md's a7p import section.
    let wind_samples = run_sampled_trajectory(
        velocity_m,
        bc,
        mass_kg,
        diameter_m,
        drag_model,
        sight_height_m,
        temperature_c,
        pressure_hpa,
        humidity,
        altitude_m,
        wind_speed_m,
        wind_direction,
        end_m * 1.1,
        sample_m,
        zero_angle,
        None,
        None,
    )?;

    // Run trajectory WITHOUT wind (for pure drop)
    let no_wind_samples = run_sampled_trajectory(
        velocity_m,
        bc,
        mass_kg,
        diameter_m,
        drag_model,
        sight_height_m,
        temperature_c,
        pressure_hpa,
        humidity,
        altitude_m,
        0.0,
        0.0,
        end_m * 1.1,
        sample_m,
        zero_angle,
        None,
        None,
    )?;

    // Build output rows
    let adj_label = match adjustment_unit {
        AdjustmentUnit::Mil => "MIL",
        AdjustmentUnit::Moa => "MOA",
        AdjustmentUnit::Smoa => "SMOA",
        AdjustmentUnit::Iphy => "IPHY",
        AdjustmentUnit::Clicks => {
            debug_assert!(false, "Clicks must be resolved via clicks_for() before display");
            "MIL"
        }
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
            (a.distance_m - range_m)
                .abs()
                .partial_cmp(&(b.distance_m - range_m).abs())
                .unwrap()
        });

        let w_closest = wind_samples.iter().min_by(|a, b| {
            (a.distance_m - range_m)
                .abs()
                .partial_cmp(&(b.distance_m - range_m).abs())
                .unwrap()
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
            let json_rows: Vec<serde_json::Value> = rows
                .iter()
                .map(|r| {
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
                })
                .collect();
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
            println!(
                "range_{},drop_{},drop_{},wind_{},wind_{},velocity_{},energy_{},tof_s",
                dist_unit,
                drop_unit,
                adj_label.to_lowercase(),
                drop_unit,
                adj_label.to_lowercase(),
                vel_unit,
                energy_unit
            );
            for r in &rows {
                println!(
                    "{:.0},{:.1},{:.3},{:.1},{:.2},{:.0},{:.0},{:.2}",
                    r.range,
                    r.drop_linear,
                    r.drop_adj,
                    r.wind_linear,
                    r.wind_adj,
                    r.velocity,
                    r.energy,
                    r.time
                );
            }
        }
        OutputFormat::Table | OutputFormat::Pdf => {
            println!();
            println!(
                "Range Table (zero: {:.0} {}, {:.0} {} crosswind)",
                zero_distance, dist_unit, wind_speed, wind_unit_label
            );
            println!(
                "┌───────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬───────┐"
            );
            println!(
                "│Range  │Drop     │Drop     │Wind     │Wind     │Vel      │Energy   │ToF    │"
            );
            println!("│({:>2})   │({:>2})     │({:>3})    │({:>2})     │({:>3})    │({:>3})    │({:>4})   │(s)    │",
                     dist_unit, drop_unit, adj_label, drop_unit, adj_label, vel_unit, energy_unit);
            println!(
                "├───────┼─────────┼─────────┼─────────┼─────────┼─────────┼─────────┼───────┤"
            );
            for r in &rows {
                println!(
                    "│{:>6.0} │{:>8.1} │{:>8.3} │{:>8.1} │{:>8.2} │{:>8.0} │{:>8.0} │{:>6.2} │",
                    r.range,
                    r.drop_linear,
                    r.drop_adj,
                    r.wind_linear,
                    r.wind_adj,
                    r.velocity,
                    r.energy,
                    r.time
                );
            }
            println!(
                "└───────┴─────────┴─────────┴─────────┴─────────┴─────────┴─────────┴───────┘"
            );
        }
    }

    Ok(())
}

/// One load in a `compare` run, in DISPLAY units per the session `--units`.
/// `bc_segments_data`/`custom_drag_table` come only from saved profiles (inline
/// `--load` specs are scalar-BC by design) and are threaded into BOTH the zero
/// solve and the sampled runs so the zero uses the same drag physics (MBA-735).
struct CompareLoad {
    name: String,
    drag_model: DragModelArg,
    bc: f64,
    mass: f64,
    velocity: f64,
    diameter: f64,
    bc_segments_data: Option<Vec<BCSegmentData>>,
    custom_drag_table: Option<ballistics_engine::drag::DragTable>,
}

/// Parse a `--load` spec: `NAME:DRAG:BC:MASS:VELOCITY[:DIAMETER]` (MBA-735).
/// Follows the `--wind-segment` colon convention; all numbers are display units.
fn parse_compare_load_spec(spec: &str, units: UnitSystem) -> Result<CompareLoad, Box<dyn Error>> {
    let parts: Vec<&str> = spec.split(':').collect();
    if parts.len() != 5 && parts.len() != 6 {
        return Err(format!(
            "invalid --load '{spec}': expected NAME:DRAG:BC:MASS:VELOCITY[:DIAMETER] \
             (5 or 6 colon-separated fields; NAME must not contain ':')"
        )
        .into());
    }
    let name = parts[0].trim();
    if name.is_empty() {
        return Err(format!("invalid --load '{spec}': NAME is empty").into());
    }
    let drag_model = match parts[1].trim().to_ascii_lowercase().as_str() {
        "g1" => DragModelArg::G1,
        "g7" => DragModelArg::G7,
        other => {
            return Err(format!("invalid --load '{spec}': DRAG '{other}' must be g1 or g7").into())
        }
    };
    let num = |i: usize, field: &str| -> Result<f64, Box<dyn Error>> {
        let v: f64 = parts[i].trim().parse().map_err(|_| {
            format!("invalid --load '{spec}': {field} '{}' is not a number", parts[i])
        })?;
        if !v.is_finite() || v <= 0.0 {
            return Err(format!("invalid --load '{spec}': {field} must be finite and > 0").into());
        }
        Ok(v)
    };
    let bc = num(2, "BC")?;
    let mass = num(3, "MASS")?;
    let velocity = num(4, "VELOCITY")?;
    let diameter = if parts.len() == 6 {
        num(5, "DIAMETER")?
    } else {
        match units {
            UnitSystem::Imperial => 0.308,
            UnitSystem::Metric => 7.82,
        }
    };
    Ok(CompareLoad {
        name: name.to_string(),
        drag_model,
        bc,
        mass,
        velocity,
        diameter,
        bc_segments_data: None,
        custom_drag_table: None,
    })
}

/// Side-by-side multi-load comparison (MBA-735): every load is zeroed independently at
/// the shared zero distance, then run through identical wind/atmosphere; rows mirror
/// `range-table` semantics (pure drop from a no-wind pass, drift from the wind pass).
/// Like range-table, saved-profile bc_segments/drag_curve are not consumed here.
#[allow(clippy::too_many_arguments)] // mirrors the sibling range-table handler
fn handle_compare(
    loads: Vec<CompareLoad>,
    zero_distance: f64,
    start: f64,
    end: f64,
    step: f64,
    wind_speed: f64,
    wind_direction: f64,
    adjustment_unit: AdjustmentUnit,
    sight_height: f64,
    temperature: f64,
    pressure: f64,
    humidity: f64,
    altitude: f64,
    units: UnitSystem,
    output: OutputFormat,
) -> Result<(), Box<dyn Error>> {
    // MBA-1355: compare is out of scope for `--adjustment-unit clicks` — reject before doing any
    // work rather than silently falling back to MIL below.
    reject_clicks_out_of_scope(adjustment_unit);
    // Shared conditions in metric
    let sight_height_m = UnitConverter::sight_height_to_metric(sight_height, units);
    let zero_distance_m = UnitConverter::distance_to_metric(zero_distance, units);
    let temperature_c = UnitConverter::temperature_to_metric(temperature, units);
    let pressure_hpa = UnitConverter::pressure_to_metric(pressure, units);
    let altitude_m = UnitConverter::altitude_to_metric(altitude, units);
    let wind_speed_m = UnitConverter::wind_to_metric(wind_speed, units);
    let end_m = UnitConverter::distance_to_metric(end, units);
    let sample_m = UnitConverter::distance_to_metric(step, units);

    let atmosphere = AtmosphericConditions {
        temperature: temperature_c,
        pressure: pressure_hpa,
        humidity,
        altitude: altitude_m,
    };

    struct LoadRow {
        drop_linear: f64,
        drop_adj: f64,
        wind_linear: f64,
        wind_adj: f64,
        velocity: f64,
        energy: f64,
        time: f64,
    }
    struct LoadResult {
        zero_angle_deg: f64,
        rows: Vec<LoadRow>,
    }

    // The shared range grid (display units)
    let mut ranges: Vec<f64> = Vec::new();
    let mut r = start;
    while r <= end + 0.1 {
        ranges.push(r);
        r += step;
    }

    let mut results: Vec<LoadResult> = Vec::new();
    for load in &loads {
        let velocity_m = UnitConverter::velocity_to_metric(load.velocity, units);
        let mass_kg = UnitConverter::mass_to_metric(load.mass, units);
        let diameter_m = UnitConverter::diameter_to_metric(load.diameter, units);
        let drag_model_enum = match load.drag_model {
            DragModelArg::G1 => DragModel::G1,
            DragModelArg::G7 => DragModel::G7,
        };

        let zero_inputs = BallisticInputs {
            bc_value: load.bc,
            bc_type: drag_model_enum,
            bullet_mass: mass_kg,
            muzzle_velocity: velocity_m,
            bullet_diameter: diameter_m,
            bullet_length: fallback_bullet_length_m(diameter_m, mass_kg),
            sight_height: sight_height_m,
            use_rk4: true,
            // The zero must use the same drag physics as the trajectory runs below.
            use_bc_segments: load.bc_segments_data.is_some(),
            bc_segments_data: load.bc_segments_data.clone(),
            custom_drag_table: load.custom_drag_table.clone(),
            ..Default::default()
        };
        let zero_angle = ballistics_engine::calculate_zero_angle_with_conditions(
            zero_inputs,
            zero_distance_m,
            sight_height_m,
            WindConditions::default(),
            atmosphere.clone(),
        )
        .map_err(|e| format!("load '{}': zeroing failed: {e}", load.name))?;

        let wind_samples = run_sampled_trajectory(
            velocity_m,
            load.bc,
            mass_kg,
            diameter_m,
            load.drag_model,
            sight_height_m,
            temperature_c,
            pressure_hpa,
            humidity,
            altitude_m,
            wind_speed_m,
            wind_direction,
            end_m * 1.1,
            sample_m,
            zero_angle,
            load.bc_segments_data.clone(),
            load.custom_drag_table.clone(),
        )
        .map_err(|e| format!("load '{}': {e}", load.name))?;
        let no_wind_samples = run_sampled_trajectory(
            velocity_m,
            load.bc,
            mass_kg,
            diameter_m,
            load.drag_model,
            sight_height_m,
            temperature_c,
            pressure_hpa,
            humidity,
            altitude_m,
            0.0,
            0.0,
            end_m * 1.1,
            sample_m,
            zero_angle,
            load.bc_segments_data.clone(),
            load.custom_drag_table.clone(),
        )
        .map_err(|e| format!("load '{}': {e}", load.name))?;

        let mut rows: Vec<LoadRow> = Vec::new();
        for &range_display in &ranges {
            let range_m = UnitConverter::distance_to_metric(range_display, units);
            let nw = no_wind_samples.iter().min_by(|a, b| {
                (a.distance_m - range_m)
                    .abs()
                    .partial_cmp(&(b.distance_m - range_m).abs())
                    .unwrap()
            });
            let w = wind_samples.iter().min_by(|a, b| {
                (a.distance_m - range_m)
                    .abs()
                    .partial_cmp(&(b.distance_m - range_m).abs())
                    .unwrap()
            });
            let (nw, w) = match (nw, w) {
                (Some(nw), Some(w)) => (nw, w),
                _ => {
                    return Err(format!(
                        "load '{}': no trajectory samples near {range_display} \
                         (bullet may not reach --end)",
                        load.name
                    )
                    .into())
                }
            };
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
            rows.push(LoadRow {
                drop_linear,
                drop_adj,
                wind_linear,
                wind_adj,
                velocity: UnitConverter::velocity_from_metric(nw.velocity_mps, units),
                energy: UnitConverter::energy_from_metric(nw.energy_j, units),
                time: nw.time_s,
            });
        }
        results.push(LoadResult {
            zero_angle_deg: zero_angle.to_degrees(),
            rows,
        });
    }

    let adj_label = match adjustment_unit {
        AdjustmentUnit::Mil => "MIL",
        AdjustmentUnit::Moa => "MOA",
        AdjustmentUnit::Smoa => "SMOA",
        AdjustmentUnit::Iphy => "IPHY",
        AdjustmentUnit::Clicks => {
            debug_assert!(false, "Clicks must be resolved via clicks_for() before display");
            "MIL"
        }
    };
    let (dist_unit, vel_unit, energy_unit, drop_unit, wind_unit_label) = match units {
        UnitSystem::Imperial => ("yd", "fps", "ft-lb", "in", "mph"),
        UnitSystem::Metric => ("m", "m/s", "J", "mm", "m/s"),
    };

    match output {
        OutputFormat::Json => {
            let json = serde_json::json!({
                "compare": {
                    "loads": loads.iter().zip(results.iter()).map(|(l, res)| {
                        serde_json::json!({
                            "name": l.name,
                            "drag_model": match l.drag_model { DragModelArg::G1 => "G1", DragModelArg::G7 => "G7" },
                            "bc": l.bc,
                            "mass": l.mass,
                            "velocity": l.velocity,
                            "diameter": l.diameter,
                            "zero_angle_deg": res.zero_angle_deg,
                            "bc_segments": l.bc_segments_data.is_some(),
                            "custom_drag_curve": l.custom_drag_table.is_some(),
                        })
                    }).collect::<Vec<_>>(),
                    "conditions": {
                        "zero_distance": zero_distance,
                        "wind_speed": wind_speed,
                        "wind_direction": wind_direction,
                        "temperature": temperature,
                        "pressure": pressure,
                        "humidity": humidity,
                        "altitude": altitude,
                        "sight_height": sight_height,
                    },
                    "units": {
                        "distance": dist_unit,
                        "drop": drop_unit,
                        "adjustment": adj_label,
                        "velocity": vel_unit,
                        "energy": energy_unit,
                        "wind_speed": wind_unit_label,
                    },
                    "rows": ranges.iter().enumerate().map(|(ri, &range)| {
                        serde_json::json!({
                            "range": range,
                            "loads": results.iter().enumerate().map(|(li, res)| {
                                let row = &res.rows[ri];
                                let baseline = &results[0].rows[ri];
                                serde_json::json!({
                                    "name": loads[li].name,
                                    "drop": row.drop_linear,
                                    "drop_adj": row.drop_adj,
                                    "drift": row.wind_linear,
                                    "drift_adj": row.wind_adj,
                                    "velocity": row.velocity,
                                    "energy": row.energy,
                                    "time": row.time,
                                    // deltas vs load #1 (zero for the baseline itself)
                                    "delta_drop": row.drop_linear - baseline.drop_linear,
                                    "delta_drift": row.wind_linear - baseline.wind_linear,
                                    "delta_velocity": row.velocity - baseline.velocity,
                                    "delta_energy": row.energy - baseline.energy,
                                })
                            }).collect::<Vec<_>>(),
                        })
                    }).collect::<Vec<_>>(),
                }
            });
            println!("{}", serde_json::to_string_pretty(&json)?);
        }
        OutputFormat::Csv => {
            // Sanitize load names for CSV headers (no commas/quotes/newlines)
            let safe: Vec<String> = loads
                .iter()
                .map(|l| {
                    l.name
                        .chars()
                        .map(|c| if c == ',' || c == '"' || c == '\n' { ' ' } else { c })
                        .collect()
                })
                .collect();
            let mut header = format!("range_{dist_unit}");
            for name in &safe {
                header.push_str(&format!(
                    ",{name}_drop_{drop_unit},{name}_drop_{adj_label},{name}_drift_{drop_unit},\
                     {name}_drift_{adj_label},{name}_velocity_{vel_unit},{name}_energy_{energy_unit},\
                     {name}_tof_s"
                ));
            }
            println!("{header}");
            for (ri, &range) in ranges.iter().enumerate() {
                let mut line = format!("{range:.0}");
                for res in &results {
                    let row = &res.rows[ri];
                    line.push_str(&format!(
                        ",{:.2},{:.2},{:.2},{:.2},{:.0},{:.0},{:.3}",
                        row.drop_linear,
                        row.drop_adj,
                        row.wind_linear,
                        row.wind_adj,
                        row.velocity,
                        row.energy,
                        row.time
                    ));
                }
                println!("{line}");
            }
        }
        OutputFormat::Pdf => {
            return Err("PDF output is not supported for compare (use the card commands)".into());
        }
        OutputFormat::Table => {
            println!("Load Comparison");
            println!(
                "Zero: {zero_distance:.0} {dist_unit} | Wind: {wind_speed:.0} {wind_unit_label} @ {wind_direction:.0}\u{b0} | Sight: {sight_height} {sh_unit} | Alt: {altitude:.0} {alt_unit}",
                sh_unit = match units { UnitSystem::Imperial => "in", UnitSystem::Metric => "mm" },
                alt_unit = match units { UnitSystem::Imperial => "ft", UnitSystem::Metric => "m" }
            );
            println!();
            for (i, l) in loads.iter().enumerate() {
                let physics_tag = match (&l.bc_segments_data, &l.custom_drag_table) {
                    (_, Some(_)) => " [custom drag curve]",
                    (Some(_), None) => " [BC segments]",
                    (None, None) => "",
                };
                println!(
                    "  #{num} {name}: {dm} BC {bc}, {mass} @ {vel} {vel_unit} (dia {dia}){physics_tag}",
                    num = i + 1,
                    name = l.name,
                    dm = match l.drag_model {
                        DragModelArg::G1 => "G1",
                        DragModelArg::G7 => "G7",
                    },
                    bc = l.bc,
                    mass = l.mass,
                    vel = l.velocity,
                    dia = l.diameter,
                    physics_tag = physics_tag,
                );
            }
            println!();

            // Two-line header: load names over their [Drop Drift Vel] groups
            let trunc = |s: &str| -> String {
                if s.chars().count() > 26 {
                    let t: String = s.chars().take(25).collect();
                    format!("{t}\u{2026}")
                } else {
                    s.to_string()
                }
            };
            let mut h1 = format!("{:>8} ", "Range");
            let mut h2 = format!("{:>8} ", format!("({dist_unit})"));
            for l in &loads {
                h1.push_str(&format!("| {:^26} ", trunc(&l.name)));
                h2.push_str(&format!("| {:>8} {:>8} {:>8} ", "Drop", "Drift", "Vel"));
            }
            println!("{h1}");
            println!("{h2}");
            let width = 9 + loads.len() * 29;
            println!("{}", "-".repeat(width));
            for (ri, &range) in ranges.iter().enumerate() {
                let mut line = format!("{range:>8.0} ");
                for res in &results {
                    let row = &res.rows[ri];
                    line.push_str(&format!(
                        "| {:>8.2} {:>8.2} {:>8.0} ",
                        row.drop_adj, row.wind_adj, row.velocity
                    ));
                }
                println!("{line}");
            }
            println!("{}", "-".repeat(width));
            println!(
                "Drop/Drift in {adj_label} (- is down/left), Vel in {vel_unit}. \
                 Use -o json or -o csv for linear drop/drift ({drop_unit}), energy \
                 ({energy_unit}), time of flight, and deltas vs load #1."
            );
        }
    }

    Ok(())
}

#[cfg(test)]
mod profile_unit_tests {
    use super::*;

    fn metric_profile() -> ProfileData {
        ProfileData {
            name: "metric-exact".to_string(),
            velocity: 762.0,
            bc: 0.243,
            mass: 11.339_809_25,
            diameter: 7.8232,
            drag_model: "G7".to_string(),
            twist_rate: Some(254.0),
            sight_height: Some(50.8),
            zero_distance: Some(91.44),
            units: "metric".to_string(),
            temperature: 15.0,
            pressure: 1_013.207_888,
            humidity: 55.0,
            altitude: 304.8,
            bullet_name: Some("175gr SMK".to_string()),
            created: Some("test".to_string()),
            wind_speed: Some(4.4704),
            wind_direction: Some(90.0),
            shooting_angle: Some(-5.0),
            auto_zero: Some(91.44),
            twist_right: Some(false),
            use_bc_segments: Some(true),
            bullet_length: Some(30.48),
            elevation_click: None,
            windage_click: None,
            bc_segments: None,
            drag_curve: None,
        }
    }

    fn assert_close(actual: f64, expected: f64) {
        let tolerance = expected.abs().max(1.0) * 1e-12;
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual={actual} expected={expected}"
        );
    }

    #[test]
    fn profile_conversion_preserves_physics_and_unitless_fields() {
        let metric = metric_profile();
        let imperial = metric.clone().converted_to(UnitSystem::Imperial).unwrap();

        assert_eq!(imperial.units, "imperial");
        assert_close(imperial.velocity, 2500.0);
        assert_close(imperial.mass, 175.0);
        assert_close(imperial.diameter, 0.308);
        assert_close(imperial.twist_rate.unwrap(), 10.0);
        assert_close(imperial.sight_height.unwrap(), 2.0);
        assert_close(imperial.zero_distance.unwrap(), 100.0);
        assert_close(imperial.temperature, 59.0);
        assert_close(imperial.pressure, 29.92);
        assert_close(imperial.altitude, 1000.0);
        assert_close(imperial.wind_speed.unwrap(), 10.0);
        assert_close(imperial.auto_zero.unwrap(), 100.0);
        assert_close(imperial.bullet_length.unwrap(), 1.2);
        assert_eq!(imperial.bc.to_bits(), metric.bc.to_bits());
        assert_eq!(imperial.humidity.to_bits(), metric.humidity.to_bits());
        assert_eq!(imperial.wind_direction, metric.wind_direction);
        assert_eq!(imperial.shooting_angle, metric.shooting_angle);
        assert_eq!(imperial.twist_right, metric.twist_right);
        assert_eq!(imperial.use_bc_segments, metric.use_bc_segments);
        assert_eq!(imperial.bullet_name, metric.bullet_name);

        let round_trip = imperial.converted_to(UnitSystem::Metric).unwrap();
        assert_close(round_trip.velocity, metric.velocity);
        assert_close(round_trip.mass, metric.mass);
        assert_close(round_trip.diameter, metric.diameter);
        assert_close(round_trip.twist_rate.unwrap(), metric.twist_rate.unwrap());
        assert_close(
            round_trip.sight_height.unwrap(),
            metric.sight_height.unwrap(),
        );
        assert_close(round_trip.temperature, metric.temperature);
        assert_close(round_trip.pressure, metric.pressure);
        assert_close(round_trip.altitude, metric.altitude);
    }

    #[test]
    fn bc_segments_and_drag_curve_survive_unit_conversion_untouched() {
        // bc_segments/drag_curve are unit-invariant/SI-pinned (see their doc comments on
        // ProfileBcSegment/ProfileDragPoint) — converted_to must be a bit-exact no-op for
        // them in both directions, never rescaling velocity_mps/mach/cd.
        let base = metric_profile();
        let with_v2_fields = ProfileData {
            bc_segments: Some(vec![
                ProfileBcSegment {
                    bc: 0.243,
                    velocity_mps: 792.0,
                },
                ProfileBcSegment {
                    bc: 0.230,
                    velocity_mps: 400.0,
                },
            ]),
            drag_curve: Some(vec![
                ProfileDragPoint { mach: 0.5, cd: 0.23 },
                ProfileDragPoint { mach: 1.2, cd: 0.45 },
                ProfileDragPoint { mach: 3.0, cd: 0.28 },
            ]),
            ..base
        };

        let imperial = with_v2_fields
            .clone()
            .converted_to(UnitSystem::Imperial)
            .unwrap();
        assert_eq!(imperial.bc_segments, with_v2_fields.bc_segments);
        assert_eq!(imperial.drag_curve, with_v2_fields.drag_curve);

        let round_trip = imperial.converted_to(UnitSystem::Metric).unwrap();
        assert_eq!(round_trip.bc_segments, with_v2_fields.bc_segments);
        assert_eq!(round_trip.drag_curve, with_v2_fields.drag_curve);

        // Same-unit short-circuit path (source == target) must also preserve them.
        let same_unit = with_v2_fields.clone().converted_to(UnitSystem::Metric).unwrap();
        assert_eq!(same_unit.bc_segments, with_v2_fields.bc_segments);
        assert_eq!(same_unit.drag_curve, with_v2_fields.drag_curve);
    }

    /// A Phase-1-era profile JSON (no bc_segments/drag_curve keys at all) must still
    /// deserialize cleanly, with both new fields defaulting to None (MBA-1323 Phase 2
    /// backward compatibility).
    #[test]
    fn phase1_profile_json_without_v2_fields_deserializes_with_none() {
        let phase1_json = r#"{
            "name": "phase1-profile",
            "velocity": 2700.0,
            "bc": 0.475,
            "mass": 175.0,
            "diameter": 0.308,
            "drag_model": "G7",
            "units": "imperial",
            "temperature": 59.0,
            "pressure": 29.92,
            "humidity": 50.0,
            "altitude": 0.0
        }"#;
        let profile: ProfileData = serde_json::from_str(phase1_json).unwrap();
        assert_eq!(profile.name, "phase1-profile");
        assert!(profile.bc_segments.is_none());
        assert!(profile.drag_curve.is_none());

        // And re-serializing must NOT introduce the new keys (skip_serializing_if), so a
        // profile untouched by Phase 2 features stays byte-for-byte compatible with tools
        // reading the old schema.
        let reserialized = serde_json::to_string(&profile).unwrap();
        assert!(!reserialized.contains("bc_segments"));
        assert!(!reserialized.contains("drag_curve"));
    }

    #[test]
    fn same_unit_profile_conversion_is_bit_preserving() {
        let profile = metric_profile();
        let before = serde_json::to_string(&profile).unwrap();
        let after = profile.converted_to(UnitSystem::Metric).unwrap();
        assert_eq!(serde_json::to_string(&after).unwrap(), before);
    }

    #[test]
    fn converted_profile_keeps_explicit_cli_precedence() {
        let profile = Some(metric_profile().converted_to(UnitSystem::Imperial).unwrap());
        assert_eq!(
            resolve_param(Some(2600.0), &profile, |p| p.velocity),
            Some(2600.0)
        );
        assert_eq!(resolve_param(None, &profile, |p| p.velocity), Some(2500.0));
        assert_eq!(
            Some(2.5).or_else(|| profile.as_ref().and_then(|p| p.sight_height)),
            Some(2.5)
        );
        assert_eq!(
            None.or_else(|| profile.as_ref().and_then(|p| p.sight_height)),
            Some(2.0)
        );
    }

    #[test]
    fn legacy_profiles_default_to_imperial_and_unknown_units_are_rejected() {
        let legacy: ProfileData = serde_json::from_str(
            r#"{"name":"legacy","velocity":2500.0,"bc":0.243,"mass":175.0,"diameter":0.308,"drag_model":"G7"}"#,
        )
        .unwrap();
        assert_eq!(legacy.units, "imperial");
        let metric = legacy.converted_to(UnitSystem::Metric).unwrap();
        assert_close(metric.velocity, 762.0);
        assert_close(metric.mass, 11.339_809_25);
        assert_close(metric.diameter, 7.8232);

        let mut invalid = metric_profile();
        invalid.units = "si".to_string();
        let error = invalid.converted_to(UnitSystem::Imperial).unwrap_err();
        assert!(error.to_string().contains("unsupported units 'si'"));
    }

    #[test]
    fn profile_click_fields_roundtrip_and_flags_override() {
        let mut profile = metric_profile();
        profile.elevation_click = Some("0.1mil".to_string());
        profile.windage_click = Some("0.2mil".to_string());

        // No flags: the profile's own fields resolve directly.
        let (el, wi) = resolve_click_values(None, None, Some(&profile))
            .unwrap()
            .unwrap();
        assert_eq!(el, parse_click_value("0.1mil").unwrap());
        assert_eq!(wi, parse_click_value("0.2mil").unwrap());

        // An explicit flag beats the saved profile.
        let (el2, wi2) = resolve_click_values(Some("0.25moa"), Some("0.5moa"), Some(&profile))
            .unwrap()
            .unwrap();
        assert_eq!(el2, parse_click_value("0.25moa").unwrap());
        assert_eq!(wi2, parse_click_value("0.5moa").unwrap());

        // Elevation-only profile: windage falls back to the profile's elevation click.
        let mut elev_only = metric_profile();
        elev_only.elevation_click = Some("0.25moa".to_string());
        let (el3, wi3) = resolve_click_values(None, None, Some(&elev_only))
            .unwrap()
            .unwrap();
        assert_eq!(el3, wi3);

        // `elevation_click`/`windage_click` are unit-invariant angular graduations, like
        // bc_segments/drag_curve — converted_to must leave them untouched.
        let imperial = profile.clone().converted_to(UnitSystem::Imperial).unwrap();
        assert_eq!(imperial.elevation_click, profile.elevation_click);
        assert_eq!(imperial.windage_click, profile.windage_click);
    }
}

#[cfg(test)]
#[cfg(feature = "profile-import")]
mod a7p_import_mapping_tests {
    use super::*;
    use ballistics_engine::profile_import::{parse_a7p, wrap_payload};

    // Spec-derived test encoders (kept deliberately independent of the lib).
    fn enc_varint(mut v: u64, out: &mut Vec<u8>) {
        loop {
            let byte = (v & 0x7f) as u8;
            v >>= 7;
            if v == 0 {
                out.push(byte);
                return;
            }
            out.push(byte | 0x80);
        }
    }
    fn enc_i32(number: u32, value: i64, out: &mut Vec<u8>) {
        enc_varint(u64::from(number) << 3, out);
        enc_varint(value as u64, out);
    }
    fn enc_str(number: u32, s: &str, out: &mut Vec<u8>) {
        enc_varint((u64::from(number) << 3) | 2, out);
        enc_varint(s.len() as u64, out);
        out.extend_from_slice(s.as_bytes());
    }
    fn enc_bytes(number: u32, payload: &[u8], out: &mut Vec<u8>) {
        enc_varint((u64::from(number) << 3) | 2, out);
        enc_varint(payload.len() as u64, out);
        out.extend_from_slice(payload);
    }

    fn barnes_338_file(bc_type: i64) -> Vec<u8> {
        let mut p = Vec::new();
        enc_str(1, "338LM/BARNES 300", &mut p); // '/' must be sanitized
        enc_str(3, "BARNES 300GR OTM", &mut p);
        enc_i32(9, 90, &mut p);
        enc_i32(10, 1000, &mut p);
        enc_i32(11, 7920, &mut p);
        enc_i32(12, 15, &mut p);
        enc_i32(13, 1000, &mut p);
        enc_i32(15, 15, &mut p);
        enc_i32(16, 10000, &mut p);
        enc_i32(17, 50, &mut p);
        enc_i32(20, 338, &mut p);
        enc_i32(21, 3000, &mut p);
        enc_i32(22, 1800, &mut p);
        enc_i32(24, bc_type, &mut p);
        let mut packed = Vec::new();
        enc_varint(10_000, &mut packed);
        enc_bytes(26, &packed, &mut p);
        // Three rows, deliberately ordered so max-velocity, max-BC, first-row,
        // and last-row selection strategies all disagree:
        //   (7500, 5000) -> BC 0.750 @ 500 m/s  (highest BC, encoded first)
        //   (7160, 7920) -> BC 0.716 @ 792 m/s  (highest velocity -- the winner)
        //   (6900, 3000) -> BC 0.690 @ 300 m/s  (lowest of everything, encoded last)
        for (bc, mv) in [(7500i64, 5000i64), (7160, 7920), (6900, 3000)] {
            let mut row = Vec::new();
            enc_i32(1, bc, &mut row);
            enc_i32(2, mv, &mut row);
            enc_bytes(27, &row, &mut p);
        }
        let mut payload = Vec::new();
        enc_bytes(1, &p, &mut payload);
        wrap_payload(&payload)
    }

    #[test]
    fn maps_g1_profile_to_metric_profiledata() {
        let doc = parse_a7p(&barnes_338_file(0)).unwrap();
        let outcome = map_a7p_to_profile(&doc, None).unwrap();
        let p = &outcome.profile;
        assert_eq!(p.name, "338LM_BARNES 300"); // sanitized
        assert_eq!(p.units, "metric");
        assert_eq!(p.drag_model, "G1");
        assert!((p.velocity - 792.0).abs() < 1e-9);
        assert!((p.bc - 0.716).abs() < 1e-9); // highest-velocity row wins
        // Independent literal (not GRAIN_TO_GRAM) so this assertion still
        // catches a corrupted grain->gram constant in production code.
        assert!((p.mass - 300.0 * 0.06479891).abs() < 1e-9); // grams
        assert!((p.diameter - 0.338 * 25.4).abs() < 1e-9); // mm
        assert!((p.bullet_length.unwrap() - 1.8 * 25.4).abs() < 1e-9); // mm
        assert!((p.twist_rate.unwrap() - 254.0).abs() < 1e-9); // mm/turn
        assert_eq!(p.twist_right, Some(true));
        assert!((p.sight_height.unwrap() - 90.0).abs() < 1e-9);
        assert!((p.zero_distance.unwrap() - 100.0).abs() < 1e-9);
        assert!((p.temperature - 15.0).abs() < 1e-9);
        assert!((p.pressure - 1000.0).abs() < 1e-9);
        assert!((p.humidity - 50.0).abs() < 1e-9);
        assert_eq!(p.bullet_name.as_deref(), Some("BARNES 300GR OTM"));
        assert!(p.created.is_some());

        // All three BC rows are mapped into bc_segments (velocity-banded, descending),
        // not warned away — MBA-1323 Phase 2.
        let segs = p.bc_segments.as_ref().expect("bc_segments populated");
        assert_eq!(segs.len(), 3);
        let velocities: Vec<f64> = segs.iter().map(|s| s.velocity_mps).collect();
        assert_eq!(velocities, vec![792.0, 500.0, 300.0]); // descending
        assert!((segs[0].bc - 0.716).abs() < 1e-9); // fastest row == primary bc
        assert!((segs[1].bc - 0.750).abs() < 1e-9);
        assert!((segs[2].bc - 0.690).abs() < 1e-9);
        assert!(outcome
            .report
            .mapped
            .iter()
            .any(|[field, ..]| field == "coef_rows[all]"));
        // powder temp sensitivity is honestly unmapped, with the converted value shown
        let tcoeff = outcome
            .report
            .unmapped
            .iter()
            .find(|(field, _)| field == "c_t_coeff")
            .expect("c_t_coeff reported");
        // 792 m/s * (1.0/100) / 15 C = 0.528 m/s per degC
        assert!(tcoeff.1.contains("0.528"));
        // zero-condition ambient temperature (field 12) is honestly unmapped too
        let zero_temp = outcome
            .report
            .unmapped
            .iter()
            .find(|(field, _)| field == "c_zero_temperature")
            .expect("c_zero_temperature reported");
        assert!(zero_temp.1.contains("15 C"));
    }

    #[test]
    fn name_override_is_applied() {
        let doc = parse_a7p(&barnes_338_file(0)).unwrap();
        let outcome = map_a7p_to_profile(&doc, Some("my-338")).unwrap();
        assert_eq!(outcome.profile.name, "my-338");
    }

    #[test]
    fn custom_bc_type_imports_a_drag_curve() {
        // bc_type=2 (CUSTOM) reinterprets the same three coef rows as (Cd, Mach) pairs
        // instead of (BC, velocity m/s) — see barnes_338_file's row comment.
        let doc = parse_a7p(&barnes_338_file(2)).unwrap();
        let outcome = map_a7p_to_profile(&doc, None).unwrap();
        let p = &outcome.profile;

        assert_eq!(p.drag_model, "CUSTOM");
        assert_eq!(p.bc, 0.0); // documented inert sentinel, not a real coefficient
        assert!(p.bc_segments.is_none());

        let curve = p.drag_curve.as_ref().expect("drag_curve populated");
        assert_eq!(curve.len(), 3);
        // sorted ascending by Mach (DragTable requirement), NOT file order.
        let machs: Vec<f64> = curve.iter().map(|pt| pt.mach).collect();
        assert_eq!(machs, vec![0.3, 0.5, 0.792]);
        assert!((curve[0].cd - 0.69).abs() < 1e-9);
        assert!((curve[1].cd - 0.75).abs() < 1e-9);
        assert!((curve[2].cd - 0.716).abs() < 1e-9);

        // The curve validates as a real engine DragTable.
        let mach_values: Vec<f64> = curve.iter().map(|pt| pt.mach).collect();
        let cd_values: Vec<f64> = curve.iter().map(|pt| pt.cd).collect();
        assert!(ballistics_engine::drag::DragTable::try_new(mach_values, cd_values).is_ok());

        assert!(outcome
            .report
            .warnings
            .iter()
            .any(|w| w.contains("inert sentinel")));
        assert!(outcome
            .report
            .mapped
            .iter()
            .any(|[field, ..]| field == "coef_rows[custom]"));
    }

    #[test]
    fn custom_bc_type_rejects_a_degenerate_curve() {
        // A single coef row cannot become a 2-point DragTable (try_new requires >= 2).
        let mut p = Vec::new();
        enc_i32(24, 2, &mut p); // CUSTOM
        let mut row = Vec::new();
        enc_i32(1, 5000, &mut row);
        enc_i32(2, 5000, &mut row);
        enc_bytes(27, &row, &mut p);
        let mut payload = Vec::new();
        enc_bytes(1, &p, &mut payload);
        let doc = parse_a7p(&wrap_payload(&payload)).unwrap();
        let err = map_a7p_to_profile(&doc, None).unwrap_err();
        assert!(err.contains("CUSTOM drag curve is invalid"), "{err}");
    }

    #[test]
    fn g7_maps_and_report_renders() {
        let doc = parse_a7p(&barnes_338_file(1)).unwrap();
        let outcome = map_a7p_to_profile(&doc, None).unwrap();
        assert_eq!(outcome.profile.drag_model, "G7");
        let rendered = render_import_report(&outcome.report);
        assert!(rendered.contains("muzzle velocity"));
        assert!(rendered.contains("NOT imported"));
    }

    #[test]
    fn missing_bc_rows_is_an_error() {
        // A file with no coef_rows cannot become a ProfileData (bc is mandatory).
        let mut p = Vec::new();
        enc_i32(11, 7920, &mut p);
        let mut payload = Vec::new();
        enc_bytes(1, &p, &mut payload);
        let doc = parse_a7p(&wrap_payload(&payload)).unwrap();
        let err = map_a7p_to_profile(&doc, None).unwrap_err();
        assert!(err.contains("no BC rows"), "{err}");
    }
}

#[cfg(test)]
mod adjustment_unit_tests {
    use super::*;

    #[test]
    fn mil_conversion_matches_geometry() {
        // 0.1 yd of drop at 100 yd = 1.0 MIL (1 MIL subtends 0.1 yd at 100 yd).
        assert!((drop_to_adjustment(0.1, 100.0, AdjustmentUnit::Mil) - 1.0).abs() < 1e-9);
        // 1.78 yd at 500 yd ~ 3.56 MIL.
        assert!((drop_to_adjustment(1.78, 500.0, AdjustmentUnit::Mil) - 3.56).abs() < 0.01);
    }

    #[test]
    fn moa_is_3438_thousandths_per_radian_ratio() {
        // MOA uses 3438 (arcminutes per radian) vs MIL's 1000, so MOA/MIL == 3.438.
        let mil = drop_to_adjustment(1.0, 100.0, AdjustmentUnit::Mil);
        let moa = drop_to_adjustment(1.0, 100.0, AdjustmentUnit::Moa);
        assert!((moa / mil - 3.438).abs() < 1e-9, "moa={moa} mil={mil}");
        // 0.1 yd at 100 yd -> 3.438 MOA.
        assert!((drop_to_adjustment(0.1, 100.0, AdjustmentUnit::Moa) - 3.438).abs() < 1e-9);
    }

    #[test]
    fn short_range_returns_zero() {
        assert_eq!(drop_to_adjustment(1.0, 0.5, AdjustmentUnit::Moa), 0.0);
        assert_eq!(drop_to_adjustment(1.0, 0.5, AdjustmentUnit::Mil), 0.0);
    }

    #[test]
    fn smoa_iphy_conversion_is_inches_per_hundred_yards() {
        // 3.6 inches of drop at 100 yd = 0.1 yd drop over 100 yd... use exact math:
        // drop_yd/range_yd * 3600 == inches per 100 yd.
        let v = drop_to_adjustment(0.1, 100.0, AdjustmentUnit::Smoa);
        assert!((v - 3.6).abs() < 1e-12, "0.1 yd @ 100 yd = 3.6 IPHY, got {v}");
        let i = drop_to_adjustment(0.1, 100.0, AdjustmentUnit::Iphy);
        assert_eq!(v, i, "smoa and iphy are the same unit");
        // sanity vs existing: TMOA value must be larger number x smaller unit
        let t = drop_to_adjustment(0.1, 100.0, AdjustmentUnit::Moa);
        assert!(t < v && (v / t - 1.047).abs() < 0.005, "tmoa {t} vs smoa {v}");
    }

    #[test]
    fn clicks_unit_requires_elevation_graduation() {
        // resolve_click_values with no flag and no profile must error, naming the flag.
        let e = resolve_click_values(None, None, None).unwrap_err();
        assert!(e.contains("--elevation-click-value"), "{e}");
    }

    #[test]
    fn windage_click_defaults_to_elevation() {
        let (el, wi) = resolve_click_values(Some("0.25moa"), None, None)
            .unwrap()
            .unwrap();
        assert_eq!(el, wi);
    }

    #[test]
    fn explicit_windage_flag_overrides_elevation_fallback() {
        let (el, wi) = resolve_click_values(Some("0.25moa"), Some("0.1mil"), None)
            .unwrap()
            .unwrap();
        assert_ne!(el, wi);
        assert_eq!(wi, parse_click_value("0.1mil").unwrap());
    }

    #[test]
    fn default_trajectory_header_is_stable() {
        // Pin the come-ups table header line for the default unit (Mil, imperial) so a
        // future change to the Clicks formatting branch can't silently reformat the
        // existing MIL/MOA/SMOA/IPHY header text.
        assert_eq!(
            come_up_header_line("yd", "MIL", "fps"),
            "│Range (yd)|Drop (MIL)|Come-Up   │ Vel (fps)│Energy    │ Time (s) │"
        );
    }
}

#[cfg(test)]
mod drag_model_arg_warning_tests {
    use super::*;

    #[test]
    fn drag_model_arg_warning_flags_silent_g1_coercion() {
        assert!(drag_model_arg_warning("g1").is_none());
        assert!(drag_model_arg_warning("G7").is_none());
        for s in ["g5", "G8", "banana"] {
            let w = drag_model_arg_warning(s).expect("must warn");
            assert!(w.contains("using G1"), "{w}");
        }
    }

    #[test]
    fn parse_drag_model_arg_still_coerces_to_g1_or_g7() {
        // The warning is purely informational (stderr) — behavior is unchanged.
        assert_eq!(parse_drag_model_arg("G7"), DragModelArg::G7);
        assert_eq!(parse_drag_model_arg("g7"), DragModelArg::G7);
        assert_eq!(parse_drag_model_arg("G1"), DragModelArg::G1);
        assert_eq!(parse_drag_model_arg("G5"), DragModelArg::G1);
        assert_eq!(parse_drag_model_arg("banana"), DragModelArg::G1);
    }
}

#[cfg(test)]
mod bc_segment_parse_tests {
    use super::*;

    #[test]
    fn resolved_bc_schedule_preserves_explicit_precedence_and_estimator_fallback() {
        let explicit = vec![BCSegmentData {
            velocity_min: 0.0,
            velocity_max: 5_000.0,
            bc_value: 0.42,
        }];
        let provided = resolve_velocity_bc_segments(
            Some(&explicit),
            true,
            0.19,
            0.224 * 0.0254,
            77.0 * GRAINS_TO_KG,
            DragModelArg::G7,
        )
        .unwrap();
        assert_eq!(provided.len(), 1);
        assert_eq!(provided[0].bc_value.to_bits(), 0.42_f64.to_bits());

        let estimated = resolve_velocity_bc_segments(
            None,
            true,
            0.19,
            0.224 * 0.0254,
            77.0 * GRAINS_TO_KG,
            DragModelArg::G7,
        )
        .unwrap();
        let expected = ballistics_engine::bc_estimation::BCSegmentEstimator::estimate_bc_segments(
            0.19, 0.224, 77.0, "", "G7",
        );
        assert_eq!(estimated.len(), expected.len());
        for (actual, expected) in estimated.iter().zip(&expected) {
            assert_eq!(
                actual.velocity_min.to_bits(),
                expected.velocity_min.to_bits()
            );
            assert_eq!(
                actual.velocity_max.to_bits(),
                expected.velocity_max.to_bits()
            );
            assert_eq!(actual.bc_value.to_bits(), expected.bc_value.to_bits());
        }

        assert!(resolve_velocity_bc_segments(
            None,
            false,
            0.19,
            0.224 * 0.0254,
            77.0 * GRAINS_TO_KG,
            DragModelArg::G7,
        )
        .is_none());
    }

    #[test]
    fn parses_valid_imperial_pair_as_fps() {
        let seg = parse_bc_segment("1500:1800:0.243", UnitSystem::Imperial).unwrap();
        assert_eq!(seg.velocity_min, 1500.0);
        assert_eq!(seg.velocity_max, 1800.0);
        assert!((seg.bc_value - 0.243).abs() < 1e-12);
    }

    #[test]
    fn converts_metric_velocity_to_fps() {
        // 305 m/s -> 1000.66 fps, 427 m/s -> 1400.92 fps; BC unchanged.
        let seg = parse_bc_segment("305:427:0.15", UnitSystem::Metric).unwrap();
        assert!((seg.velocity_min - 305.0 * 3.280_839_895).abs() < 1e-6);
        assert!((seg.velocity_max - 427.0 * 3.280_839_895).abs() < 1e-6);
        assert!((seg.bc_value - 0.15).abs() < 1e-12);
    }

    #[test]
    fn rejects_wrong_field_count() {
        assert!(parse_bc_segment("1000:1400", UnitSystem::Imperial).is_err());
        assert!(parse_bc_segment("1000:1400:0.2:extra", UnitSystem::Imperial).is_err());
    }

    #[test]
    fn rejects_vmin_not_less_than_vmax() {
        assert!(parse_bc_segment("1400:1000:0.2", UnitSystem::Imperial).is_err());
        assert!(parse_bc_segment("1400:1400:0.2", UnitSystem::Imperial).is_err());
    }

    #[test]
    fn rejects_nonpositive_bc_and_nonnumeric() {
        assert!(parse_bc_segment("1000:1400:0", UnitSystem::Imperial).is_err());
        assert!(parse_bc_segment("1000:1400:-0.1", UnitSystem::Imperial).is_err());
        assert!(parse_bc_segment("abc:1400:0.2", UnitSystem::Imperial).is_err());
        assert!(parse_bc_segment("1000:xyz:0.2", UnitSystem::Imperial).is_err());
        assert!(parse_bc_segment("1000:1400:bc", UnitSystem::Imperial).is_err());
    }
}

#[cfg(test)]
mod profile_bc_segments_and_drag_curve_consumption_tests {
    use super::*;

    /// bc_segments_from_profile must produce a gapless, non-overlapping partition of
    /// [0, +inf) — the failure mode of an earlier draft was an overlap between the fastest
    /// and slowest bands (each half-open band's [min, max) must exactly abut its neighbors).
    #[test]
    fn bands_are_gapless_and_non_overlapping_for_two_rows() {
        let rows = vec![
            ProfileBcSegment {
                bc: 0.5,
                velocity_mps: 900.0,
            },
            ProfileBcSegment {
                bc: 0.2,
                velocity_mps: 300.0,
            },
        ];
        let segs = bc_segments_from_profile(&rows);
        assert_eq!(segs.len(), 2);

        let fast = &segs[0];
        let slow = &segs[1];
        assert!((fast.bc_value - 0.5).abs() < 1e-9);
        assert!((slow.bc_value - 0.2).abs() < 1e-9);

        // Slowest band starts at absolute zero, not at its own stated velocity.
        assert_eq!(slow.velocity_min, 0.0);
        // The two bands share exactly one boundary (no gap, no overlap).
        assert_eq!(slow.velocity_max, fast.velocity_min);
        // Fastest band's own velocity is its floor (extends up to/above muzzle), not a
        // midpoint or the other row's velocity.
        let expected_floor_fps = 900.0 * 3.280_839_895;
        assert!((fast.velocity_min - expected_floor_fps).abs() < 1e-6);
        assert!(fast.velocity_max > expected_floor_fps);

        // Sample a velocity from each band and confirm exactly one segment claims it.
        for probe in [fast.velocity_min + 10.0, slow.velocity_max - 10.0] {
            let hits = segs
                .iter()
                .filter(|s| probe >= s.velocity_min && probe < s.velocity_max)
                .count();
            assert_eq!(hits, 1, "velocity {probe} fps must be claimed by exactly one band");
        }
    }

    #[test]
    fn single_row_covers_the_entire_velocity_axis() {
        let rows = vec![ProfileBcSegment {
            bc: 0.3,
            velocity_mps: 800.0,
        }];
        let segs = bc_segments_from_profile(&rows);
        assert_eq!(segs.len(), 1);
        assert_eq!(segs[0].velocity_min, 0.0);
        assert!(segs[0].velocity_max > 800.0 * 3.280_839_895);
    }

    #[test]
    fn unsorted_input_rows_are_sorted_by_velocity() {
        // File order deliberately reversed (slowest first) — output must still be
        // fastest-first with correct banding, independent of input order.
        let rows = vec![
            ProfileBcSegment {
                bc: 0.2,
                velocity_mps: 300.0,
            },
            ProfileBcSegment {
                bc: 0.5,
                velocity_mps: 900.0,
            },
        ];
        let segs = bc_segments_from_profile(&rows);
        assert!((segs[0].bc_value - 0.5).abs() < 1e-9);
        assert!((segs[1].bc_value - 0.2).abs() < 1e-9);
    }

    #[test]
    fn drag_table_from_profile_builds_a_valid_table() {
        let points = vec![
            ProfileDragPoint { mach: 0.5, cd: 0.23 },
            ProfileDragPoint { mach: 1.2, cd: 0.45 },
            ProfileDragPoint { mach: 3.0, cd: 0.28 },
        ];
        let table = drag_table_from_profile(&points).expect("valid table");
        assert_eq!(table.mach_values, vec![0.5, 1.2, 3.0]);
        assert_eq!(table.cd_values, vec![0.23, 0.45, 0.28]);
    }

    #[test]
    fn drag_table_from_profile_rejects_degenerate_input() {
        let points = vec![ProfileDragPoint { mach: 0.5, cd: 0.23 }];
        assert!(drag_table_from_profile(&points).is_err());
    }
}

#[cfg(test)]
mod wind_angle_unit_tests {
    use super::*;

    #[test]
    fn trajectory_components_store_wind_direction_in_radians() {
        let (inputs, wind, _) = build_trajectory_components(
            800.0,
            0.4,
            0.01,
            0.00762,
            DragModelArg::G1,
            0.05,
            15.0,
            1013.25,
            50.0,
            0.0,
            10.0,
            90.0,
            1000.0,
            100.0,
            None,
            None,
        );

        assert!((inputs.wind_angle - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
        assert_eq!(inputs.wind_angle.to_bits(), wind.direction.to_bits());
    }
}
