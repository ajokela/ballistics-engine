// CLI API module - provides simplified interfaces for command-line tool
use crate::cluster_bc::ClusterBCDegradation;
use crate::pitch_damping::{calculate_pitch_damping_coefficient, PitchDampingCoefficients};
use crate::precession_nutation::{
    calculate_combined_angular_motion, projectile_moments_of_inertia, AngularState,
    PrecessionNutationParams,
};
use crate::trajectory_sampling::{
    projected_sample_count, sample_trajectory, TrajectoryData, TrajectoryOutputs,
    TrajectorySample,
};
use crate::trajectory_observation::TrajectoryTermination;
use crate::wind_shear::WindShearModel;
use crate::DragModel;
use nalgebra::{Vector3, Vector6};
use std::error::Error;
use std::fmt;

/// Unit system for CLI-style inputs and outputs.
///
/// The single crate-wide unit-system selector, shared by the CLI binary
/// (`--units`), the truing core ([`crate::truing`]) and the WEZ sweep core
/// ([`crate::wez`]). It only selects how user-facing quantities are
/// interpreted and displayed; the solver itself always works in SI.
///
/// Variant order (Metric first) is load-bearing for the CLI: clap lists
/// `--units` possible values in declaration order.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
pub enum UnitSystem {
    /// Metric units (velocity in m/s, mass in grams, distance in meters, diameter in mm, Celsius)
    Metric,
    /// Imperial units (velocity in fps, mass in grains, distance in yards, diameter in inches, Fahrenheit)
    Imperial,
}

// Output format for results
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OutputFormat {
    Table,
    Json,
    Csv,
}

/// Which standard atmosphere a ballistic coefficient's numeric value is referenced to
/// (MBA-1365).
///
/// A BC is not a pure property of the bullet — it also encodes an assumed reference air
/// density. This engine's retardation constant ([`crate::constants::CD_TO_RETARD`]) is
/// calibrated against the ICAO Standard Atmosphere's sea-level density
/// ([`crate::constants::ICAO_DENSITY_LB_FT3`]), which most modern published BCs already
/// assume. Some vendors (notably Sierra, Hornady, and Barnes for many bullets) instead
/// publish BCs referenced to the older, denser Army Standard Metro atmosphere
/// ([`crate::constants::ASM_DENSITY_LB_FT3`]). Feeding an ASM-referenced BC into this
/// engine as if it were ICAO-referenced under-predicts drag by about 1.8% — declaring
/// `ArmyStandardMetro` here corrects for that exactly once, at input normalization (see
/// [`crate::constants::ASM_TO_ICAO_BC`] and `TrajectorySolver::new`).
///
/// Does not apply when a custom drag table (`BallisticInputs::custom_drag_table`) is
/// active: that path divides by sectional density, not a BC value, so no BC reference
/// conversion is physically meaningful there (see
/// `BallisticInputs::bc_reference_standard_inert_warning`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
pub enum BcReferenceStandard {
    /// The ICAO Standard Atmosphere reference density — matches this engine's
    /// retardation constant. The default; byte-identical to pre-MBA-1365 behavior.
    #[default]
    Icao,
    /// The (older) Army Standard Metro reference density used by some vendor-published
    /// BCs (e.g. many Sierra/Hornady/Barnes bullets).
    ArmyStandardMetro,
}

/// Shared text for the MBA-1365 "`--bc-reference army-standard-metro` is inert" warning —
/// single source of truth for every surface that can trigger it (native CLI, WASM,
/// [`BallisticInputs::bc_reference_standard_inert_warning`]) so the wording can't drift.
pub const BC_REFERENCE_STANDARD_INERT_WARNING: &str =
    "warning: --bc-reference army-standard-metro has no effect together with a custom drag \
     table (--drag-table): the deck's Cd is divided by sectional density, not a BC value, so \
     no BC-reference conversion applies";

// Error type for CLI operations
#[derive(Debug)]
pub struct BallisticsError {
    message: String,
}

impl fmt::Display for BallisticsError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl Error for BallisticsError {}

impl From<String> for BallisticsError {
    fn from(msg: String) -> Self {
        BallisticsError { message: msg }
    }
}

impl From<&str> for BallisticsError {
    fn from(msg: &str) -> Self {
        BallisticsError {
            message: msg.to_string(),
        }
    }
}

/// Which plane sampled drop values are referenced to (MBA-1403).
///
/// See [`BallisticInputs::drops_reference`]. This is an output-mode toggle for the
/// trajectory sampler, not new plane machinery: the solver implements incline as gravity
/// rotation into the shot frame, so sampled drop is already perpendicular to the LOS and
/// the JBM-equivalent "target" mode is a `1 / cos(shooting_angle)` reference transform
/// plus relabeling.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum DropsReference {
    /// Drop measured perpendicular to the line of sight (the historical default).
    #[default]
    Los,
    /// Drop measured vertically in the target plane (JBM's "target plane" reference):
    /// the LOS-perpendicular drop scaled by `1 / cos(shooting_angle)`.
    Target,
}

// Ballistic input parameters - MBA-151 Reconciled Structure
// Unified structure used by both ballistics-engine and ballistics_rust
// Duplicates removed, all necessary fields included
#[derive(Debug, Clone)]
pub struct BallisticInputs {
    // Core ballistics parameters (using intuitive names)
    pub bc_value: f64,        // Ballistic coefficient (G1, G7, etc.)
    pub bc_type: DragModel,   // Drag model (G1, G7, G8, etc.)
    /// Which standard atmosphere `bc_value`/`bc_segments`/`bc_segments_data` are
    /// referenced to (MBA-1365). `Icao` (the default) is a no-op; `ArmyStandardMetro`
    /// is converted to the ICAO reference exactly once, in `TrajectorySolver::new`,
    /// before any retardation math runs. See [`BcReferenceStandard`].
    pub bc_reference_standard: BcReferenceStandard,
    pub bullet_mass: f64,     // kg
    pub muzzle_velocity: f64, // m/s
    pub bullet_diameter: f64, // meters
    pub bullet_length: f64,   // meters

    // Targeting and positioning
    pub muzzle_angle: f64,    // radians (launch angle)
    pub target_distance: f64, // meters
    pub azimuth_angle: f64, // horizontal aiming angle in radians (small aim offset within the shot frame)
    /// Compass bearing the shot is fired ALONG, radians, 0 = North, π/2 = East.
    /// Used only by the Coriolis model (Earth-rotation depends on which way downrange
    /// points relative to true North). Distinct from `azimuth_angle`, which is the
    /// small horizontal *aiming* offset and rotates the launch velocity.
    pub shot_azimuth: f64,
    pub shooting_angle: f64,   // uphill/downhill angle in radians
    /// Rifle cant angle in radians about the line of sight — positive = clockwise from the
    /// shooter's view (top of the scope tips right). Rotates the sight-frame aim offsets
    /// (`muzzle_angle`, `azimuth_angle`) about the LOS and swings the bore's sight-height
    /// offset laterally, producing the classic canted-rifle POI error (right-and-low for
    /// clockwise cant with an upward zero). Zeroing always solves un-canted ("zero level,
    /// fire canted"). NOTE: treats `muzzle_angle` as a sight-frame offset — the standard
    /// zero-then-fire usage; a raw gravity-frame launch angle would not rotate physically.
    /// 0.0 = level rifle (bit-identical to pre-cant behavior). (MBA-1286)
    pub cant_angle: f64,
    pub sight_height: f64,     // meters above bore
    /// Lateral offset between the sight axis and the bore axis, meters (MBA-1396;
    /// offset-mounted optics): positive = the sight sits RIGHT of the bore, so the bore
    /// starts LEFT of the line of sight (initial lateral position `z -= offset`). When a
    /// zero is solved, a windage-zero convergence term of `offset / zero_distance` is
    /// added to `azimuth_angle` (see [`BallisticInputs::windage_zero_bias_rad`]) so the
    /// trajectory crosses the LOS laterally at the zero range — matching AB/JBM/Shooter,
    /// not a constant parallel offset. Without a zero solve (an explicit muzzle angle),
    /// only the physical displacement applies. Distinct from `zero_poi_horizontal_m`,
    /// which is an angular ZERO-STATE bias, not mount geometry. 0.0 (the default) is
    /// byte-identical to pre-MBA-1396 behavior.
    pub sight_offset_lateral_m: f64,
    pub muzzle_height: f64,    // meters above ground
    pub target_height: f64,    // meters above ground for zeroing
    /// Deliberate vertical point-of-impact offset AT THE ZERO RANGE, meters (MBA-1359;
    /// Kestrel "zero height" semantics): positive = the rifle is deliberately zeroed to
    /// impact HIGH by this much at the zero distance. Applied POST-solve by
    /// `calculate_and_set_zero_angle` as a constant angular bias of
    /// `zero_poi_vertical_m / zero_distance` on the solved elevation — the zero trials
    /// themselves always solve perfect convergence, and the small-angle addition is
    /// identical in both zero target frames. 0.0 (the default) is byte-identical to
    /// pre-MBA-1359 behavior.
    pub zero_poi_vertical_m: f64,
    /// Deliberate horizontal point-of-impact offset AT THE ZERO RANGE, meters (MBA-1359;
    /// Kestrel "zero offset" semantics): positive = impacts RIGHT by this much at the
    /// zero distance. Applied POST-solve as an azimuth bias of
    /// `zero_poi_horizontal_m / zero_distance` when a zero is solved — see
    /// [`BallisticInputs::windage_zero_bias_rad`]. 0.0 (the default) is neutral.
    pub zero_poi_horizontal_m: f64,
    pub ground_threshold: f64, // meters below which to stop

    // Environmental conditions
    pub altitude: f64,    // meters
    pub temperature: f64, // Celsius
    pub pressure: f64,    // millibars/hPa
    /// Relative humidity as a FRACTION in `[0, 1]` (e.g. 0.5 = 50%). NOTE the scale
    /// differs from [`AtmosphericConditions::humidity`], which is a PERCENT in `[0, 100]`.
    /// The atmosphere helpers (`calculate_air_density_*`) expect percent, so convert via
    /// [`BallisticInputs::humidity_percent`] before passing this value to them (MBA-722).
    pub humidity: f64,
    pub latitude: Option<f64>, // degrees

    // Wind conditions
    pub wind_speed: f64, // m/s
    pub wind_angle: f64, // radians (0=headwind, PI/2=from right)

    // Bullet characteristics
    pub twist_rate: f64,               // inches per turn
    pub is_twist_right: bool,          // right-hand twist
    pub caliber_inches: f64,           // diameter in inches
    pub weight_grains: f64,            // mass in grains
    pub manufacturer: Option<String>,  // Bullet manufacturer
    pub bullet_model: Option<String>,  // Bullet model name
    pub bullet_id: Option<String>,     // Unique bullet identifier
    pub bullet_cluster: Option<usize>, // BC cluster ID for cluster_bc module

    // Integration method selection
    pub use_rk4: bool,           // Use RK4 integration instead of Euler
    pub use_adaptive_rk45: bool, // Use RK45 adaptive step size integration

    // Advanced effects flags
    pub enable_advanced_effects: bool,
    pub enable_magnus: bool,   // Magnus force (independent of Coriolis)
    pub enable_coriolis: bool, // Coriolis deflection (requires latitude)
    pub use_powder_sensitivity: bool,
    pub powder_temp_sensitivity: f64, // m/s per degree Celsius
    pub powder_temp: f64,           // Celsius
    /// Optional measured powder-temperature -> muzzle-velocity curve, as
    /// (temperature_celsius, muzzle_velocity_m_s) points sorted ascending by
    /// temperature. When present it supersedes the linear `powder_temp_sensitivity`
    /// model: the muzzle velocity is interpolated from this table at the ambient
    /// `temperature` (clamped to the endpoints — no extrapolation beyond measured
    /// data). This is the data-driven, non-linear alternative to the constant slope.
    pub powder_temp_curve: Option<Vec<(f64, f64)>>,
    /// Temperature (Celsius) at which to interpolate `powder_temp_curve` — the POWDER
    /// temperature, which may differ from the ambient `temperature` (air). `None` uses
    /// `temperature`. Decouples the velocity lookup from the air-density temperature.
    pub powder_curve_temp_c: Option<f64>,
    pub tipoff_yaw: f64,            // radians
    pub tipoff_decay_distance: f64, // meters
    /// Enables velocity-keyed `bc_segments_data`. Explicit Mach-keyed `bc_segments` retain their
    /// legacy behavior and remain active when this flag is false.
    pub use_bc_segments: bool,
    pub bc_segments: Option<Vec<(f64, f64)>>, // Mach-BC pairs
    pub bc_segments_data: Option<Vec<crate::BCSegmentData>>, // Velocity-BC segments
    pub use_enhanced_spin_drift: bool,
    /// Legacy compatibility flag. Name-derived "form factors" are intentionally not multiplied
    /// into reference Cd when `bc_value` is already the retardation denominator (MBA-1184).
    pub use_form_factor: bool,
    pub enable_wind_shear: bool,
    pub wind_shear_model: String,
    pub enable_trajectory_sampling: bool,
    pub sample_interval: f64, // meters
    /// Which plane sampled drop values are referenced to (MBA-1403). The default
    /// [`DropsReference::Los`] keeps the historical behavior (drop perpendicular to the
    /// line of sight) byte-identical; [`DropsReference::Target`] reports drop as vertical
    /// in the target plane — the LOS-perpendicular drop scaled by
    /// `1 / cos(shooting_angle)` (JBM's "target plane" checkbox) — and, when a non-zero
    /// `target_height` is supplied, slopes the sampler's LOS datum toward it. This is an
    /// OUTPUT-mode toggle: the solved trajectory itself is unchanged, and zero-solve
    /// `target_height` semantics are untouched.
    pub drops_reference: DropsReference,
    pub enable_pitch_damping: bool,
    pub enable_precession_nutation: bool,
    // MBA-959: apply aerodynamic jump as a muzzle launch-angle perturbation.
    // EXPERIMENTAL — the underlying model is heuristic and not yet validated; default OFF.
    pub enable_aerodynamic_jump: bool,
    pub use_cluster_bc: bool, // Use cluster-based BC degradation

    // Custom drag model support
    pub custom_drag_table: Option<crate::drag::DragTable>,
    /// Whole-curve multiplier applied to the custom deck's interpolated Cd (MBA-1356):
    /// `Cd_used = table.interpolate(mach) * cd_scale`. Meaningful only alongside
    /// `custom_drag_table` — it is read at the three custom-deck interpolation sites in
    /// cli_api.rs/derivatives.rs/fast_trajectory.rs and left untouched on the standard
    /// G-model/BC path (those users true their drag via `--bc-adjustment` instead; scaling
    /// a reference table too would double-count). Default `1.0` is neutral (byte-identical
    /// to the pre-MBA-1356 behavior). Validated finite and > 0 by `validate_for_solve`.
    pub cd_scale: f64,

    // Legacy field for compatibility
    pub bc_type_str: Option<String>,
}

impl BallisticInputs {
    /// `humidity` as a PERCENT in `[0, 100]`, clamped — the scale the atmosphere
    /// density helpers expect. Centralizes the 0–1 → 0–100 conversion so callers don't
    /// re-derive it (and can't accidentally feed the raw 0–1 fraction as a percentage).
    /// See the field doc on [`BallisticInputs::humidity`] (MBA-722).
    pub fn humidity_percent(&self) -> f64 {
        (self.humidity * 100.0).clamp(0.0, 100.0)
    }

    /// Azimuth bias (radians) a zero solved at `zero_distance_m` applies to the launch
    /// direction. Two independent terms share this one convergence point:
    /// * MBA-1359 — the deliberate horizontal POI offset at the zero range
    ///   (`zero_poi_horizontal_m`), an angular zero-state bias;
    /// * MBA-1396 — the lateral sight-mount offset (`sight_offset_lateral_m`): the bore
    ///   starts `offset` left of the LOS (see `initial_position`), so the windage zero
    ///   steers `offset / zero_distance` right to cross the LOS at the zero range.
    ///
    /// `calculate_and_set_zero_angle` applies this to its OWN solver's `azimuth_angle`;
    /// callers that copy a solved zero angle onto separate flight inputs (the CLI/WASM
    /// auto-zero paths) must add this to their flight `azimuth_angle` themselves — the
    /// returned elevation angle cannot carry it. Returns exactly `0.0` when both offsets
    /// are `0.0` or `zero_distance_m` is not positive, so default inputs stay
    /// byte-identical.
    pub fn windage_zero_bias_rad(&self, zero_distance_m: f64) -> f64 {
        if zero_distance_m > 0.0 {
            (self.zero_poi_horizontal_m + self.sight_offset_lateral_m) / zero_distance_m
        } else {
            0.0
        }
    }

    /// Sectional density in lb/in²: `weight_grains / 7000 / diameter_in²`.
    ///
    /// Derived from the imperial mirror fields (`weight_grains` / `caliber_inches`), falling
    /// back to the SI `bullet_mass` (kg) / `bullet_diameter` (meters) for SI-only callers
    /// (mirrors the fallbacks in derivatives.rs). `None` when neither source is usable.
    pub fn sectional_density_lb_in2(&self) -> Option<f64> {
        let weight_gr = if self.weight_grains > 0.0 {
            self.weight_grains
        } else {
            self.bullet_mass / crate::constants::GRAINS_TO_KG // kg -> grains
        };
        let diameter_in = if self.caliber_inches > 0.0 {
            self.caliber_inches
        } else {
            self.bullet_diameter / 0.0254 // meters -> inches
        };
        if weight_gr > 0.0 && diameter_in > 0.0 {
            Some(weight_gr / 7000.0 / (diameter_in * diameter_in))
        } else {
            None
        }
    }

    /// Retardation denominator to use when `custom_drag_table` is active.
    ///
    /// A custom drag table supplies the projectile's ACTUAL drag coefficient, so the
    /// point-mass retardation formula must divide it by the projectile's SECTIONAL DENSITY
    /// (lb/in²), not by a ballistic coefficient: BC = SD / i (form factor i vs the reference
    /// projectile), and with the projectile's own curve i == 1, so Cd_own / SD == Cd_ref / BC.
    /// Dividing the curve's Cd by `bc_value` made custom-table trajectories wrongly scale
    /// with whatever BC happened to be set.
    ///
    /// Falls back to `fallback_bc` (with a one-time stderr warning) when mass/diameter are
    /// unavailable, so degenerate inputs degrade to the old behavior instead of panicking.
    pub fn custom_drag_denominator(&self, fallback_bc: f64) -> f64 {
        match self.sectional_density_lb_in2() {
            Some(sd) => sd,
            None => {
                static WARN_ONCE: std::sync::Once = std::sync::Once::new();
                WARN_ONCE.call_once(|| {
                    eprintln!(
                        "Warning: custom drag table active but bullet mass/diameter are \
                         unavailable; falling back to bc_value for the retardation denominator"
                    );
                });
                fallback_bc
            }
        }
    }

    /// Whether a declared `bc_reference_standard` is inert on this configuration
    /// (MBA-1365), and if so, the message explaining why.
    ///
    /// A custom drag table (`custom_drag_table`) supplies the projectile's actual Cd and
    /// divides by sectional density, not by a BC value (see `custom_drag_denominator`),
    /// so an `ArmyStandardMetro` reference has no physical effect once a table is active.
    /// Pure and side-effect-free — same shape as `main.rs`'s `adjustment_unit_noop_warning`
    /// (MBA-1414) — so each surface decides how to surface it: the native CLI prints it to
    /// stderr once per run, while the WASM browser terminal (no visible stderr) splices it
    /// into its table-output text instead.
    pub fn bc_reference_standard_inert_warning(&self) -> Option<&'static str> {
        if self.custom_drag_table.is_some()
            && matches!(self.bc_reference_standard, BcReferenceStandard::ArmyStandardMetro)
        {
            Some(BC_REFERENCE_STANDARD_INERT_WARNING)
        } else {
            None
        }
    }

    /// Apply every input-conditioning step integration requires, exactly once (MBA-1415).
    ///
    /// **Every public entry point into integration must call this**, not just
    /// [`TrajectorySolver::new`]. `fast_trajectory::fast_integrate` and
    /// `fast_integrate_with_segments` are `pub` and are called directly by the Python
    /// binding, so conditioning that lives only in the solver constructor is invisible to
    /// those consumers: the caller sets a field, gets no error, and the value is silently
    /// ignored. That exact shape has already cost this project twice — MBA-1296 (a dropped
    /// field that zeroed Coriolis in production) and the MBA-1356 review catch
    /// (`cd_scale` dropped by the segmented fast path). Adding a step here rather than at a
    /// call site is what keeps the next one from repeating it.
    ///
    /// **Idempotent**, and it must stay that way: Monte Carlo builds a solver per sample, and
    /// the fast entry points may receive inputs that already passed through a solver. Each
    /// step below is either an absolute override (derived fields, powder-resolved velocity) or
    /// self-disarming (the BC conversion rewrites `bc_reference_standard` to the reference it
    /// just converted into, so a second call finds nothing to do). A step that scales a value
    /// in place without disarming itself would double-apply silently — do not add one.
    pub fn normalize_for_solve(&mut self) {
        // MBA-1365: normalize a declared Army-Standard-Metro BC reference to the ICAO
        // reference this engine's retardation constant (CD_TO_RETARD) is calibrated against,
        // before any retardation math runs. Every BC representation on these inputs (the
        // scalar, explicit Mach segments, and velocity segments) is a real BC and gets the
        // same treatment. `Icao` (the default) takes this branch never, so callers that never
        // set the field are byte-identical to pre-MBA-1365 behavior. A custom drag table
        // divides by sectional density rather than a BC, so the conversion is physically inert
        // there — see `bc_reference_standard_inert_warning`, which callers surface at their own
        // display sites (native stderr / WASM table text) from their OWN copy of the inputs,
        // before this runs.
        if matches!(
            self.bc_reference_standard,
            BcReferenceStandard::ArmyStandardMetro
        ) {
            self.bc_value *= crate::constants::ASM_TO_ICAO_BC;
            if let Some(segments) = self.bc_segments.as_mut() {
                for (_mach, bc) in segments.iter_mut() {
                    *bc *= crate::constants::ASM_TO_ICAO_BC;
                }
            }
            if let Some(segments) = self.bc_segments_data.as_mut() {
                for segment in segments.iter_mut() {
                    segment.bc_value *= crate::constants::ASM_TO_ICAO_BC;
                }
            }
            // Disarm: these BC values are now ICAO-referenced, which is simply true, and it
            // makes a second call a no-op instead of a silent second scaling.
            self.bc_reference_standard = BcReferenceStandard::Icao;
        }

        // Derived imperial fields, recomputed from the canonical SI values. These are an
        // absolute override, NOT a fallback: a caller that sets only `caliber_inches` /
        // `weight_grains` and leaves the SI fields at their defaults has those imperial values
        // overwritten here. SI is authoritative on these inputs.
        self.caliber_inches = self.bullet_diameter / 0.0254;
        self.weight_grains = self.bullet_mass / crate::constants::GRAINS_TO_KG;

        // Resolve the muzzle velocity for the ambient temperature before integration. A
        // measured powder-temperature -> velocity curve (data-driven, non-linear) takes
        // precedence when supplied; otherwise fall back to the linear powder-temperature-
        // sensitivity model (MBA-963). Both operate in canonical SI (Celsius, m/s), so every
        // solver built from these inputs — the main trajectory AND the zero-angle search —
        // sees the same temperature-resolved velocity. In particular, when a zero solve passes
        // the zero-day temperature, the curve automatically yields the zero-day velocity.
        // (The curve interpolates at the POWDER temperature — powder_curve_temp_c, falling
        // back to ambient. Air temperature still drives density separately; this only sets the
        // velocity. Absolute override, so applying it twice is a no-op.)
        self.muzzle_velocity = resolve_powder_adjusted_velocity(
            self.muzzle_velocity,
            self.temperature,
            self.use_powder_sensitivity,
            self.powder_temp_sensitivity,
            self.powder_temp,
            self.powder_temp_curve.as_deref(),
            self.powder_curve_temp_c,
        );
    }
}

impl Default for BallisticInputs {
    fn default() -> Self {
        let mass_kg = 0.01;
        let diameter_m = 0.00762;
        let bc = 0.5;
        let muzzle_angle_rad = 0.0;
        let bc_type = DragModel::G1;

        Self {
            // Core ballistics parameters
            bc_value: bc,
            bc_type,
            bc_reference_standard: BcReferenceStandard::Icao,
            bullet_mass: mass_kg,
            muzzle_velocity: 800.0,
            bullet_diameter: diameter_m,
            // MBA-1135: mass-based length estimate so the default is self-consistent with the
            // default mass/diameter (was a mass-blind 4.5-caliber literal). The twist default below
            // stays a fixed 1:12" per the ticket (a constant is a sensible velocity-agnostic default).
            bullet_length: crate::stability::estimate_bullet_length_m(diameter_m, mass_kg),

            // Targeting and positioning
            muzzle_angle: muzzle_angle_rad,
            target_distance: 100.0,
            azimuth_angle: 0.0,
            shot_azimuth: 0.0,
            shooting_angle: 0.0,
            cant_angle: 0.0,
            sight_height: 0.05,
            sight_offset_lateral_m: 0.0, // Sight directly above the bore (MBA-1396)
            muzzle_height: 0.0,       // Default 0 - height is in sight_height
            target_height: 0.0,       // Target at ground level by default
            zero_poi_vertical_m: 0.0, // No deliberate POI offset at the zero range (MBA-1359)
            zero_poi_horizontal_m: 0.0,
            ground_threshold: -100.0, // Effectively disable ground detection (allow bullet to drop 100m below start)

            // Environmental conditions
            altitude: 0.0,
            temperature: 15.0,
            pressure: 1013.25, // Standard sea level pressure (millibars)
            humidity: 0.5,     // 50% relative humidity
            latitude: None,

            // Wind conditions
            wind_speed: 0.0,
            wind_angle: 0.0,

            // Bullet characteristics
            twist_rate: 12.0, // 1:12" typical
            is_twist_right: true,
            caliber_inches: diameter_m / 0.0254, // Convert to inches
            weight_grains: mass_kg / crate::constants::GRAINS_TO_KG, // Convert to grains
            manufacturer: None,
            bullet_model: None,
            bullet_id: None,
            bullet_cluster: None,

            // Integration method selection
            use_rk4: true,           // Use Runge-Kutta methods by default
            use_adaptive_rk45: true, // Default to RK45 adaptive for best accuracy

            // Advanced effects (disabled by default)
            enable_advanced_effects: false,
            enable_magnus: false,
            enable_coriolis: false,
            use_powder_sensitivity: false,
            powder_temp_sensitivity: 0.0,
            powder_temp: 15.0,
            powder_temp_curve: None,
            powder_curve_temp_c: None,
            tipoff_yaw: 0.0,
            tipoff_decay_distance: 50.0,
            use_bc_segments: false,
            bc_segments: None,
            bc_segments_data: None,
            use_enhanced_spin_drift: false,
            use_form_factor: false,
            enable_wind_shear: false,
            wind_shear_model: "none".to_string(),
            enable_trajectory_sampling: false,
            sample_interval: 10.0, // Default 10 meter intervals
            drops_reference: DropsReference::Los, // historical LOS-perpendicular drops
            enable_pitch_damping: false,
            enable_precession_nutation: false,
            enable_aerodynamic_jump: false,
            use_cluster_bc: false, // Disabled by default for backward compatibility

            // Custom drag model support
            custom_drag_table: None,
            cd_scale: 1.0,

            // Legacy field for compatibility
            bc_type_str: None,
        }
    }
}

/// Interpolate a muzzle velocity (m/s) from a measured powder-temperature curve at
/// `temp_c` (Celsius). `curve` is `(temperature_celsius, velocity_m_s)` points; it is
/// sorted ascending by temperature before use. Values below the first point or above
/// the last are CLAMPED to the endpoint velocity (no extrapolation beyond measured
/// data), and segments are linearly interpolated. A single point yields a constant.
pub fn interpolate_powder_temp_curve(curve: &[(f64, f64)], temp_c: f64) -> f64 {
    debug_assert!(!curve.is_empty());
    if curve.is_empty() {
        return 0.0;
    }
    // Defensive: accept unsorted input by sorting a local copy only when needed.
    // Callers (CLI/WASM parsers) already sort, so the common path is a no-op scan.
    let mut sorted;
    let pts: &[(f64, f64)] = if curve.windows(2).all(|w| w[0].0 <= w[1].0) {
        curve
    } else {
        sorted = curve.to_vec();
        sorted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        &sorted
    };
    let n = pts.len();
    if temp_c <= pts[0].0 {
        return pts[0].1; // clamp below the coldest measured point
    }
    if temp_c >= pts[n - 1].0 {
        return pts[n - 1].1; // clamp above the hottest measured point
    }
    for i in 1..n {
        let (t0, v0) = pts[i - 1];
        let (t1, v1) = pts[i];
        if temp_c <= t1 {
            let span = t1 - t0;
            if span.abs() < f64::EPSILON {
                return v1; // coincident temps: avoid divide-by-zero, take the upper
            }
            let f = (temp_c - t0) / span;
            return v0 + f * (v1 - v0);
        }
    }
    pts[n - 1].1
}

/// Parse a `--sweep START:END:STEP` temperature range (display units) for the
/// `powder` command — ONE parser shared by the native CLI and WASM so their
/// validation cannot drift. Returns the
/// inclusive row temperatures. Guarded: STEP must be positive, END >= START, and the
/// row count is capped so a typo can't emit an unbounded table.
pub fn parse_powder_sweep(s: &str) -> Result<Vec<f64>, String> {
    const MAX_SWEEP_ROWS: usize = 500;
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() != 3 {
        return Err(format!(
            "Invalid --sweep '{}': expected START:END:STEP (e.g. \"20:110:10\")",
            s
        ));
    }
    let parse = |p: &str, name: &str| -> Result<f64, String> {
        p.trim()
            .parse::<f64>()
            .map_err(|_| format!("Invalid --sweep {}: '{}' is not a number", name, p.trim()))
    };
    let start = parse(parts[0], "START")?;
    let end = parse(parts[1], "END")?;
    let step = parse(parts[2], "STEP")?;
    if !step.is_finite() || step <= 0.0 {
        return Err(format!("Invalid --sweep STEP {}: must be positive", step));
    }
    if !start.is_finite() || !end.is_finite() || end < start {
        return Err(format!(
            "Invalid --sweep range {}:{}: END must be >= START",
            start, end
        ));
    }
    // Row count computed and bounds-checked in f64 BEFORE the usize cast: a huge
    // range would saturate the cast and overflow the `+ 1` (panic under
    // overflow-checks, wrap-to-zero in release — silently bypassing this cap).
    // The epsilon keeps a fractional STEP whose quotient lands a hair below an
    // integer (e.g. 0:0.3:0.1 -> 2.9999999999999996) from dropping the END row.
    let n_f = ((end - start) / step + 1e-9).floor();
    if !n_f.is_finite() || n_f + 1.0 > MAX_SWEEP_ROWS as f64 {
        return Err(format!(
            "--sweep would produce more than {} rows; use a larger STEP",
            MAX_SWEEP_ROWS
        ));
    }
    let n = n_f as usize + 1;
    // Index-multiplied (not accumulated) so float drift can't drop the END row.
    Ok((0..n).map(|i| start + step * i as f64).collect())
}

/// Resolve the powder-temperature-adjusted muzzle velocity (m/s) — the velocity the
/// solver actually flies. A non-empty measured `powder_temp_curve` takes precedence
/// (interpolated at `powder_curve_temp_c`, falling back to the ambient
/// `ambient_temperature_c`; clamped at the endpoints) and REPLACES the nominal
/// velocity; otherwise, when `use_powder_sensitivity` is set, the linear model adds
/// `sensitivity x (ambient - reference)` to it. All inputs canonical SI (Celsius,
/// m/s). This is the single shared implementation behind the trajectory solve, the
/// `powder` CLI/WASM command (MBA-737), and anything else that must agree with them.
pub fn resolve_powder_adjusted_velocity(
    nominal_velocity_mps: f64,
    ambient_temperature_c: f64,
    use_powder_sensitivity: bool,
    powder_temp_sensitivity_mps_per_c: f64,
    powder_reference_temp_c: f64,
    powder_temp_curve: Option<&[(f64, f64)]>,
    powder_curve_temp_c: Option<f64>,
) -> f64 {
    if let Some(curve) = powder_temp_curve {
        if !curve.is_empty() {
            let lookup_c = powder_curve_temp_c.unwrap_or(ambient_temperature_c);
            return interpolate_powder_temp_curve(curve, lookup_c);
        }
        // A supplied-but-empty curve has always suppressed the linear fallback
        // (the historical `else if`); preserve that exactly.
        return nominal_velocity_mps;
    }
    if use_powder_sensitivity {
        let temp_delta_c = ambient_temperature_c - powder_reference_temp_c;
        return nominal_velocity_mps + powder_temp_sensitivity_mps_per_c * temp_delta_c;
    }
    nominal_velocity_mps
}

// Wind conditions
#[derive(Debug, Clone)]
pub struct WindConditions {
    pub speed: f64, // m/s
    // radians, wind-FROM convention: 0 = headwind, PI/2 = from the right,
    // PI = tailwind, 3*PI/2 = from the left (matches WindSock / the bindings).
    pub direction: f64,
    /// Vertical wind component, m/s. POSITIVE = UPDRAFT (raises POI downrange); negative =
    /// downdraft. Default 0.0. Enters the wind vector via [`crate::wind::wind_vector`]'s third
    /// argument (MBA-728). Boundary-layer shear scales horizontal wind only — vertical passes
    /// through unscaled wherever shear is applied on top of this. This scalar field (like
    /// [`WindConditions::speed`]/[`WindConditions::direction`]) is ignored once downrange wind
    /// segments are set on the solver — each [`crate::wind::WindSegment`] carries its own
    /// `vertical_mps` instead.
    pub vertical_speed: f64,
}

impl Default for WindConditions {
    fn default() -> Self {
        Self {
            speed: 0.0,
            direction: 0.0,
            vertical_speed: 0.0,
        }
    }
}

// Atmospheric conditions
#[derive(Debug, Clone)]
pub struct AtmosphericConditions {
    pub temperature: f64, // Celsius
    pub pressure: f64,    // hPa
    /// Relative humidity as a PERCENT in `[0, 100]`. NOTE: [`BallisticInputs::humidity`]
    /// uses a 0–1 FRACTION instead — convert with `BallisticInputs::humidity_percent` when
    /// crossing between them (MBA-722).
    pub humidity: f64,
    pub altitude: f64, // meters
}

impl Default for AtmosphericConditions {
    fn default() -> Self {
        Self {
            temperature: 15.0,
            pressure: 1013.25,
            humidity: 50.0,
            altitude: 0.0,
        }
    }
}

// Trajectory point data
#[derive(Debug, Clone)]
pub struct TrajectoryPoint {
    pub time: f64,
    pub position: Vector3<f64>,
    pub velocity_magnitude: f64,
    pub kinetic_energy: f64,
    /// The projectile's own drag coefficient at this point (MBA-1423).
    ///
    /// Filled in one pass after integration finishes, not at the sites that build this struct,
    /// so a solver family added later cannot silently omit it. `None` when sectional density is
    /// unavailable — see [`TrajectorySolver::effective_drag_coefficient`] for the definition and
    /// why this is the projectile's Cd rather than the reference table's.
    pub drag_coefficient: Option<f64>,
}

impl TrajectoryPoint {
    /// The value `--with-drag-coefficient` emits for this point, or `None` when the JSON key
    /// must be ABSENT — flag off, or sectional density unknown (MBA-1427).
    ///
    /// Host-compilable on purpose. The browser terminal's emit site lives in `src/wasm.rs`,
    /// which is `cfg(target_arch = "wasm32")`-gated out of every native build — CI never
    /// executes a line of it. Keeping the gating rule here means a native test can pin the
    /// absent-not-null contract even though the surface that applies it only compiles for wasm.
    pub fn drag_coefficient_json_value(&self, with_drag_coefficient: bool) -> Option<f64> {
        if with_drag_coefficient {
            self.drag_coefficient
        } else {
            None
        }
    }
}

// Trajectory result
#[derive(Debug, Clone)]
pub struct TrajectoryResult {
    pub max_range: f64,
    pub max_height: f64,
    pub time_of_flight: f64,
    pub impact_velocity: f64,
    pub impact_energy: f64,
    /// Projectile mass used to derive full-state observation energy.
    pub projectile_mass_kg: f64,
    /// Height of the horizontal line of sight in the solver's ground-referenced frame.
    pub line_of_sight_height_m: f64,
    /// Station speed of sound used for Mach observations and transition flags.
    pub station_speed_of_sound_mps: f64,
    /// Explicit reason the integration stopped; consumers must not infer this from the endpoint.
    pub termination: TrajectoryTermination,
    pub points: Vec<TrajectoryPoint>,
    pub sampled_points: Option<Vec<TrajectorySample>>, // Trajectory samples at regular intervals
    pub min_pitch_damping: Option<f64>, // Minimum pitch damping coefficient (for stability warning)
    pub transonic_mach: Option<f64>,    // Mach number when entering transonic regime
    pub angular_state: Option<AngularState>, // Final angular state if precession/nutation enabled
    pub max_yaw_angle: Option<f64>,     // Maximum yaw angle during flight (radians)
    pub max_precession_angle: Option<f64>, // Maximum precession angle (radians)
    // MBA-959: aerodynamic-jump components applied at the muzzle (None unless
    // enable_aerodynamic_jump). EXPERIMENTAL.
    pub aerodynamic_jump: Option<crate::aerodynamic_jump::AerodynamicJumpComponents>,
    /// Downrange distance (m) of the downward Mach 1.2 crossing (station speed of sound),
    /// populated identically by all three solver paths (Euler/RK4/RK45). `None` if the
    /// trajectory never crosses 1.2 while descending (e.g. launched already below 1.2, or the
    /// solve terminates while still above it). MBA-1405: feeds `mv_calibration_window`.
    pub mach_1_2_distance_m: Option<f64>,
    /// Downrange distance (m) of the downward Mach 1.0 crossing. `None` if the trajectory
    /// never goes subsonic within the solve. MBA-1405: the far edge of the MV calibration
    /// window (`mv_calibration_window`) is this distance.
    pub mach_1_0_distance_m: Option<f64>,
    /// Downrange distance (m) of the downward Mach 0.9 crossing. `None` if the trajectory
    /// never crosses 0.9 while descending. MBA-1405: feeds `dsf_window_start`. NOTE: unlike
    /// the 1.2/1.0 crossings, this one is intentionally NOT reflected in the historical flat
    /// `transonic_distances` Vec threaded through
    /// [`crate::trajectory_sampling::TrajectoryData`] — existing consumers of that Vec
    /// interpret it strictly as `{1.2, 1.0}`, so this field is the only place the 0.9
    /// crossing is recorded.
    pub mach_0_9_distance_m: Option<f64>,
}

const RK45_TOLERANCE: f64 = 1e-6;
const RK45_SAFETY_FACTOR: f64 = 0.9;
const RK45_MAX_DT: f64 = 0.01;
const RK45_MIN_DT: f64 = 1e-6;
const TRAJECTORY_TIME_LIMIT_S: f64 = 100.0;

/// Hard ceiling for points retained by a single [`TrajectorySolver`] result.
///
/// The cap applies across Euler, fixed RK4, and adaptive RK45, including the exact terminal
/// endpoint. Solves that would exceed it return [`BallisticsError`] instead of
/// truncating or growing their point buffer without bound.
pub const MAX_TRAJECTORY_POINTS: usize = 250_000;

/// Pack the CLI solver's split position/velocity vectors into the shared six-component RK45 norm.
fn cli_rk45_error_norm(
    position: &Vector3<f64>,
    velocity: &Vector3<f64>,
    fifth_position: &Vector3<f64>,
    fifth_velocity: &Vector3<f64>,
    fourth_position: &Vector3<f64>,
    fourth_velocity: &Vector3<f64>,
) -> f64 {
    let pack_state = |position: &Vector3<f64>, velocity: &Vector3<f64>| {
        Vector6::new(
            position.x, position.y, position.z, velocity.x, velocity.y, velocity.z,
        )
    };
    let state = pack_state(position, velocity);
    let fifth_order = pack_state(fifth_position, fifth_velocity);
    let fourth_order = pack_state(fourth_position, fourth_velocity);

    crate::trajectory_integration::rk45_error_norm(&state, &fifth_order, &fourth_order)
}

struct Rk45Trial {
    position: Vector3<f64>,
    velocity: Vector3<f64>,
    suggested_dt: f64,
    error: f64,
}

struct Rk45AcceptedStep {
    position: Vector3<f64>,
    velocity: Vector3<f64>,
    used_dt: f64,
    next_dt: f64,
    error: f64,
}

/// Tracks the downward (decelerating) Mach crossings of a solved trajectory.
///
/// `record_downward_crossings` is called once per integration step from all three solver
/// paths (Euler/RK4/RK45). It both (a) appends to the historical flat `distances` Vec — used
/// unchanged by [`crate::trajectory_sampling`] to flag sampled points — and (b) records each
/// threshold's exact crossing distance on `self`, for callers that want the labeled value
/// directly (MBA-1405) rather than positionally decoding the flat Vec.
///
/// IMPORTANT: only the historical 1.2/1.0 thresholds are appended to `distances`. The 0.9
/// threshold added for MBA-1405 is recorded ONLY on `self.mach_0_9_distance_m` — existing
/// consumers of the flat Vec interpret it strictly as `{mach_1_2, mach_1_0}` (in that order,
/// when present), so appending a third entry would silently corrupt their positional
/// assumptions.
#[derive(Default)]
struct MachTransitionTracker {
    previous_mach: Option<f64>,
    crossed_transonic: bool,
    crossed_subsonic: bool,
    crossed_narrow: bool,
    /// Downrange distance (m) of the downward Mach 1.2 crossing, once seen. Mirrors the first
    /// flat-Vec entry.
    mach_1_2_distance_m: Option<f64>,
    /// Downrange distance (m) of the downward Mach 1.0 crossing, once seen. Mirrors the second
    /// flat-Vec entry.
    mach_1_0_distance_m: Option<f64>,
    /// Downrange distance (m) of the downward Mach 0.9 crossing, once seen. NOT reflected in
    /// the flat Vec (see struct docs).
    mach_0_9_distance_m: Option<f64>,
}

impl MachTransitionTracker {
    fn record_downward_crossings(&mut self, mach: f64, downrange_m: f64, distances: &mut Vec<f64>) {
        if !mach.is_finite() {
            self.previous_mach = None;
            return;
        }

        if let Some(previous_mach) = self.previous_mach {
            if !self.crossed_transonic && previous_mach >= 1.2 && mach < 1.2 {
                self.crossed_transonic = true;
                distances.push(downrange_m);
                self.mach_1_2_distance_m = Some(downrange_m);
            }
            if !self.crossed_subsonic && previous_mach >= 1.0 && mach < 1.0 {
                self.crossed_subsonic = true;
                distances.push(downrange_m);
                self.mach_1_0_distance_m = Some(downrange_m);
            }
            if !self.crossed_narrow && previous_mach >= 0.9 && mach < 0.9 {
                self.crossed_narrow = true;
                // Intentionally NOT pushed to `distances` — see struct docs.
                self.mach_0_9_distance_m = Some(downrange_m);
            }
        }
        self.previous_mach = Some(mach);
    }
}

impl TrajectoryResult {
    /// Interpolate position at a given downrange distance (X coordinate, McCoy).
    /// Returns the interpolated (x, y, z) position at that range.
    /// If the target range exceeds the trajectory, returns the last point.
    pub fn position_at_range(&self, target_range: f64) -> Option<Vector3<f64>> {
        if self.points.is_empty() {
            return None;
        }

        // Find the two points that bracket the target range
        for i in 0..self.points.len() - 1 {
            let p1 = &self.points[i];
            let p2 = &self.points[i + 1];

            // Check if target range is between these two points (X is downrange)
            if p1.position.x <= target_range && p2.position.x >= target_range {
                // Linear interpolation factor
                let dx = p2.position.x - p1.position.x;
                if dx.abs() < 1e-10 {
                    return Some(p1.position);
                }
                let t = (target_range - p1.position.x) / dx;

                // Interpolate Y and Z, use exact target_range for X
                return Some(Vector3::new(
                    target_range,
                    p1.position.y + t * (p2.position.y - p1.position.y),
                    p1.position.z + t * (p2.position.z - p1.position.z),
                ));
            }
        }

        // Target range is beyond trajectory - return last point
        self.points.last().map(|p| p.position)
    }
}

// Trajectory solver
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum StationAtmosphereResolution {
    /// Preserve the historical CLI/FFI convention: sea-level standard values at a nonzero
    /// altitude are treated as omitted and resolved from the ICAO atmosphere.
    LegacyDefaultSentinels,
    /// Temperature and pressure have already been resolved by a presence-aware caller and must
    /// remain authoritative even when they equal the historical sentinel values.
    Authoritative,
}

#[derive(Clone)]
pub struct TrajectorySolver {
    inputs: BallisticInputs,
    wind: WindConditions,
    atmosphere: AtmosphericConditions,
    station_atmosphere_resolution: StationAtmosphereResolution,
    max_range: f64,
    time_step: f64,
    max_trajectory_points: usize,
    cluster_bc: Option<ClusterBCDegradation>,
    /// Geometry-derived `(longitudinal, transverse)` moments used by angular diagnostics.
    precession_nutation_inertias: (f64, f64),
    /// Optional downrange-segmented wind. When `Some`, the per-step wind vector is
    /// looked up by downrange distance from this `WindSock` and the scalar `wind`
    /// field is ignored. When `None`, the constant `wind` vector is used (default),
    /// so a non-segmented solve is numerically identical to pre-feature behavior.
    wind_sock: Option<crate::wind::WindSock>,
    /// Optional downrange-segmented atmosphere (MBA-1137). When `Some`, the per-substep local
    /// atmosphere recompute samples the base (station-referenced) temperature/pressure/humidity by
    /// downrange distance from this `AtmoSock`, then feeds them through the SAME altitude-lapse
    /// pipeline as a single-station solve — so the downrange zone and the vertical altitude lapse
    /// compose without double-counting. When `None` (default), the resolved single-station
    /// conditions are used.
    atmo_sock: Option<crate::atmosphere::AtmoSock>,
}

/// Which frame a zero solve's `target_height_m` lives in (MBA-1412).
///
/// `SightLine`: the classic zero contract — the height is sight geometry (typically the sight
/// height) and the rifle zeroes LEVEL, per the MBA-1286 doctrine that a zero is a property of a
/// level rifle; `shooting_angle` applies to the subsequent shot only. Used by
/// `calculate_zero_angle_with_conditions` and every legacy surface behind it (bindings, WASM,
/// FFI, CLI).
///
/// `WorldVertical`: solve-json v1's documented contract (docs/SOLVE_JSON_V1.md) — the height is
/// an absolute world-vertical height above the ground datum; inclined zeroing projects the
/// shot-frame trajectory into the world frame (MBA-1302 behavior, unchanged).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ZeroTargetFrame {
    SightLine,
    WorldVertical,
}

/// Both line-of-sight crossings a single fixed bore angle can produce (Tier 2 whole-branch
/// review, C2): a rifle whose sight sits above its bore generally crosses a level line of
/// sight TWICE on a rising shot — an ascending "near zero" close to the muzzle, and a
/// descending "far zero" past the apex. The classic 25/300-yard battle-zero pairing is one
/// bore angle producing both. `TrajectorySolver`'s internal `find_zero_range` and the public
/// `calculate_zero_range_from_angle_*` functions return both rather than silently choosing
/// one, because `find_zero_angle`'s own choice of root depends on the target distance it was
/// asked to solve for — the forward and inverse solvers are NOT a single-valued pair.
///
/// Either field may be `None` when the solved envelope doesn't reach that crossing (e.g. a
/// short `--max-range` that stops before the far crossing, or a target height high enough
/// that the trajectory never rises to it at all). At least one is always `Some` on `Ok`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ZeroCrossings {
    /// Downrange distance (meters) of the near, ascending crossing, if the solved envelope
    /// reaches it.
    pub near_m: Option<f64>,
    /// Downrange distance (meters) of the far, descending crossing -- what "sighted in at D
    /// yards" conventionally means -- if the solved envelope reaches it.
    pub far_m: Option<f64>,
}

impl TrajectorySolver {
    pub fn new(
        inputs: BallisticInputs,
        wind: WindConditions,
        atmosphere: AtmosphericConditions,
    ) -> Self {
        Self::new_with_station_atmosphere_resolution(
            inputs,
            wind,
            atmosphere,
            StationAtmosphereResolution::LegacyDefaultSentinels,
        )
    }

    /// Construct a solver from station temperature and pressure that a presence-aware service
    /// has already resolved. Unlike [`Self::new`], exact sea-level standard values remain
    /// authoritative at nonzero altitude rather than acting as legacy omission sentinels.
    ///
    /// `pub` (MBA-1397; was `pub(crate)`): a QNH-aware caller (native CLI, WASM) that has
    /// already reduced a declared altimeter setting to station pressure via
    /// [`crate::atmosphere::resolve_station_conditions_with_pressure_mode`] needs this
    /// constructor too, and those callers live in a different crate (`src/main.rs` is a
    /// separate binary crate over the library). Widening visibility only, no signature
    /// change: every existing `pub(crate)` caller is unaffected.
    pub fn new_with_resolved_station_atmosphere(
        inputs: BallisticInputs,
        wind: WindConditions,
        atmosphere: AtmosphericConditions,
    ) -> Self {
        Self::new_with_station_atmosphere_resolution(
            inputs,
            wind,
            atmosphere,
            StationAtmosphereResolution::Authoritative,
        )
    }

    fn new_with_station_atmosphere_resolution(
        mut inputs: BallisticInputs,
        wind: WindConditions,
        atmosphere: AtmosphericConditions,
        station_atmosphere_resolution: StationAtmosphereResolution,
    ) -> Self {
        // MBA-1415: every input-conditioning step lives in one shared, idempotent helper that
        // the fast_integrate* entry points call too. Keeping it out of this constructor is the
        // point: conditioning that lives only here is invisible to the external bindings that
        // call the fast path directly.
        inputs.normalize_for_solve();

        // Initialize cluster BC if enabled
        let cluster_bc = if inputs.use_cluster_bc {
            Some(ClusterBCDegradation::new())
        } else {
            None
        };
        let precession_nutation_inertias = projectile_moments_of_inertia(
            inputs.bullet_mass,
            inputs.bullet_diameter,
            inputs.bullet_length,
        );

        Self {
            inputs,
            wind,
            atmosphere,
            station_atmosphere_resolution,
            max_range: 1000.0,
            time_step: 0.001,
            max_trajectory_points: MAX_TRAJECTORY_POINTS,
            cluster_bc,
            precession_nutation_inertias,
            wind_sock: None,
            atmo_sock: None,
        }
    }

    pub fn set_max_range(&mut self, range: f64) {
        self.max_range = range;
    }

    pub fn set_time_step(&mut self, step: f64) {
        self.time_step = step;
    }

    /// Calculate a level-rifle zero with this solver's configured atmosphere, wind (including
    /// downrange segments), effects, integration method, and time step, then install the resulting
    /// muzzle angle on this solver. The solver is mutated only after a zero has converged.
    pub(crate) fn calculate_and_set_zero_angle(
        &mut self,
        target_distance_m: f64,
        target_height_m: f64,
        frame: ZeroTargetFrame,
    ) -> Result<f64, BallisticsError> {
        let angle = self.find_zero_angle(target_distance_m, target_height_m, frame)?;
        // MBA-1359: a deliberate POI offset at the zero range (Kestrel ZH/ZO) is an angular
        // bias ON the solved angle, applied post-solve so the zero trials above still solve
        // perfect convergence. The small-angle addition is frame-independent (identical in
        // SightLine and WorldVertical), and with both offsets at their 0.0 defaults the
        // additions below are exact no-ops (byte-identical). The returned angle carries the
        // vertical bias so auto-zero callers that copy it onto flight inputs inherit it; the
        // horizontal bias lands on THIS solver's azimuth_angle — separate flight inputs must
        // add `windage_zero_bias_rad` themselves.
        let angle = if target_distance_m > 0.0 {
            angle + self.inputs.zero_poi_vertical_m / target_distance_m
        } else {
            angle
        };
        self.inputs.muzzle_angle = angle;
        self.apply_windage_zero_bias(target_distance_m);
        Ok(angle)
    }

    /// Apply just the windage-zero convergence bias (see
    /// [`BallisticInputs::windage_zero_bias_rad`]) to this solver's azimuth, for a zero at
    /// `target_distance_m`, without searching for or installing an elevation.
    ///
    /// Used by [`Self::calculate_and_set_zero_angle`] above as part of a full zero search, and
    /// by solve_v1's round-trip re-solve path (0.33.0 decision-support Task 2) on its own: when
    /// an explicit `muzzle_angle_rad` already carries the correct elevation (round-tripped from
    /// a previous `resolved_request`, which bakes any `zero_poi_vertical_m` bias into that
    /// angle already), the elevation search does not run, but the windage term is a SEPARATE
    /// quantity that depends only on `zero_distance_m` -- not on whether the elevation was
    /// searched for or supplied directly -- so it must still be applied here or an
    /// offset-mounted sight / deliberate horizontal zero bias would silently stop converging
    /// once a request round-trips.
    pub(crate) fn apply_windage_zero_bias(&mut self, target_distance_m: f64) {
        self.inputs.azimuth_angle += self.inputs.windage_zero_bias_rad(target_distance_m);
    }

    fn find_zero_angle(
        &self,
        target_distance_m: f64,
        target_height_m: f64,
        frame: ZeroTargetFrame,
    ) -> Result<f64, BallisticsError> {
        // Binary search for the angle that hits the target. Use only positive angles to ensure a
        // proper upward ballistic arc.
        let mut low_angle = 0.0;
        let mut high_angle = 0.2; // about 11 degrees
        let tolerance = 1e-7;
        let max_iterations = 60;

        // MBA-194: validate the initial bracket before starting the binary search.
        let low_height = self.zero_trial_height_at(low_angle, target_distance_m, frame)?;
        let high_height = self.zero_trial_height_at(high_angle, target_distance_m, frame)?;

        match (low_height, high_height) {
            (Some(low_height), Some(high_height)) => {
                let low_error = low_height - target_height_m;
                let high_error = high_height - target_height_m;

                if low_error > 0.0 && high_error > 0.0 {
                    // Both angles overshoot. Zero degrees is the lowest supported launch angle;
                    // retain the historical behavior and let the search choose its best result.
                } else if low_error < 0.0 && high_error < 0.0 {
                    // Both angles undershoot. Preserve the historical expansion up to 45 degrees.
                    let mut expanded = false;
                    for multiplier in [2.0, 3.0, 4.0] {
                        let new_high = (high_angle * multiplier).min(0.785);
                        if let Ok(Some(height)) =
                            self.zero_trial_height_at(new_high, target_distance_m, frame)
                        {
                            if height - target_height_m > 0.0 {
                                high_angle = new_high;
                                expanded = true;
                                break;
                            }
                        }
                        if new_high >= 0.785 {
                            break;
                        }
                    }
                    if !expanded {
                        return Err("Cannot find zero angle: target beyond effective range even at maximum angle".into());
                    }
                }
            }
            (None, Some(_)) => {
                // The low angle does not reach the target while the high angle does; the search
                // will raise the low end until it reaches a valid trajectory.
            }
            (Some(_), None) => {
                return Err(
                    "Cannot find zero angle: high angle trajectory doesn't reach target distance"
                        .into(),
                );
            }
            (None, None) => {
                return Err(
                    "Cannot find zero angle: trajectory cannot reach target distance at any angle"
                        .into(),
                );
            }
        }

        for _ in 0..max_iterations {
            let mid_angle = (low_angle + high_angle) / 2.0;
            match self.zero_trial_height_at(mid_angle, target_distance_m, frame)? {
                Some(height) => {
                    let error = height - target_height_m;
                    // MBA-193: height accuracy is the primary convergence criterion. At 0.1 mm,
                    // short-range zero-day atmosphere differences remain observable.
                    if error.abs() < 0.0001 {
                        return Ok(mid_angle);
                    }

                    // Only use angle tolerance after precision is exhausted and the remaining
                    // height error is still practically acceptable.
                    if (high_angle - low_angle).abs() < tolerance {
                        if error.abs() < 0.01 {
                            return Ok(mid_angle);
                        }
                        return Err("Zero angle did not converge: residual height error too large (target not reachable / not bracketed)".into());
                    }

                    if error > 0.0 {
                        high_angle = mid_angle;
                    } else {
                        low_angle = mid_angle;
                    }
                }
                None => {
                    low_angle = mid_angle;
                    if (high_angle - low_angle).abs() < tolerance {
                        return Err("Trajectory cannot reach target distance - angle converged without valid solution".into());
                    }
                }
            }
        }

        Err("Failed to find zero angle".into())
    }

    /// Solve one zero-angle trial without losing any solver configuration. Only the trial clone's
    /// launch angle, level-rifle convention, and integration range differ from the final solve.
    fn zero_trial_height_at(
        &self,
        angle_rad: f64,
        target_distance_m: f64,
        frame: ZeroTargetFrame,
    ) -> Result<Option<f64>, BallisticsError> {
        let mut trial = self.clone();
        trial.inputs.muzzle_angle = angle_rad;
        // MBA-959: zero on the bare bore so aerodynamic jump remains an additive fire-time POI
        // shift rather than being silently absorbed by the zero search.
        trial.inputs.enable_aerodynamic_jump = false;
        // MBA-1286: a zero is a property of a level rifle's sight geometry. Cant is applied only
        // to the subsequent shot.
        trial.inputs.cant_angle = 0.0;
        // MBA-1412: in the SightLine contract the incline is likewise a shot-time condition,
        // not sight geometry. MBA-1302 left it in the trial while shot_frame_altitude lifted
        // heights into the gravity frame, adding d*sin(shooting_angle) (~9 m at 5.71 deg /
        // 100 yd) against a sight-frame target height, so no inclined legacy caller could
        // bracket. Zero level; shoot inclined. solve-json v1's WorldVertical contract keeps
        // the MBA-1302 projection (its callers supply world-frame heights).
        if frame == ZeroTargetFrame::SightLine {
            trial.inputs.shooting_angle = 0.0;
        }
        trial.set_max_range(target_distance_m * 2.0);
        let result = trial.solve()?;

        for (index, point) in result.points.iter().enumerate() {
            if point.position.x >= target_distance_m {
                let shot_y_m = if index == 0 {
                    point.position.y
                } else {
                    let previous = &result.points[index - 1];
                    let span = point.position.x - previous.position.x;
                    let fraction = (target_distance_m - previous.position.x) / span;
                    previous.position.y + fraction * (point.position.y - previous.position.y)
                };
                return Ok(Some(crate::atmosphere::shot_frame_altitude(
                    0.0,
                    target_distance_m,
                    shot_y_m,
                    trial.inputs.shooting_angle,
                )));
            }
        }
        Ok(None)
    }

    /// Inverse of [`Self::find_zero_angle`] (MBA-1402), corrected by the Tier 2 whole-branch
    /// review (C2): given a FIXED bore angle rather than searching for one, run that
    /// trajectory and locate the downrange distance(s) where it crosses `target_height_m` —
    /// the zero RANGE(S) a stored/portable bore angle implies under THESE conditions. Hornady
    /// and Kestrel 4DOF treat the zero angle as the portable quantity (captured once, then
    /// reusable independent of the day it was solved); this is what makes that value usable
    /// again.
    ///
    /// A rifle whose sight sits above the bore generally crosses a level line of sight TWICE
    /// on a rising shot: once ascending, close to the muzzle (the "near zero"), and once
    /// descending past the apex (the "far zero"). Both are equally real zero ranges for the
    /// SAME bore angle — this is the classic 25/300-yard battle-zero relationship: a rifle
    /// zeroed at 25 yd (the near crossing of its bore angle) is, for the same angle, ALSO
    /// zeroed again around 300 yd (the far crossing). [`Self::find_zero_angle`]'s own choice
    /// of root depends on which target distance it was asked to solve for, so this is not a
    /// single-valued inverse of that function — it is this function's job to report both
    /// roots it can find rather than silently pick one and call it "the" zero range.
    ///
    /// Returns [`ZeroCrossings`] with the near (first ascending, error negative → positive)
    /// and far (descending, error positive → negative) crossings independently, either of
    /// which may be absent if the solved envelope doesn't reach it (e.g. `--max-range` too
    /// short to reach the far crossing, or an angle so shallow / a target height so high that
    /// the near crossing never happens). Errors only when NEITHER crossing is found. Each
    /// crossing distance is linearly interpolated between its two bracketing sample points,
    /// mirroring `zero_trial_height_at`'s own height-at-a-known-distance interpolation in the
    /// other direction.
    fn find_zero_range(
        &self,
        angle_rad: f64,
        target_height_m: f64,
        frame: ZeroTargetFrame,
    ) -> Result<ZeroCrossings, BallisticsError> {
        let mut trial = self.clone();
        trial.inputs.muzzle_angle = angle_rad;
        // Mirrors zero_trial_height_at's trial setup exactly (MBA-959 / MBA-1286 / MBA-1412),
        // so both directions apply the same conventions.
        trial.inputs.enable_aerodynamic_jump = false;
        trial.inputs.cant_angle = 0.0;
        if frame == ZeroTargetFrame::SightLine {
            trial.inputs.shooting_angle = 0.0;
        }
        let result = trial.solve()?;

        // The height-above-target-height error is, for a normal flat-fire trajectory,
        // unimodal: it starts below the line of sight (bore below sight), rises through it
        // ONCE (ascending -- the near zero), continues to the apex, then falls back through
        // it ONCE more (descending -- the far zero) before diverging away for good. So at
        // most one ascending and one descending crossing exist; classify each sign change by
        // direction rather than by position (first/last) so a near-only trajectory (solved
        // range too short to reach the far crossing) is never mislabeled as "far".
        let mut near_crossing: Option<f64> = None;
        let mut far_crossing: Option<f64> = None;
        let mut previous: Option<(f64, f64)> = None; // (downrange_m, height_error_m)
        for point in &result.points {
            let height = crate::atmosphere::shot_frame_altitude(
                0.0,
                point.position.x,
                point.position.y,
                trial.inputs.shooting_angle,
            );
            let error = height - target_height_m;
            if let Some((prev_x, prev_error)) = previous {
                if prev_error == 0.0 {
                    // An exact touch is itself a crossing; assign it to whichever slot is
                    // still open (near first, since it can only occur before far).
                    if near_crossing.is_none() {
                        near_crossing = Some(prev_x);
                    } else {
                        far_crossing = Some(prev_x);
                    }
                }
                if prev_error * error < 0.0 {
                    let fraction = prev_error / (prev_error - error);
                    let crossing = prev_x + fraction * (point.position.x - prev_x);
                    if prev_error < 0.0 && error > 0.0 {
                        // Ascending: the near zero. Keep only the first one found.
                        if near_crossing.is_none() {
                            near_crossing = Some(crossing);
                        }
                    } else {
                        // Descending: the far zero.
                        far_crossing = Some(crossing);
                    }
                }
            }
            previous = Some((point.position.x, error));
        }
        // The loop above only resolves a point once a later point makes it `previous`; catch
        // the trajectory's very last sampled point landing exactly on the line too.
        if let Some((last_x, last_error)) = previous {
            if last_error == 0.0 {
                if near_crossing.is_none() {
                    near_crossing = Some(last_x);
                } else {
                    far_crossing = Some(last_x);
                }
            }
        }

        if near_crossing.is_none() && far_crossing.is_none() {
            return Err(BallisticsError::from(
                "Cannot find zero range: this angle never crosses the target height within the \
                 solved range (angle too shallow to reach it, or both crossings lie beyond \
                 the solver's max range)."
                    .to_string(),
            ));
        }

        Ok(ZeroCrossings {
            near_m: near_crossing,
            far_m: far_crossing,
        })
    }

    /// Equivalent horizontal range for an inclined shot (MBA-1395): the flat-fire range
    /// whose ANGULAR elevation correction — relative to the SAME dialed zero — equals the
    /// inclined solution's angular correction at `target_range_m`. This is the
    /// "shoot-to" range SIG BDX (AMR), Leica (EHR), and Gunwerks BR2 report so users of
    /// fixed BDC turrets/reticles can dial as if the shot were flat, computed by an
    /// inverse lookup over one flat re-solve (McDonald's Sierra inclined-fire treatment /
    /// Litz) rather than the rifleman's-rule cosine approximation, which matches LINEAR
    /// drop and accumulates error at long range.
    ///
    /// The solver's current state IS the inclined solution: `muzzle_angle` carries the
    /// solved zero and `shooting_angle` the look angle. The flat reference clones the
    /// solver and zeroes only `shooting_angle` (the WASM auto-zero flat-clone precedent),
    /// keeping the zero, wind, atmosphere, and integration configuration identical, so
    /// both angular corrections are measured against the same zero.
    ///
    /// Returns `None` where the inverse is ill-defined instead of guessing:
    /// - `target_range_m <= zero_distance_m` (inside the zero, no BDC correction exists);
    /// - the inclined correction at the target is not positive (bullet at/above the LOS);
    /// - either solve fails, or the flat solve cannot bracket the correction.
    ///
    /// Past the zero range the flat correction is monotone in range, so a plain bisection
    /// over the flat solve's interpolated corrections converges; the answer is refined to
    /// centimeter level, far below any display precision.
    pub fn equivalent_horizontal_range(
        &self,
        target_range_m: f64,
        zero_distance_m: f64,
    ) -> Option<f64> {
        if !target_range_m.is_finite() || !zero_distance_m.is_finite() {
            return None;
        }
        if target_range_m <= zero_distance_m || target_range_m <= 0.0 {
            return None;
        }

        // Height of the bullet path at a downrange distance, linearly interpolated
        // between integration points (mirrors zero_trial_height_at's interpolation).
        // In the shot frame the LOS is the x-axis, so drop below the LOS is
        // `los_height - y` for the inclined and flat solves alike.
        fn path_y_at(points: &[TrajectoryPoint], distance_m: f64) -> Option<f64> {
            let mut previous: Option<(f64, f64)> = None;
            for point in points {
                if point.position.x >= distance_m {
                    return Some(match previous {
                        None => point.position.y,
                        Some((prev_x, prev_y)) => {
                            let span = point.position.x - prev_x;
                            if span <= 0.0 {
                                point.position.y
                            } else {
                                let fraction = (distance_m - prev_x) / span;
                                prev_y + fraction * (point.position.y - prev_y)
                            }
                        }
                    });
                }
                previous = Some((point.position.x, point.position.y));
            }
            None
        }

        // Inclined solution's angular correction at the target (radians, small-angle:
        // drop below the LOS divided by range). Sampling is display-side only — turn it
        // off in both internal solves.
        let mut inclined = self.clone();
        inclined.inputs.enable_trajectory_sampling = false;
        let inclined_result = inclined.solve().ok()?;
        let los_height = inclined_result.line_of_sight_height_m;
        let inclined_drop = los_height - path_y_at(&inclined_result.points, target_range_m)?;
        let correction = inclined_drop / target_range_m;
        if correction <= 0.0 {
            return None;
        }

        // One flat re-solve against the same zero.
        let mut flat = self.clone();
        flat.inputs.enable_trajectory_sampling = false;
        flat.inputs.shooting_angle = 0.0;
        let flat_result = flat.solve().ok()?;
        let flat_correction_at = |range_m: f64| -> Option<f64> {
            Some((los_height - path_y_at(&flat_result.points, range_m)?) / range_m)
        };

        // Bracket, then bisect. Gravity's along-LOS component only shrinks under a look
        // angle, so the flat correction at the true range bounds the inclined one from
        // above and the equivalent range lies in (zero_distance, target_range]. The flat
        // clone can strike the ground SHORT of the inclined terminal (an uphill shot
        // stays airborne longer), so the upper bracket is additionally capped at the
        // flat solve's own terminal distance — the shoot-to range is shorter than the
        // true range, so the cap does not exclude the root.
        let flat_terminal_m = flat_result.points.last().map(|p| p.position.x)?;
        let mut low = zero_distance_m.max(1.0);
        let mut high = target_range_m.min(flat_terminal_m);
        if high <= low {
            return None;
        }
        if flat_correction_at(low)? - correction > 0.0 {
            return None; // no bracket below
        }
        if flat_correction_at(high)? - correction < 0.0 {
            return None; // no bracket above (correction not reachable flat)
        }
        for _ in 0..60 {
            let mid = 0.5 * (low + high);
            let error = flat_correction_at(mid)? - correction;
            if error.abs() == 0.0 {
                return Some(mid);
            }
            if error < 0.0 {
                low = mid;
            } else {
                high = mid;
            }
            if high - low < 0.01 {
                break;
            }
        }
        Some(0.5 * (low + high))
    }

    /// Reject malformed state before it reaches an integration loop.
    ///
    /// `new` resolves powder-temperature velocity overrides and refreshes the imperial mirror
    /// fields, so validation belongs here: it sees the effective muzzle velocity, covers values
    /// changed through solver setters, and applies uniformly to Euler, RK4, and RK45.
    fn validate_for_solve(&self) -> Result<(), BallisticsError> {
        let require_finite = |name: &str, value: f64| {
            if value.is_finite() {
                Ok(())
            } else {
                Err(BallisticsError::from(format!("{name} must be finite")))
            }
        };
        let require_positive = |name: &str, value: f64| {
            if value.is_finite() && value > 0.0 {
                Ok(())
            } else {
                Err(BallisticsError::from(format!(
                    "{name} must be finite and greater than zero"
                )))
            }
        };

        // These four quantities are required by every point-mass solve. In particular, validate
        // muzzle_velocity after `new` has applied a measured curve or linear powder correction.
        // A custom drag table supplies the actual Cd and divides by sectional density, so bc_value
        // is physically ignored (see custom_drag_denominator). Require it only in the no-table case;
        // mass + diameter are always required (they drive the SD denominator when a table is set).
        if self.inputs.custom_drag_table.is_none() {
            require_positive("bc_value", self.inputs.bc_value)?;
        }
        require_positive("bullet_mass", self.inputs.bullet_mass)?;
        require_positive("bullet_diameter", self.inputs.bullet_diameter)?;
        require_positive("muzzle_velocity", self.inputs.muzzle_velocity)?;
        // MBA-1356: cd_scale multiplies the custom-deck Cd; a non-finite or non-positive value
        // would zero/invert drag or poison integration with NaN. Required unconditionally (not
        // gated on custom_drag_table) since the default (1.0) is always valid and a caller could
        // set an invalid scale without ever setting a table. The [0.5, 2.0] "unusual" range is a
        // CLI-level warning (Task 2), not a hard engine rule.
        require_positive("cd_scale", self.inputs.cd_scale)?;

        require_finite("muzzle_angle", self.inputs.muzzle_angle)?;
        require_finite("azimuth_angle", self.inputs.azimuth_angle)?;
        require_finite("shooting_angle", self.inputs.shooting_angle)?;
        require_finite("cant_angle", self.inputs.cant_angle)?;
        require_finite("muzzle_height", self.inputs.muzzle_height)?;

        // MBA-1359: a deliberate zero POI offset is a small linear offset at the zero range
        // (fractions of an inch to a few inches). |1 m| is far beyond any plausible
        // deliberate offset and almost certainly a unit error (inches/cm passed as meters).
        for (name, value) in [
            ("zero_poi_vertical_m", self.inputs.zero_poi_vertical_m),
            ("zero_poi_horizontal_m", self.inputs.zero_poi_horizontal_m),
        ] {
            require_finite(name, value)?;
            if value.abs() >= 1.0 {
                return Err(BallisticsError::from(format!(
                    "{name} must be smaller than 1.0 m in magnitude (it is a linear POI \
                     offset at the zero range, in meters)"
                )));
            }
        }

        // MBA-1396: a lateral sight-mount offset is physically bounded by rail/mount
        // geometry (an inch or two). |0.5 m| is almost certainly a unit error
        // (inches/mm passed as meters).
        require_finite(
            "sight_offset_lateral_m",
            self.inputs.sight_offset_lateral_m,
        )?;
        if self.inputs.sight_offset_lateral_m.abs() >= 0.5 {
            return Err(BallisticsError::from(
                "sight_offset_lateral_m must be smaller than 0.5 m in magnitude (it is \
                 the lateral sight-to-bore mount offset, in meters)",
            ));
        }

        // Negative infinity is the documented ignore-ground sentinel. NaN and positive infinity
        // make the loop condition meaningless and are rejected.
        if !(self.inputs.ground_threshold.is_finite()
            || self.inputs.ground_threshold == f64::NEG_INFINITY)
        {
            return Err(BallisticsError::from(
                "ground_threshold must be finite or negative infinity",
            ));
        }

        match &self.wind_sock {
            Some(wind_sock) => wind_sock
                .validate_segments()
                .map_err(BallisticsError::from)?,
            None => {
                require_finite("wind.speed", self.wind.speed)?;
                require_finite("wind.direction", self.wind.direction)?;
                require_finite("wind.vertical_speed", self.wind.vertical_speed)?;
            }
        }

        require_finite("atmosphere.temperature", self.atmosphere.temperature)?;
        require_finite("atmosphere.pressure", self.atmosphere.pressure)?;
        require_finite("atmosphere.humidity", self.atmosphere.humidity)?;
        require_finite("atmosphere.altitude", self.atmosphere.altitude)?;

        require_positive("max_range", self.max_range)?;
        // Adaptive RK45 owns its step size; the caller-provided fixed step is used only by Euler
        // and fixed RK4.
        if !self.inputs.use_rk4 || !self.inputs.use_adaptive_rk45 {
            require_positive("time_step", self.time_step)?;
        }

        if self.inputs.enable_trajectory_sampling {
            require_finite("sight_height", self.inputs.sight_height)?;
            require_positive("sample_interval", self.inputs.sample_interval)?;
            projected_sample_count(self.max_range, self.inputs.sample_interval)?;
        }

        // MBA-1403: the target-plane drops reference divides sampled drop by
        // cos(shooting_angle); at or beyond 90 degrees the transform is undefined.
        // Gated on the non-default mode only, so LOS-mode validation is unchanged.
        if self.inputs.drops_reference == DropsReference::Target {
            require_finite("target_height", self.inputs.target_height)?;
            if self.inputs.shooting_angle.cos() <= 1e-9 {
                return Err(BallisticsError::from(
                    "drops reference 'target' is undefined for shooting angles at or beyond 90 degrees",
                ));
            }
        }

        if self.inputs.enable_coriolis {
            require_finite("shot_azimuth", self.inputs.shot_azimuth)?;
            if let Some(latitude) = self.inputs.latitude {
                require_finite("latitude", latitude)?;
            }
        }

        Ok(())
    }

    /// Public solve results must never report success with NaN or infinity, nor with values a
    /// physical trajectory cannot produce: a negative terminal downrange distance, time of
    /// flight, speed, or energy (MBA-1293 — a stiff-input integration explosion reported
    /// `Ok(max_range: -50.59)`). The input gate catches malformed scalar state; this
    /// postcondition also covers overflow and malformed optional tables/segments without
    /// imposing arbitrary upper bounds on otherwise finite inputs.
    fn validate_result_sanity(&self, result: &TrajectoryResult) -> Result<(), BallisticsError> {
        let require_finite = |name: &str, value: f64| {
            if value.is_finite() {
                Ok(())
            } else {
                Err(BallisticsError::from(format!(
                    "trajectory result contains non-finite {name}"
                )))
            }
        };
        let require_non_negative = |name: &str, value: f64| {
            if value >= 0.0 {
                Ok(())
            } else {
                Err(BallisticsError::from(format!(
                    "trajectory result contains non-physical negative {name} ({value})"
                )))
            }
        };
        let require_indexed_finite = |collection: &str, index: usize, field: &str, value: f64| {
            if value.is_finite() {
                Ok(())
            } else {
                Err(BallisticsError::from(format!(
                    "trajectory result contains non-finite {collection}[{index}].{field}"
                )))
            }
        };
        let require_indexed_non_negative =
            |collection: &str, index: usize, field: &str, value: f64| {
                if value >= 0.0 {
                    Ok(())
                } else {
                    Err(BallisticsError::from(format!(
                        "trajectory result contains non-physical negative {collection}[{index}].{field} ({value})"
                    )))
                }
            };

        require_finite("max_range", result.max_range)?;
        require_finite("max_height", result.max_height)?;
        require_finite("time_of_flight", result.time_of_flight)?;
        require_finite("impact_velocity", result.impact_velocity)?;
        require_finite("impact_energy", result.impact_energy)?;
        require_finite("projectile_mass_kg", result.projectile_mass_kg)?;
        require_finite(
            "line_of_sight_height_m",
            result.line_of_sight_height_m,
        )?;
        require_finite(
            "station_speed_of_sound_mps",
            result.station_speed_of_sound_mps,
        )?;

        // The solve starts at x = 0 and only ever fires downrange, so these scalars are
        // non-negative for every physically meaningful trajectory. (max_height is exempt:
        // points can legitimately sit below y = 0 with an elevated muzzle.)
        require_non_negative("max_range", result.max_range)?;
        require_non_negative("time_of_flight", result.time_of_flight)?;
        require_non_negative("impact_velocity", result.impact_velocity)?;
        require_non_negative("impact_energy", result.impact_energy)?;
        require_non_negative("projectile_mass_kg", result.projectile_mass_kg)?;
        require_non_negative(
            "station_speed_of_sound_mps",
            result.station_speed_of_sound_mps,
        )?;

        for (index, point) in result.points.iter().enumerate() {
            require_indexed_finite("points", index, "time", point.time)?;
            require_indexed_finite("points", index, "position.x", point.position.x)?;
            require_indexed_finite("points", index, "position.y", point.position.y)?;
            require_indexed_finite("points", index, "position.z", point.position.z)?;
            require_indexed_finite(
                "points",
                index,
                "velocity_magnitude",
                point.velocity_magnitude,
            )?;
            require_indexed_finite("points", index, "kinetic_energy", point.kinetic_energy)?;
            require_indexed_non_negative("points", index, "time", point.time)?;
            require_indexed_non_negative(
                "points",
                index,
                "velocity_magnitude",
                point.velocity_magnitude,
            )?;
            require_indexed_non_negative("points", index, "kinetic_energy", point.kinetic_energy)?;
        }

        if let Some(samples) = &result.sampled_points {
            for (index, sample) in samples.iter().enumerate() {
                require_indexed_finite("sampled_points", index, "distance_m", sample.distance_m)?;
                require_indexed_finite("sampled_points", index, "drop_m", sample.drop_m)?;
                require_indexed_finite(
                    "sampled_points",
                    index,
                    "wind_drift_m",
                    sample.wind_drift_m,
                )?;
                require_indexed_finite(
                    "sampled_points",
                    index,
                    "velocity_mps",
                    sample.velocity_mps,
                )?;
                require_indexed_finite("sampled_points", index, "energy_j", sample.energy_j)?;
                require_indexed_finite("sampled_points", index, "time_s", sample.time_s)?;
            }
        }

        for (name, value) in [
            ("min_pitch_damping", result.min_pitch_damping),
            ("transonic_mach", result.transonic_mach),
            ("max_yaw_angle", result.max_yaw_angle),
            ("max_precession_angle", result.max_precession_angle),
        ] {
            if let Some(value) = value {
                require_finite(name, value)?;
            }
        }

        if let Some(state) = result.angular_state {
            for (name, value) in [
                ("angular_state.pitch_angle", state.pitch_angle),
                ("angular_state.yaw_angle", state.yaw_angle),
                ("angular_state.pitch_rate", state.pitch_rate),
                ("angular_state.yaw_rate", state.yaw_rate),
                ("angular_state.precession_angle", state.precession_angle),
                ("angular_state.nutation_phase", state.nutation_phase),
            ] {
                require_finite(name, value)?;
            }
        }

        if let Some(jump) = result.aerodynamic_jump {
            for (name, value) in [
                ("aerodynamic_jump.vertical_jump_moa", jump.vertical_jump_moa),
                (
                    "aerodynamic_jump.horizontal_jump_moa",
                    jump.horizontal_jump_moa,
                ),
                ("aerodynamic_jump.jump_angle_rad", jump.jump_angle_rad),
                (
                    "aerodynamic_jump.magnus_component_moa",
                    jump.magnus_component_moa,
                ),
                ("aerodynamic_jump.yaw_component_moa", jump.yaw_component_moa),
                (
                    "aerodynamic_jump.stabilization_factor",
                    jump.stabilization_factor,
                ),
            ] {
                require_finite(name, value)?;
            }
        }

        Ok(())
    }

    /// Integration methods store the pre-step state in `points`. Validate each newly accepted
    /// state as well, otherwise a poisoned final step could terminate the loop and leave only the
    /// previous finite point in an apparently successful result.
    ///
    /// Beyond finiteness, an accepted state must respect the physical speed budget: drag is
    /// dissipative (it drives the projectile toward the wind frame, never past it), Magnus and
    /// Coriolis act perpendicular to the velocity and do no work, and gravity adds at most g*t.
    /// Ground-frame speed therefore cannot legitimately exceed muzzle speed + strongest wind +
    /// g*t. Exceeding that budget means the integrator itself diverged — for stiff inputs the
    /// minimum-step RK45 acceptance can multiply speed by orders of magnitude in one step
    /// (MBA-1293: 13x and a sign reversal in a single 1 microsecond step) — so the solve must
    /// fail rather than report the garbage as `Ok`.
    fn validate_integration_state(
        &self,
        position: &Vector3<f64>,
        velocity: &Vector3<f64>,
        time: f64,
    ) -> Result<(), BallisticsError> {
        if !(position.iter().all(|value| value.is_finite())
            && velocity.iter().all(|value| value.is_finite())
            && time.is_finite())
        {
            return Err(BallisticsError::from(
                "trajectory integration produced a non-finite state (often from physically \
                 extreme inputs — e.g. an absurd bore/muzzle height placing the launch far \
                 from sea level, or a degenerate atmosphere; check those inputs, or set \
                 --altitude explicitly)",
            ));
        }

        let speed = velocity.magnitude();
        let budget = self.speed_budget(time);
        if speed > budget {
            return Err(BallisticsError::from(format!(
                "trajectory integration diverged: speed {speed:.3e} m/s at t={time:.6}s exceeds \
                 the physical budget of {budget:.3e} m/s"
            )));
        }
        Ok(())
    }

    /// Ceiling on ground-frame speed a physical trajectory can reach by time `t` (see
    /// [`Self::validate_integration_state`]). The factor-2 slack absorbs boundary-layer
    /// wind-shear amplification and integrator transients; genuine divergence clears the
    /// budget by orders of magnitude.
    fn speed_budget(&self, time: f64) -> f64 {
        let scalar_wind = self.wind.speed.abs() + self.wind.vertical_speed.abs();
        let wind_bound = match &self.wind_sock {
            Some(sock) => scalar_wind.max(sock.max_speed_mps()),
            None => scalar_wind,
        };
        2.0 * (self.inputs.muzzle_velocity + wind_bound + 10.0)
            + crate::constants::G_ACCEL_MPS2 * time
    }

    /// Store one public trajectory point without exceeding the per-solve resource budget.
    fn push_trajectory_point(
        &self,
        points: &mut Vec<TrajectoryPoint>,
        point: TrajectoryPoint,
    ) -> Result<(), BallisticsError> {
        if points.len() >= self.max_trajectory_points {
            return Err(BallisticsError::from(format!(
                "trajectory point limit of {} exceeded",
                self.max_trajectory_points
            )));
        }
        points.push(point);
        Ok(())
    }

    /// Supply downrange-segmented wind. Each segment is `(speed_kmh, angle_deg,
    /// until_distance_m)`; the wind for a given downrange distance is the first
    /// segment whose `until_distance_m` exceeds it (a step function), and wind is
    /// zero beyond the last segment. An empty list clears segmented wind (reverts
    /// to the scalar `wind`). The angle convention matches `WindConditions`
    /// (0 = headwind, 90 = from the right).
    pub fn set_wind_segments(&mut self, segments: Vec<crate::wind::WindSegment>) {
        self.wind_sock = if segments.is_empty() {
            None
        } else {
            Some(crate::wind::WindSock::new(segments))
        };
    }

    /// Supply downrange-segmented atmosphere (MBA-1137). Each segment is
    /// `(temp_c, pressure_hpa, humidity_percent, until_distance_m)`, defined at the shooter base
    /// altitude; the per-substep local-atmosphere recompute selects the active zone by downrange
    /// distance (first zone whose `until_distance_m` exceeds it; the last zone is held beyond the
    /// final threshold). The zone's base conditions are composed with the vertical altitude lapse
    /// via `get_local_atmosphere_humid`, so a steeply-arcing shot still sees the y-lapse on top of
    /// the zone base. An empty list clears segmented atmosphere (reverts to the resolved
    /// single-station conditions).
    pub fn set_atmo_segments(&mut self, segments: Vec<crate::atmosphere::AtmoSegment>) {
        self.atmo_sock = if segments.is_empty() {
            None
        } else {
            Some(crate::atmosphere::AtmoSock::new(segments))
        };
    }

    /// Effective initial launch direction `(elevation, azimuth)` in radians, including
    /// the aerodynamic-jump muzzle perturbation when `enable_aerodynamic_jump` is set.
    ///
    /// Aerodynamic jump is the fixed angular departure imparted as the projectile
    /// transitions from the constrained bore to free flight; applying it as an initial
    /// launch-angle offset is the physically correct integration point. Returns the bare
    /// `(muzzle_angle, azimuth_angle)` when the flag is off, so a default solve is
    /// numerically identical to pre-feature behavior. (MBA-959)
    fn launch_angles_from(
        &self,
        aj: Option<&crate::aerodynamic_jump::AerodynamicJumpComponents>,
    ) -> (f64, f64) {
        let (mut elev, mut azim) = (self.inputs.muzzle_angle, self.inputs.azimuth_angle);
        // MBA-1286: cant rotates the sight-frame aim offsets about the line of sight.
        // Positive = clockwise from the shooter: the upward zero correction leaks right
        // (+z) and shrinks by cos(cant) -> POI right and low. Exactly-0.0 skips all float
        // ops so un-canted solves stay bit-identical. Aerodynamic jump is added AFTER the
        // rotation: it arises at bore exit from crosswind/spin in the ground frame, not
        // from the rifle's sight geometry.
        if self.inputs.cant_angle != 0.0 {
            let (sin_c, cos_c) = self.inputs.cant_angle.sin_cos();
            let (e0, a0) = (elev, azim);
            elev = e0 * cos_c - a0 * sin_c;
            azim = a0 * cos_c + e0 * sin_c;
        }
        match aj {
            Some(c) => {
                // vertical_/horizontal_jump_moa ARE the jump angles expressed in MOA.
                const MOA_PER_RAD: f64 = 3437.7467707849;
                (
                    elev + c.vertical_jump_moa / MOA_PER_RAD,
                    azim + c.horizontal_jump_moa / MOA_PER_RAD,
                )
            }
            None => (elev, azim),
        }
    }

    /// Compute the aerodynamic-jump components for the current inputs, or `None` when the
    /// feature is disabled / inputs are degenerate.
    ///
    /// Uses Bryan Litz's crosswind aerodynamic-jump estimator
    /// (`Y = 0.01*Sg - 0.0024*L + 0.032` MOA/mph) fed by the engine's own Miller Sg.
    /// Aerodynamic jump is a vertical effect, so only the elevation is perturbed.
    /// The estimator is a regression best near Sg ~ 1.75 — see MBA-959.
    fn aerodynamic_jump_components(
        &self,
    ) -> Option<crate::aerodynamic_jump::AerodynamicJumpComponents> {
        if !self.inputs.enable_aerodynamic_jump {
            return None;
        }
        // Reject degenerate/non-finite inputs before they can reach the launch angle.
        // A bare `<= 0.0` test lets NaN through (NaN comparisons are always false), and a
        // NaN/Inf here would poison the muzzle angle and collapse the whole trajectory.
        let diameter_m = self.inputs.bullet_diameter;
        if !(self.inputs.twist_rate.is_finite()
            && self.inputs.twist_rate != 0.0
            && diameter_m.is_finite()
            && diameter_m > 0.0
            && self.inputs.bullet_length.is_finite()
            && self.inputs.bullet_length > 0.0
            && self.inputs.muzzle_velocity.is_finite())
        {
            return None;
        }

        // Engine's own gyroscopic (Miller) stability factor — same Sg shown elsewhere.
        let (_, _, temp_c, pressure_hpa) = self.resolved_atmosphere();
        let sg = crate::stability::compute_stability_coefficient(
            &self.inputs,
            (self.atmosphere.altitude, temp_c, pressure_hpa, 0.0),
        );
        if !(sg.is_finite() && sg > 0.0) {
            return None;
        }
        let length_calibers = self.inputs.bullet_length / diameter_m;

        // Crosswind-from-the-right (mph) for Litz's estimator. Wind direction uses the
        // wind-FROM convention (0 = headwind, +90deg = from the right), matching the
        // fast-integrate path (fast_trajectory::aerodynamic_jump_launch_offset_rad) and
        // the lateral windage sign, so a from-the-right wind on a right-twist barrel
        // jumps the impact UP and drifts it left.
        const MS_TO_MPH: f64 = 2.236_936_292_054_4;
        let crosswind_from_right_mps = if let Some(sock) = &self.wind_sock {
            -sock.vector_for_range_stateless(0.0)[2]
        } else {
            self.wind.speed * self.wind.direction.sin()
        };
        let crosswind_from_right_mph = crosswind_from_right_mps * MS_TO_MPH;

        let vertical_jump_moa = crate::aerodynamic_jump::litz_crosswind_jump_moa(
            sg,
            length_calibers,
            crosswind_from_right_mph,
            self.inputs.is_twist_right,
        );
        if !vertical_jump_moa.is_finite() {
            return None;
        }

        const MOA_PER_RAD: f64 = 3437.7467707849;
        Some(crate::aerodynamic_jump::AerodynamicJumpComponents {
            vertical_jump_moa,
            // Aerodynamic jump is a vertical effect; the Litz estimator has no horizontal term.
            horizontal_jump_moa: 0.0,
            jump_angle_rad: vertical_jump_moa.abs() / MOA_PER_RAD,
            magnus_component_moa: 0.0,
            yaw_component_moa: 0.0,
            stabilization_factor: (sg / 1.5).clamp(0.0, 1.0),
        })
    }

    fn resolved_atmosphere(&self) -> (f64, f64, f64, f64) {
        let (temp_c, pressure_hpa) = match self.station_atmosphere_resolution {
            StationAtmosphereResolution::LegacyDefaultSentinels => {
                crate::atmosphere::resolve_station_conditions(
                    self.atmosphere.temperature,
                    self.atmosphere.pressure,
                    self.atmosphere.altitude,
                )
            }
            StationAtmosphereResolution::Authoritative => {
                (self.atmosphere.temperature, self.atmosphere.pressure)
            }
        };
        let (density, speed_of_sound) = crate::atmosphere::calculate_atmosphere(
            self.atmosphere.altitude,
            Some(temp_c),
            Some(pressure_hpa),
            self.atmosphere.humidity,
        );
        (density, speed_of_sound, temp_c, pressure_hpa)
    }

    fn precession_nutation_params(
        &self,
        velocity_mps: f64,
        air_density_kg_m3: f64,
        speed_of_sound_mps: f64,
    ) -> PrecessionNutationParams {
        let (spin_inertia, transverse_inertia) = self.precession_nutation_inertias;
        let spin_rate_rad_s = if self.inputs.twist_rate > 0.0 {
            let velocity_fps = velocity_mps * 3.28084;
            let twist_rate_ft = self.inputs.twist_rate / 12.0;
            (velocity_fps / twist_rate_ft) * 2.0 * std::f64::consts::PI
        } else {
            0.0
        };

        PrecessionNutationParams {
            mass_kg: self.inputs.bullet_mass,
            caliber_m: self.inputs.bullet_diameter,
            length_m: self.inputs.bullet_length,
            spin_rate_rad_s,
            spin_inertia,
            transverse_inertia,
            velocity_mps,
            air_density_kg_m3,
            mach: velocity_mps / speed_of_sound_mps,
            pitch_damping_coeff: PitchDampingCoefficients::default().subsonic,
            nutation_damping_factor: 0.05,
        }
    }

    /// Append the exact state at the earliest boundary crossed by the final integration step.
    ///
    /// Each solver stores its pre-step state. Keeping only that point makes early ground and time
    /// exits indistinguishable from ordinary integration knots, and historically left the
    /// reported endpoint one step short. Interpolating all supported boundaries here gives every
    /// solver one explicit terminal point and one authoritative termination reason.
    fn append_terminal_endpoint(
        &self,
        points: &mut Vec<TrajectoryPoint>,
        post_position: Vector3<f64>,
        post_velocity: Vector3<f64>,
        post_time: f64,
        max_height: &mut f64,
    ) -> Result<TrajectoryTermination, BallisticsError> {
        let previous = points
            .last()
            .cloned()
            .ok_or_else(|| BallisticsError::from("No trajectory points generated"))?;

        let mut crossings = Vec::with_capacity(3);
        if previous.position.x < self.max_range && post_position.x >= self.max_range {
            let span = post_position.x - previous.position.x;
            if span.is_finite() && span > 0.0 {
                crossings.push((
                    (self.max_range - previous.position.x) / span,
                    TrajectoryTermination::MaxRange,
                ));
            }
        }
        if self.inputs.ground_threshold.is_finite()
            && previous.position.y > self.inputs.ground_threshold
            && post_position.y <= self.inputs.ground_threshold
        {
            let span = post_position.y - previous.position.y;
            if span.is_finite() && span < 0.0 {
                crossings.push((
                    (self.inputs.ground_threshold - previous.position.y) / span,
                    TrajectoryTermination::GroundThreshold,
                ));
            }
        }
        if previous.time < TRAJECTORY_TIME_LIMIT_S && post_time >= TRAJECTORY_TIME_LIMIT_S {
            let span = post_time - previous.time;
            if span.is_finite() && span > 0.0 {
                crossings.push((
                    (TRAJECTORY_TIME_LIMIT_S - previous.time) / span,
                    TrajectoryTermination::TimeLimit,
                ));
            }
        }

        let (fraction, termination) = crossings
            .into_iter()
            .filter(|(fraction, _)| fraction.is_finite() && (0.0..=1.0).contains(fraction))
            .min_by(|left, right| {
                let priority = |termination: TrajectoryTermination| match termination {
                    TrajectoryTermination::GroundThreshold => 0,
                    TrajectoryTermination::MaxRange => 1,
                    TrajectoryTermination::TimeLimit => 2,
                    TrajectoryTermination::VelocityFloor => 3,
                };
                left.0
                    .total_cmp(&right.0)
                    .then_with(|| priority(left.1).cmp(&priority(right.1)))
            })
            .ok_or_else(|| {
                BallisticsError::from(
                    "trajectory integration stopped without crossing a supported boundary",
                )
            })?;

        let mut position = previous.position + (post_position - previous.position) * fraction;
        match termination {
            TrajectoryTermination::MaxRange => position.x = self.max_range,
            TrajectoryTermination::GroundThreshold => {
                position.y = self.inputs.ground_threshold;
            }
            TrajectoryTermination::TimeLimit | TrajectoryTermination::VelocityFloor => {}
        }
        let velocity_magnitude = previous.velocity_magnitude
            + (post_velocity.magnitude() - previous.velocity_magnitude) * fraction;
        let mut time = previous.time + (post_time - previous.time) * fraction;
        if termination == TrajectoryTermination::TimeLimit {
            time = TRAJECTORY_TIME_LIMIT_S;
        }
        let kinetic_energy =
            0.5 * self.inputs.bullet_mass * velocity_magnitude * velocity_magnitude;

        if position.y > *max_height {
            *max_height = position.y;
        }
        let terminal_point = TrajectoryPoint {
            time,
            position,
            velocity_magnitude,
            kinetic_energy,
            drag_coefficient: None,
        };
        if terminal_point.position.x < previous.position.x {
            return Err(BallisticsError::from(
                "trajectory terminal state reversed downrange before the crossed boundary",
            ));
        }
        if terminal_point.position.x == previous.position.x {
            // A very early ground/time crossing can be distinct in time but less than one ULP
            // downrange. There is no representable range at which to retain both states, so make
            // the terminal state authoritative instead of creating a duplicate-X trajectory that
            // the checked observation API must reject.
            let last = points.last_mut().ok_or_else(|| {
                BallisticsError::from("trajectory points disappeared during terminal finalization")
            })?;
            *last = terminal_point;
        } else {
            self.push_trajectory_point(points, terminal_point)?;
        }
        Ok(termination)
    }

    fn gravity_acceleration(&self) -> Vector3<f64> {
        let theta = self.inputs.shooting_angle;
        Vector3::new(
            -crate::constants::G_ACCEL_MPS2 * theta.sin(),
            -crate::constants::G_ACCEL_MPS2 * theta.cos(),
            0.0,
        )
    }

    fn get_wind_at_altitude(&self, altitude_m: f64) -> Vector3<f64> {
        // Scale the operative surface wind by the boundary-layer multiplier. `altitude_m` is the
        // bullet's height relative to the muzzle (McCoy Y). The multiplier is floored at 1.0, so
        // flat-fire trajectories keep ~full wind and only high-arcing shots see increased wind.
        //
        // We build the vector with THIS solver's non-shear sign convention (X=-cos, Z=-sin; see
        // the `wind_vector` used in solve_rk4/solve_euler, matching WindSock) and scale it, so that
        // "shear on" equals "shear off" * ratio (ratio == 1.0 for flat fire). An earlier revision
        // attenuated the wind near the line of sight and flipped its sign relative to the non-shear
        // path; this keeps them sign-consistent.
        // Map the requested model name to the boundary-layer model (MBA-965).
        // Names match wind_shear::get_wind_at_position. Unknown strings should
        // never reach here (the CLI parses an enum), but default to PowerLaw to
        // preserve the historical "exponential" behaviour for any caller that
        // forwards an unexpected value.
        let model = match self.inputs.wind_shear_model.as_str() {
            "logarithmic" => WindShearModel::Logarithmic,
            "power_law" | "powerlaw" | "exponential" => WindShearModel::PowerLaw,
            "ekman_spiral" | "ekman" => WindShearModel::EkmanSpiral,
            "custom_layers" | "custom" => WindShearModel::CustomLayers,
            _ => WindShearModel::PowerLaw,
        };
        let speed_ratio = crate::wind_shear::boundary_layer_speed_ratio(altitude_m, model);

        // 0deg = headwind, 90deg = from the right (McCoy wind-FROM convention, matching
        // WindConditions / WindSock); wind enters drag via velocity - wind.
        //
        // MBA-728: the horizontal vector is built with vertical=0.0 and scaled by speed_ratio,
        // then wind.vertical_speed is added back UNSCALED — boundary-layer shear scales
        // horizontal wind only, vertical passes through as-is.
        crate::wind::wind_vector(self.wind.speed, self.wind.direction, 0.0) * speed_ratio
            + Vector3::new(0.0, self.wind.vertical_speed, 0.0)
    }

    pub fn solve(&self) -> Result<TrajectoryResult, BallisticsError> {
        self.validate_for_solve()?;
        let mut result = if self.inputs.use_rk4 {
            if self.inputs.use_adaptive_rk45 {
                self.solve_rk45()?
            } else {
                self.solve_rk4()?
            }
        } else {
            self.solve_euler()?
        };
        self.apply_spin_drift(&mut result);
        self.validate_result_sanity(&result)?;
        Ok(result)
    }

    /// Gyroscopic spin drift via the empirical Litz model, applied in the engine
    /// (not the WASM formatter) so it covers Euler/RK4/RK45 and all consumers.
    /// Uses the canonical SI fields and converts to grains/inches correctly,
    /// avoiding the kg/m-vs-grains/in unit bug in `calculate_enhanced_spin_drift`.
    /// Frame (McCoy): Z = lateral (windage), so drift adds to `position.z`.
    fn apply_spin_drift(&self, result: &mut TrajectoryResult) {
        if !self.inputs.use_enhanced_spin_drift {
            return;
        }
        let d_in = self.inputs.bullet_diameter / 0.0254; // m -> in
        let m_gr = self.inputs.bullet_mass / crate::constants::GRAINS_TO_KG; // kg -> grains
        let twist_in = self.inputs.twist_rate; // inches/turn
        if d_in <= 0.0 || m_gr <= 0.0 || twist_in <= 0.0 {
            return;
        }

        // MBA-1134 (rank 31): single source of truth for the muzzle Sg —
        // stability::compute_stability_coefficient via spin_drift::effective_sg_from_inputs. This
        // ADDS the (v/2800)^(1/3) muzzle-velocity term the bare miller_stability() lacked, so the
        // spin-drift Sg now matches the reported SG and the aerodynamic-jump Sg. The linear Miller
        // density correction ((T/T0)*(P0/P), a no-op at sea-level standard) and the 4.5-caliber
        // length fallback are handled inside effective_sg_from_inputs.
        let sg = self.effective_spin_drift_sg();

        for p in result.points.iter_mut() {
            if p.time <= 0.0 {
                continue;
            }
            // Canonical Litz drift, shared with the fast / Monte-Carlo path (spin_drift::litz_*).
            p.position.z +=
                crate::spin_drift::litz_drift_meters(sg, p.time, self.inputs.is_twist_right);
        }

        // sampled_points are snapshotted from the PRE-drift trajectory inside each solver, so the
        // sampled wind_drift_m column would omit the spin drift that result.points carry. Apply
        // the same canonical Litz drift to keep the two user-facing outputs consistent.
        if let Some(samples) = result.sampled_points.as_mut() {
            for s in samples.iter_mut() {
                if s.time_s <= 0.0 {
                    continue;
                }
                s.wind_drift_m +=
                    crate::spin_drift::litz_drift_meters(sg, s.time_s, self.inputs.is_twist_right);
            }
        }
    }

    /// Muzzle gyroscopic stability Sg used by the empirical Litz spin-drift post-process
    /// (MBA-1134). Extracted so the exact value is unit-testable and provably identical to the Sg
    /// the fast / Monte-Carlo path uses — both go through
    /// [`crate::spin_drift::effective_sg_from_inputs`] with the resolved muzzle atmosphere.
    fn effective_spin_drift_sg(&self) -> f64 {
        let (_, _, temp_c, press_hpa) = self.resolved_atmosphere();
        crate::spin_drift::effective_sg_from_inputs(&self.inputs, temp_c, press_hpa)
    }

    /// Bore muzzle position at t=0 (bore-origin frame, `muzzle_height` above ground).
    /// With cant the rifle rotates about the LINE OF SIGHT, so the bore — sight_height
    /// below the sight — swings laterally by `-sight_height*sin(cant)` (left of the aim
    /// plane for clockwise cant) and rises by `sight_height*(1-cos(cant))` toward the
    /// pivot. (MBA-1286)
    ///
    /// A lateral sight-mount offset (MBA-1396) additionally displaces the bore
    /// `sight_offset_lateral_m` from the LOS (a sight mounted right of the bore puts the
    /// bore left of the sight line). Cant rotates the WHOLE bore-to-sight displacement
    /// vector rigidly about the LOS, so the vertical drop `-sight_height` and the lateral
    /// offset `-sight_offset_lateral_m` couple: rotating `(u, w) = (-sh, -off)` by the
    /// cant angle gives `y += sh*(1-cos) + off*sin` and `z = -sh*sin - off*cos`. Both the
    /// cant-only case (off = 0) and the offset-only case (cos = 1, sin = 0) collapse to
    /// the historical terms, so only a simultaneously canted AND offset rifle sees the
    /// coupling. The windage-zero convergence that makes the trajectory cross the LOS at
    /// the zero range lives in `windage_zero_bias_rad`, not here. Exactly-0.0 cant AND
    /// offset return the historical position (bit-identical).
    fn initial_position(&self) -> Vector3<f64> {
        if self.inputs.cant_angle == 0.0 && self.inputs.sight_offset_lateral_m == 0.0 {
            return Vector3::new(0.0, self.inputs.muzzle_height, 0.0);
        }
        let (sin_c, cos_c) = self.inputs.cant_angle.sin_cos();
        let sh = self.inputs.sight_height;
        let off = self.inputs.sight_offset_lateral_m;
        Vector3::new(
            0.0,
            self.inputs.muzzle_height + sh * (1.0 - cos_c) + off * sin_c,
            -sh * sin_c - off * cos_c,
        )
    }

    /// Shared post-integration sampling used by all three integrators (MBA-1403).
    ///
    /// Builds the `TrajectoryData`/`TrajectoryOutputs` pair and runs
    /// [`sample_trajectory`] when `enable_trajectory_sampling` is set. This is the
    /// single place the sampler's LOS datum (`target_vertical_height_m`) is chosen,
    /// replacing three per-integrator copies that had to be edited in parallel.
    fn build_sampled_points(
        &self,
        points: &[TrajectoryPoint],
        max_height: f64,
        transonic_distances: Vec<f64>,
        mach_transitions: &MachTransitionTracker,
    ) -> Result<Option<Vec<TrajectorySample>>, BallisticsError> {
        if !self.inputs.enable_trajectory_sampling {
            return Ok(None);
        }

        let last_point = points.last().ok_or("No trajectory points generated")?;
        let trajectory_data = TrajectoryData {
            times: points.iter().map(|p| p.time).collect(),
            positions: points.iter().map(|p| p.position).collect(),
            velocities: points
                .iter()
                .map(|p| {
                    // Reconstruct velocity vectors from magnitude (approximate)
                    Vector3::new(0.0, 0.0, p.velocity_magnitude)
                })
                .collect(),
            transonic_distances, // populated by the integrator at each Mach-threshold crossing
            mach_1_2_distance_m: mach_transitions.mach_1_2_distance_m,
            mach_1_0_distance_m: mach_transitions.mach_1_0_distance_m,
            mach_0_9_distance_m: mach_transitions.mach_0_9_distance_m,
        };

        // For LOS calculation in ground-referenced coordinates:
        // sight_position_m is the sight's actual y-position above ground
        // (muzzle_height + sight_height, not just sight_height)
        // For flat shots, target is at same height as the sight (horizontal LOS)
        let sight_position_m = self.inputs.muzzle_height + self.inputs.sight_height;
        let target_reference = self.inputs.drops_reference == DropsReference::Target;
        // MBA-1403: only the (new, non-default) target mode reads `target_height` here,
        // and only when a real (non-zero) height was supplied — the sampler's LOS then
        // slopes toward the actual target instead of staying at the sight height. LOS
        // mode keeps the historical datum byte-identically, and zero-solve
        // `target_height` semantics (an explicit function argument, never this field)
        // are untouched.
        let target_vertical_height_m = if target_reference && self.inputs.target_height != 0.0 {
            self.inputs.target_height
        } else {
            sight_position_m
        };
        let outputs = TrajectoryOutputs {
            target_distance_horiz_m: last_point.position.x, // X is downrange
            target_vertical_height_m,
            time_of_flight_s: last_point.time,
            max_ord_dist_horiz_m: max_height,
            sight_height_m: sight_position_m,
        };

        // Sample at specified intervals
        let mut samples = sample_trajectory(
            &trajectory_data,
            &outputs,
            self.inputs.sample_interval,
            self.inputs.bullet_mass,
        )?;

        // MBA-1403 target-plane reference: scale the LOS-perpendicular drop to vertical
        // in the target plane. Scaling by the positive constant 1/cos preserves signs
        // and ordering, so the sampler's zero-crossing/apex flags are unaffected.
        // validate_for_solve rejects target mode at |shooting_angle| >= 90 degrees, so
        // cos_theta is strictly positive here. Downstream mil/moa conversions derive
        // from drop_m, so they inherit the transform without further edits.
        if target_reference {
            let cos_theta = self.inputs.shooting_angle.cos();
            for sample in &mut samples {
                sample.drop_m /= cos_theta;
            }
        }
        Ok(Some(samples))
    }

    fn solve_euler(&self) -> Result<TrajectoryResult, BallisticsError> {
        // Simple trajectory integration using Euler method
        let mut time = 0.0;
        // Bullet starts at the BORE position, which is muzzle_height above ground
        // The sight is sight_height ABOVE the bore, so we don't add sight_height here
        // cant-adjusted via initial_position (MBA-1286)
        let mut position = self.initial_position();
        // Calculate initial velocity components with both elevation and azimuth
        // McCoy coordinate system: X=downrange, Y=vertical, Z=lateral (right)
        // Launch direction includes the aerodynamic-jump muzzle perturbation when enabled
        // (a no-op returning the bare muzzle/azimuth angles otherwise). MBA-959. Computed
        // once here and reused for the result so it isn't evaluated twice per solve.
        let aj_components = self.aerodynamic_jump_components();
        let (launch_elev, launch_azim) = self.launch_angles_from(aj_components.as_ref());
        let horizontal_velocity = self.inputs.muzzle_velocity * launch_elev.cos();
        let mut velocity = Vector3::new(
            horizontal_velocity * launch_azim.cos(), // X: downrange (forward)
            self.inputs.muzzle_velocity * launch_elev.sin(), // Y: vertical component
            horizontal_velocity * launch_azim.sin(), // Z: lateral (side deviation)
        );

        let mut points = Vec::new();
        let mut max_height = position.y;
        let mut min_pitch_damping = f64::INFINITY; // Track minimum pitch damping coefficient
        let mut transonic_mach = None; // Track when we enter transonic
                                       // Downrange distances where the projectile crosses Mach 1.2 (transonic) then Mach 1.0
                                       // (subsonic), so the sampled trajectory output can flag those transitions
                                       // (trajectory_sampling::add_trajectory_flags consumes this).
        let mut transonic_distances: Vec<f64> = Vec::new();
        let mut mach_transitions = MachTransitionTracker::default();

        // Initialize angular state for precession/nutation tracking
        let mut angular_state = if self.inputs.enable_precession_nutation {
            Some(AngularState {
                pitch_angle: 0.001, // Small initial disturbance
                yaw_angle: 0.001,
                pitch_rate: 0.0,
                yaw_rate: 0.0,
                precession_angle: 0.0,
                nutation_phase: 0.0,
            })
        } else {
            None
        };
        let mut max_yaw_angle = 0.0;
        let mut max_precession_angle = 0.0;

        // Calculate air density
        let (air_density, speed_of_sound, resolved_temp_c, resolved_press_hpa) =
            self.resolved_atmosphere();
        // MBA-1136 (rank 30): base density RATIO for the local-altitude atmosphere recompute done
        // per-substep inside calculate_acceleration. The `air_density` / `speed_of_sound` above
        // stay the frozen station values, still used for the Mach-transition, pitch-damping and
        // precession/nutation diagnostics (which are intentionally referenced to station Mach).
        let base_ratio = air_density / 1.225;

        // Wind vector (McCoy): X=downrange (head/tail wind), Y=0, Z=lateral (crosswind)
        // 0deg = headwind, 90deg = from the right (McCoy wind-FROM convention, matching
        // WindSock); wind enters drag via velocity - wind. Used when no segmented wind.
        // MBA-728: no shear/no segments here, so vertical_speed passes straight through
        // (there is no horizontal-only scaling step on this path).
        let wind_vector =
            crate::wind::wind_vector(self.wind.speed, self.wind.direction, self.wind.vertical_speed);

        // Pitch-damping coefficients depend only on the (constant) bullet_model; compute once
        // instead of re-deriving them (with a to_lowercase alloc) every integration step.
        let pitch_coeffs = PitchDampingCoefficients::from_bullet_type(
            self.inputs.bullet_model.as_deref().unwrap_or("default"),
        );

        // Main integration loop (X is downrange)
        while position.x < self.max_range
            && position.y > self.inputs.ground_threshold
            && time < TRAJECTORY_TIME_LIMIT_S
        {
            // Store trajectory point
            let velocity_magnitude = velocity.magnitude();
            let kinetic_energy =
                0.5 * self.inputs.bullet_mass * velocity_magnitude * velocity_magnitude;

            self.push_trajectory_point(
                &mut points,
                TrajectoryPoint {
                    time,
                    position,
                    velocity_magnitude,
                    kinetic_energy,
                    drag_coefficient: None,
                },
            )?;

            // Record Mach-transition distances (constant sea-level speed of sound, matching the
            // transonic_mach tracking). Each threshold is recorded once, in descending order.
            {
                let mach_here = if speed_of_sound > 0.0 {
                    velocity_magnitude / speed_of_sound
                } else {
                    0.0
                };
                mach_transitions.record_downward_crossings(
                    mach_here,
                    position.x,
                    &mut transonic_distances,
                );
            }

            // Track max height
            if position.y > max_height {
                max_height = position.y;
            }

            // Calculate pitch damping if enabled
            if self.inputs.enable_pitch_damping {
                let mach = velocity_magnitude / speed_of_sound;

                // Track when we enter transonic
                if transonic_mach.is_none() && mach < 1.2 && mach > 0.8 {
                    transonic_mach = Some(mach);
                }

                // Calculate pitch damping coefficient
                let pitch_damping = calculate_pitch_damping_coefficient(mach, &pitch_coeffs);

                // Track minimum (most critical for stability)
                if pitch_damping < min_pitch_damping {
                    min_pitch_damping = pitch_damping;
                }
            }

            // Calculate precession/nutation if enabled
            if self.inputs.enable_precession_nutation {
                if let Some(ref mut state) = angular_state {
                    let velocity_magnitude = velocity.magnitude();
                    let params = self.precession_nutation_params(
                        velocity_magnitude,
                        air_density,
                        speed_of_sound,
                    );

                    // Update angular state
                    *state = calculate_combined_angular_motion(
                        &params,
                        state,
                        time,
                        self.time_step,
                        0.001, // Initial disturbance
                    );

                    // Track maximums
                    if state.yaw_angle.abs() > max_yaw_angle {
                        max_yaw_angle = state.yaw_angle.abs();
                    }
                    if state.precession_angle.abs() > max_precession_angle {
                        max_precession_angle = state.precession_angle.abs();
                    }
                }
            }

            // Use the same acceleration kernel as RK4/RK45 so all three solvers share ONE drag
            // model. solve_euler previously used a bespoke frontal-area drag (0.5*rho*Cd*A*v^2/m)
            // that IGNORED the ballistic coefficient entirely (diverging up to ~2.3x from the
            // BC-retardation RK4/RK45 path), and also omitted the Magnus/Coriolis terms.
            // calculate_acceleration applies BC-retardation drag, gravity, Coriolis, Magnus, wind
            // shear, and the zero-relative-velocity gravity-only guard.
            let acceleration = self.calculate_acceleration(
                &position,
                &velocity,
                &wind_vector,
                (resolved_temp_c, resolved_press_hpa, base_ratio),
            );

            // Update state
            velocity += acceleration * self.time_step;
            position += velocity * self.time_step;
            time += self.time_step;
            self.validate_integration_state(&position, &velocity, time)?;
        }

        let termination =
            self.append_terminal_endpoint(&mut points, position, velocity, time, &mut max_height)?;

        // Get final values
        // MBA-1423: fill the per-point Cd here, once, rather than inside each integrator loop.
        // Every solver family reaches this line, so none can leave the field empty while the
        // others populate it. Must precede the `last_point` borrow below.
        self.annotate_drag_coefficients(&mut points, speed_of_sound);

        let last_point = points.last().ok_or("No trajectory points generated")?;

        // Create trajectory sampling data if enabled (shared helper, MBA-1403)
        let sampled_points = self.build_sampled_points(
            &points,
            max_height,
            transonic_distances,
            &mach_transitions,
        )?;

        Ok(TrajectoryResult {
            max_range: last_point.position.x, // X is downrange
            max_height,
            time_of_flight: last_point.time,
            impact_velocity: last_point.velocity_magnitude,
            impact_energy: last_point.kinetic_energy,
            projectile_mass_kg: self.inputs.bullet_mass,
            line_of_sight_height_m: self.inputs.muzzle_height + self.inputs.sight_height,
            station_speed_of_sound_mps: speed_of_sound,
            termination,
            points,
            sampled_points,
            min_pitch_damping: if self.inputs.enable_pitch_damping {
                Some(min_pitch_damping)
            } else {
                None
            },
            transonic_mach,
            angular_state,
            max_yaw_angle: if self.inputs.enable_precession_nutation {
                Some(max_yaw_angle)
            } else {
                None
            },
            max_precession_angle: if self.inputs.enable_precession_nutation {
                Some(max_precession_angle)
            } else {
                None
            },
            aerodynamic_jump: aj_components,
            mach_1_2_distance_m: mach_transitions.mach_1_2_distance_m,
            mach_1_0_distance_m: mach_transitions.mach_1_0_distance_m,
            mach_0_9_distance_m: mach_transitions.mach_0_9_distance_m,
        })
    }

    fn solve_rk4(&self) -> Result<TrajectoryResult, BallisticsError> {
        // RK4 trajectory integration for better accuracy
        let mut time = 0.0;
        // Bullet starts at the BORE position, which is muzzle_height above ground
        // The sight is sight_height ABOVE the bore, so we don't add sight_height here
        // The sight_height affects the LOS calculation and zero angle, not the starting position
        // cant-adjusted via initial_position (MBA-1286)
        let mut position = self.initial_position();

        // Calculate initial velocity components with both elevation and azimuth
        // McCoy coordinate system: X=downrange, Y=vertical, Z=lateral (right)
        // Launch direction includes the aerodynamic-jump muzzle perturbation when enabled
        // (a no-op returning the bare muzzle/azimuth angles otherwise). MBA-959. Computed
        // once here and reused for the result so it isn't evaluated twice per solve.
        let aj_components = self.aerodynamic_jump_components();
        let (launch_elev, launch_azim) = self.launch_angles_from(aj_components.as_ref());
        let horizontal_velocity = self.inputs.muzzle_velocity * launch_elev.cos();
        let mut velocity = Vector3::new(
            horizontal_velocity * launch_azim.cos(), // X: downrange (forward)
            self.inputs.muzzle_velocity * launch_elev.sin(), // Y: vertical component
            horizontal_velocity * launch_azim.sin(), // Z: lateral (side deviation)
        );

        let mut points = Vec::new();
        let mut max_height = position.y;
        let mut min_pitch_damping = f64::INFINITY; // Track minimum pitch damping coefficient
        let mut transonic_mach = None; // Track when we enter transonic
                                       // Downrange distances where the projectile crosses Mach 1.2 (transonic) then Mach 1.0
                                       // (subsonic), so the sampled trajectory output can flag those transitions
                                       // (trajectory_sampling::add_trajectory_flags consumes this).
        let mut transonic_distances: Vec<f64> = Vec::new();
        let mut mach_transitions = MachTransitionTracker::default();

        // Initialize angular state for precession/nutation tracking
        let mut angular_state = if self.inputs.enable_precession_nutation {
            Some(AngularState {
                pitch_angle: 0.001, // Small initial disturbance
                yaw_angle: 0.001,
                pitch_rate: 0.0,
                yaw_rate: 0.0,
                precession_angle: 0.0,
                nutation_phase: 0.0,
            })
        } else {
            None
        };
        let mut max_yaw_angle = 0.0;
        let mut max_precession_angle = 0.0;

        // Calculate air density
        let (air_density, speed_of_sound, resolved_temp_c, resolved_press_hpa) =
            self.resolved_atmosphere();
        // MBA-1136 (rank 30): base density RATIO for the local-altitude atmosphere recompute done
        // per-substep inside calculate_acceleration. The `air_density` / `speed_of_sound` above
        // stay the frozen station values, still used for the Mach-transition, pitch-damping and
        // precession/nutation diagnostics (which are intentionally referenced to station Mach).
        let base_ratio = air_density / 1.225;

        // Wind vector (McCoy): X=downrange (head/tail wind), Y=0, Z=lateral (crosswind)
        // 0deg = headwind, 90deg = from the right (McCoy wind-FROM convention, matching
        // WindSock); wind enters drag via velocity - wind. Used when no segmented wind.
        // MBA-728: no shear/no segments here, so vertical_speed passes straight through
        // (there is no horizontal-only scaling step on this path).
        let wind_vector =
            crate::wind::wind_vector(self.wind.speed, self.wind.direction, self.wind.vertical_speed);

        // Pitch-damping coefficients depend only on the (constant) bullet_model; compute once
        // instead of re-deriving them (with a to_lowercase alloc) every integration step.
        let pitch_coeffs = PitchDampingCoefficients::from_bullet_type(
            self.inputs.bullet_model.as_deref().unwrap_or("default"),
        );

        // Main RK4 integration loop (X is downrange)
        while position.x < self.max_range
            && position.y > self.inputs.ground_threshold
            && time < TRAJECTORY_TIME_LIMIT_S
        {
            // Store trajectory point
            let velocity_magnitude = velocity.magnitude();
            let kinetic_energy =
                0.5 * self.inputs.bullet_mass * velocity_magnitude * velocity_magnitude;

            self.push_trajectory_point(
                &mut points,
                TrajectoryPoint {
                    time,
                    position,
                    velocity_magnitude,
                    kinetic_energy,
                    drag_coefficient: None,
                },
            )?;

            // Record Mach-transition distances (constant sea-level speed of sound, matching the
            // transonic_mach tracking). Each threshold is recorded once, in descending order.
            {
                let mach_here = if speed_of_sound > 0.0 {
                    velocity_magnitude / speed_of_sound
                } else {
                    0.0
                };
                mach_transitions.record_downward_crossings(
                    mach_here,
                    position.x,
                    &mut transonic_distances,
                );
            }

            if position.y > max_height {
                max_height = position.y;
            }

            // Calculate pitch damping if enabled (RK4 solver)
            if self.inputs.enable_pitch_damping {
                let mach = velocity_magnitude / speed_of_sound;

                // Track when we enter transonic
                if transonic_mach.is_none() && mach < 1.2 && mach > 0.8 {
                    transonic_mach = Some(mach);
                }

                // Calculate pitch damping coefficient
                let pitch_damping = calculate_pitch_damping_coefficient(mach, &pitch_coeffs);

                // Track minimum (most critical for stability)
                if pitch_damping < min_pitch_damping {
                    min_pitch_damping = pitch_damping;
                }
            }

            // Calculate precession/nutation if enabled (RK4 solver)
            if self.inputs.enable_precession_nutation {
                if let Some(ref mut state) = angular_state {
                    let velocity_magnitude = velocity.magnitude();
                    let params = self.precession_nutation_params(
                        velocity_magnitude,
                        air_density,
                        speed_of_sound,
                    );

                    // Update angular state
                    *state = calculate_combined_angular_motion(
                        &params,
                        state,
                        time,
                        self.time_step,
                        0.001, // Initial disturbance
                    );

                    // Track maximums
                    if state.yaw_angle.abs() > max_yaw_angle {
                        max_yaw_angle = state.yaw_angle.abs();
                    }
                    if state.precession_angle.abs() > max_precession_angle {
                        max_precession_angle = state.precession_angle.abs();
                    }
                }
            }

            // RK4 method
            let dt = self.time_step;

            // k1
            let acc1 = self.calculate_acceleration(
                &position,
                &velocity,
                &wind_vector,
                (resolved_temp_c, resolved_press_hpa, base_ratio),
            );

            // k2
            let pos2 = position + velocity * (dt * 0.5);
            let vel2 = velocity + acc1 * (dt * 0.5);
            let acc2 = self.calculate_acceleration(
                &pos2,
                &vel2,
                &wind_vector,
                (resolved_temp_c, resolved_press_hpa, base_ratio),
            );

            // k3
            let pos3 = position + vel2 * (dt * 0.5);
            let vel3 = velocity + acc2 * (dt * 0.5);
            let acc3 = self.calculate_acceleration(
                &pos3,
                &vel3,
                &wind_vector,
                (resolved_temp_c, resolved_press_hpa, base_ratio),
            );

            // k4
            let pos4 = position + vel3 * dt;
            let vel4 = velocity + acc3 * dt;
            let acc4 = self.calculate_acceleration(
                &pos4,
                &vel4,
                &wind_vector,
                (resolved_temp_c, resolved_press_hpa, base_ratio),
            );

            // Update position and velocity
            position += (velocity + vel2 * 2.0 + vel3 * 2.0 + vel4) * (dt / 6.0);
            velocity += (acc1 + acc2 * 2.0 + acc3 * 2.0 + acc4) * (dt / 6.0);
            time += dt;
            self.validate_integration_state(&position, &velocity, time)?;
        }

        let termination =
            self.append_terminal_endpoint(&mut points, position, velocity, time, &mut max_height)?;

        // Get final values
        // MBA-1423: fill the per-point Cd here, once, rather than inside each integrator loop.
        // Every solver family reaches this line, so none can leave the field empty while the
        // others populate it. Must precede the `last_point` borrow below.
        self.annotate_drag_coefficients(&mut points, speed_of_sound);

        let last_point = points.last().ok_or("No trajectory points generated")?;

        // Create trajectory sampling data if enabled (shared helper, MBA-1403)
        let sampled_points = self.build_sampled_points(
            &points,
            max_height,
            transonic_distances,
            &mach_transitions,
        )?;

        Ok(TrajectoryResult {
            max_range: last_point.position.x, // X is downrange
            max_height,
            time_of_flight: last_point.time,
            impact_velocity: last_point.velocity_magnitude,
            impact_energy: last_point.kinetic_energy,
            projectile_mass_kg: self.inputs.bullet_mass,
            line_of_sight_height_m: self.inputs.muzzle_height + self.inputs.sight_height,
            station_speed_of_sound_mps: speed_of_sound,
            termination,
            points,
            sampled_points,
            min_pitch_damping: if self.inputs.enable_pitch_damping {
                Some(min_pitch_damping)
            } else {
                None
            },
            transonic_mach,
            angular_state,
            max_yaw_angle: if self.inputs.enable_precession_nutation {
                Some(max_yaw_angle)
            } else {
                None
            },
            max_precession_angle: if self.inputs.enable_precession_nutation {
                Some(max_precession_angle)
            } else {
                None
            },
            aerodynamic_jump: aj_components,
            mach_1_2_distance_m: mach_transitions.mach_1_2_distance_m,
            mach_1_0_distance_m: mach_transitions.mach_1_0_distance_m,
            mach_0_9_distance_m: mach_transitions.mach_0_9_distance_m,
        })
    }

    fn solve_rk45(&self) -> Result<TrajectoryResult, BallisticsError> {
        // RK45 adaptive step size integration (Dormand-Prince method)
        let mut time = 0.0;
        // Bullet starts at the BORE position, which is muzzle_height above ground
        // The sight is sight_height ABOVE the bore, so we don't add sight_height here
        // cant-adjusted via initial_position (MBA-1286)
        let mut position = self.initial_position();

        // Calculate initial velocity components
        // McCoy coordinate system: X=downrange, Y=vertical, Z=lateral (right)
        // Launch direction includes the aerodynamic-jump muzzle perturbation when enabled
        // (a no-op returning the bare muzzle/azimuth angles otherwise). MBA-959. Computed
        // once here and reused for the result so it isn't evaluated twice per solve.
        let aj_components = self.aerodynamic_jump_components();
        let (launch_elev, launch_azim) = self.launch_angles_from(aj_components.as_ref());
        let horizontal_velocity = self.inputs.muzzle_velocity * launch_elev.cos();
        let mut velocity = Vector3::new(
            horizontal_velocity * launch_azim.cos(), // X: downrange (forward)
            self.inputs.muzzle_velocity * launch_elev.sin(), // Y: vertical component
            horizontal_velocity * launch_azim.sin(), // Z: lateral (side deviation)
        );

        let mut points = Vec::new();
        let mut max_height = position.y;
        let mut dt = 0.001; // Initial step size

        // Air density and wind are constant for the whole solve (self.atmosphere / self.wind
        // are immutable); compute once instead of every iteration (mirrors solve_rk4).
        let (air_density, speed_of_sound, resolved_temp_c, resolved_press_hpa) =
            self.resolved_atmosphere();
        // MBA-1136 (rank 30): base density RATIO for the local-altitude atmosphere recompute done
        // per-substep inside calculate_acceleration. The `air_density` / `speed_of_sound` above
        // stay the frozen station values, still used for the Mach-transition, pitch-damping and
        // precession/nutation diagnostics (which are intentionally referenced to station Mach).
        let base_ratio = air_density / 1.225;
        // 0deg = headwind, 90deg = from the right (McCoy wind-FROM convention, matching
        // WindSock); wind enters drag via velocity - wind. Used when no segmented wind.
        // MBA-728: no shear/no segments here, so vertical_speed passes straight through
        // (there is no horizontal-only scaling step on this path).
        let wind_vector =
            crate::wind::wind_vector(self.wind.speed, self.wind.direction, self.wind.vertical_speed);

        // Mach-transition distances for the sampled-output flags (see solve_euler/solve_rk4).
        let mut transonic_distances: Vec<f64> = Vec::new();
        let mut mach_transitions = MachTransitionTracker::default();

        // Pitch-damping / precession diagnostics (MBA-966). Previously only the
        // Euler and fixed-RK4 solvers tracked these, so the default adaptive
        // RK45 path always reported null even with --enable-pitch-damping /
        // --enable-precession set. Mirror the RK4 tracking here.
        let mut min_pitch_damping = f64::INFINITY;
        let mut transonic_mach: Option<f64> = None;
        let pitch_coeffs = PitchDampingCoefficients::from_bullet_type(
            self.inputs.bullet_model.as_deref().unwrap_or("default"),
        );
        let mut angular_state = if self.inputs.enable_precession_nutation {
            Some(AngularState {
                pitch_angle: 0.001,
                yaw_angle: 0.001,
                pitch_rate: 0.0,
                yaw_rate: 0.0,
                precession_angle: 0.0,
                nutation_phase: 0.0,
            })
        } else {
            None
        };
        let mut max_yaw_angle = 0.0;
        let mut max_precession_angle = 0.0;

        while position.x < self.max_range
            && position.y > self.inputs.ground_threshold
            && time < TRAJECTORY_TIME_LIMIT_S
        {
            // Store current point
            let velocity_magnitude = velocity.magnitude();
            let kinetic_energy = 0.5 * self.inputs.bullet_mass * velocity_magnitude.powi(2);

            self.push_trajectory_point(
                &mut points,
                TrajectoryPoint {
                    time,
                    position,
                    velocity_magnitude,
                    kinetic_energy,
                    drag_coefficient: None,
                },
            )?;

            // Record Mach-transition distances (constant sea-level speed of sound, matching the
            // transonic_mach tracking). Each threshold is recorded once, in descending order.
            {
                let mach_here = if speed_of_sound > 0.0 {
                    velocity_magnitude / speed_of_sound
                } else {
                    0.0
                };
                mach_transitions.record_downward_crossings(
                    mach_here,
                    position.x,
                    &mut transonic_distances,
                );
            }

            if position.y > max_height {
                max_height = position.y;
            }

            // Pitch damping (RK45 solver) — track the minimum coefficient and the
            // Mach at which the projectile enters the transonic band (MBA-966).
            if self.inputs.enable_pitch_damping {
                let mach = velocity_magnitude / speed_of_sound;
                if transonic_mach.is_none() && mach < 1.2 && mach > 0.8 {
                    transonic_mach = Some(mach);
                }
                let pitch_damping = calculate_pitch_damping_coefficient(mach, &pitch_coeffs);
                if pitch_damping < min_pitch_damping {
                    min_pitch_damping = pitch_damping;
                }
            }

            // Retry the same state until the embedded error estimate accepts the
            // candidate. No trajectory or angular state advances on rejection.
            let accepted_step = self.adaptive_rk45_step(
                &position,
                &velocity,
                dt,
                &wind_vector,
                (resolved_temp_c, resolved_press_hpa, base_ratio),
            );
            debug_assert!(
                accepted_step.error <= RK45_TOLERANCE || accepted_step.used_dt <= RK45_MIN_DT
            );

            // Precession / nutation advances only after the translational step
            // is accepted, using that accepted interval rather than a rejected
            // trial's dt.
            if self.inputs.enable_precession_nutation {
                if let Some(ref mut state) = angular_state {
                    let params = self.precession_nutation_params(
                        velocity_magnitude,
                        air_density,
                        speed_of_sound,
                    );

                    *state = calculate_combined_angular_motion(
                        &params,
                        state,
                        time,
                        accepted_step.used_dt,
                        0.001,
                    );

                    if state.yaw_angle.abs() > max_yaw_angle {
                        max_yaw_angle = state.yaw_angle.abs();
                    }
                    if state.precession_angle.abs() > max_precession_angle {
                        max_precession_angle = state.precession_angle.abs();
                    }
                }
            }

            position = accepted_step.position;
            velocity = accepted_step.velocity;
            time += accepted_step.used_dt;
            self.validate_integration_state(&position, &velocity, time)?;

            // Adapt the step size for the NEXT iteration.
            dt = accepted_step.next_dt;
        }

        // Ensure we have at least one point
        if points.is_empty() {
            return Err(BallisticsError::from("No trajectory points calculated"));
        }

        // Shared MBA-968/MBA-1218 range-crossing interpolation for all solver modes.
        let termination =
            self.append_terminal_endpoint(&mut points, position, velocity, time, &mut max_height)?;

        // MBA-1423: see the sibling solvers — fill the per-point Cd before the borrow below.
        self.annotate_drag_coefficients(&mut points, speed_of_sound);

        let last_point = points.last().unwrap();

        // Generate sampled trajectory points if enabled (shared helper, MBA-1403)
        let sampled_points = self.build_sampled_points(
            &points,
            max_height,
            transonic_distances,
            &mach_transitions,
        )?;

        Ok(TrajectoryResult {
            max_range: last_point.position.x, // X is downrange
            max_height,
            time_of_flight: last_point.time,
            impact_velocity: last_point.velocity_magnitude,
            impact_energy: last_point.kinetic_energy,
            projectile_mass_kg: self.inputs.bullet_mass,
            line_of_sight_height_m: self.inputs.muzzle_height + self.inputs.sight_height,
            station_speed_of_sound_mps: speed_of_sound,
            termination,
            points,
            sampled_points,
            min_pitch_damping: if self.inputs.enable_pitch_damping {
                Some(min_pitch_damping)
            } else {
                None
            },
            transonic_mach,
            angular_state,
            max_yaw_angle: if self.inputs.enable_precession_nutation {
                Some(max_yaw_angle)
            } else {
                None
            },
            max_precession_angle: if self.inputs.enable_precession_nutation {
                Some(max_precession_angle)
            } else {
                None
            },
            aerodynamic_jump: aj_components,
            mach_1_2_distance_m: mach_transitions.mach_1_2_distance_m,
            mach_1_0_distance_m: mach_transitions.mach_1_0_distance_m,
            mach_0_9_distance_m: mach_transitions.mach_0_9_distance_m,
        })
    }

    fn adaptive_rk45_step(
        &self,
        position: &Vector3<f64>,
        velocity: &Vector3<f64>,
        initial_dt: f64,
        wind_vector: &Vector3<f64>,
        resolved_atmo: (f64, f64, f64),
    ) -> Rk45AcceptedStep {
        let mut trial_dt = initial_dt;

        loop {
            let trial = self.rk45_step(
                position,
                velocity,
                trial_dt,
                wind_vector,
                RK45_TOLERANCE,
                resolved_atmo,
            );
            // A finite-but-extreme input or malformed optional curve can overflow an embedded
            // trial. Do not let a NaN suggested step poison `trial_dt` and retry forever: shrink
            // to the minimum step, return the non-finite trial there, and let the immediate
            // integration-state check turn it into a clean Err.
            let next_dt = if trial.suggested_dt.is_finite() {
                (RK45_SAFETY_FACTOR * trial.suggested_dt).clamp(RK45_MIN_DT, RK45_MAX_DT)
            } else {
                RK45_MIN_DT
            };

            if trial.error <= RK45_TOLERANCE || trial_dt <= RK45_MIN_DT {
                return Rk45AcceptedStep {
                    position: trial.position,
                    velocity: trial.velocity,
                    used_dt: trial_dt,
                    next_dt,
                    error: trial.error,
                };
            }

            trial_dt = next_dt;
        }
    }

    fn rk45_step(
        &self,
        position: &Vector3<f64>,
        velocity: &Vector3<f64>,
        dt: f64,
        wind_vector: &Vector3<f64>,
        tolerance: f64,
        resolved_atmo: (f64, f64, f64), // (base_temp_c, base_press_hpa, base_ratio)
    ) -> Rk45Trial {
        // Dormand-Prince coefficients
        const A21: f64 = 1.0 / 5.0;
        const A31: f64 = 3.0 / 40.0;
        const A32: f64 = 9.0 / 40.0;
        const A41: f64 = 44.0 / 45.0;
        const A42: f64 = -56.0 / 15.0;
        const A43: f64 = 32.0 / 9.0;
        const A51: f64 = 19372.0 / 6561.0;
        const A52: f64 = -25360.0 / 2187.0;
        const A53: f64 = 64448.0 / 6561.0;
        const A54: f64 = -212.0 / 729.0;
        const A61: f64 = 9017.0 / 3168.0;
        const A62: f64 = -355.0 / 33.0;
        const A63: f64 = 46732.0 / 5247.0;
        const A64: f64 = 49.0 / 176.0;
        const A65: f64 = -5103.0 / 18656.0;
        const A71: f64 = 35.0 / 384.0;
        const A73: f64 = 500.0 / 1113.0;
        const A74: f64 = 125.0 / 192.0;
        const A75: f64 = -2187.0 / 6784.0;
        const A76: f64 = 11.0 / 84.0;

        // 5th order coefficients
        const B1: f64 = 35.0 / 384.0;
        const B3: f64 = 500.0 / 1113.0;
        const B4: f64 = 125.0 / 192.0;
        const B5: f64 = -2187.0 / 6784.0;
        const B6: f64 = 11.0 / 84.0;

        // 4th order coefficients for error estimation
        const B1_ERR: f64 = 5179.0 / 57600.0;
        const B3_ERR: f64 = 7571.0 / 16695.0;
        const B4_ERR: f64 = 393.0 / 640.0;
        const B5_ERR: f64 = -92097.0 / 339200.0;
        const B6_ERR: f64 = 187.0 / 2100.0;
        const B7_ERR: f64 = 1.0 / 40.0;

        // Compute RK45 stages
        let k1_v = self.calculate_acceleration(position, velocity, wind_vector, resolved_atmo);
        let k1_p = *velocity;

        let p2 = position + dt * A21 * k1_p;
        let v2 = velocity + dt * A21 * k1_v;
        let k2_v = self.calculate_acceleration(&p2, &v2, wind_vector, resolved_atmo);
        let k2_p = v2;

        let p3 = position + dt * (A31 * k1_p + A32 * k2_p);
        let v3 = velocity + dt * (A31 * k1_v + A32 * k2_v);
        let k3_v = self.calculate_acceleration(&p3, &v3, wind_vector, resolved_atmo);
        let k3_p = v3;

        let p4 = position + dt * (A41 * k1_p + A42 * k2_p + A43 * k3_p);
        let v4 = velocity + dt * (A41 * k1_v + A42 * k2_v + A43 * k3_v);
        let k4_v = self.calculate_acceleration(&p4, &v4, wind_vector, resolved_atmo);
        let k4_p = v4;

        let p5 = position + dt * (A51 * k1_p + A52 * k2_p + A53 * k3_p + A54 * k4_p);
        let v5 = velocity + dt * (A51 * k1_v + A52 * k2_v + A53 * k3_v + A54 * k4_v);
        let k5_v = self.calculate_acceleration(&p5, &v5, wind_vector, resolved_atmo);
        let k5_p = v5;

        let p6 = position + dt * (A61 * k1_p + A62 * k2_p + A63 * k3_p + A64 * k4_p + A65 * k5_p);
        let v6 = velocity + dt * (A61 * k1_v + A62 * k2_v + A63 * k3_v + A64 * k4_v + A65 * k5_v);
        let k6_v = self.calculate_acceleration(&p6, &v6, wind_vector, resolved_atmo);
        let k6_p = v6;

        let p7 = position + dt * (A71 * k1_p + A73 * k3_p + A74 * k4_p + A75 * k5_p + A76 * k6_p);
        let v7 = velocity + dt * (A71 * k1_v + A73 * k3_v + A74 * k4_v + A75 * k5_v + A76 * k6_v);
        let k7_v = self.calculate_acceleration(&p7, &v7, wind_vector, resolved_atmo);
        let k7_p = v7;

        // 5th order solution
        let new_pos = position + dt * (B1 * k1_p + B3 * k3_p + B4 * k4_p + B5 * k5_p + B6 * k6_p);
        let new_vel = velocity + dt * (B1 * k1_v + B3 * k3_v + B4 * k4_v + B5 * k5_v + B6 * k6_v);

        // 4th order solution for error estimate
        let pos_err = position
            + dt * (B1_ERR * k1_p
                + B3_ERR * k3_p
                + B4_ERR * k4_p
                + B5_ERR * k5_p
                + B6_ERR * k6_p
                + B7_ERR * k7_p);
        let vel_err = velocity
            + dt * (B1_ERR * k1_v
                + B3_ERR * k3_v
                + B4_ERR * k4_v
                + B5_ERR * k5_v
                + B6_ERR * k6_v
                + B7_ERR * k7_v);

        // Estimate error
        let error = cli_rk45_error_norm(position, velocity, &new_pos, &new_vel, &pos_err, &vel_err);

        // Calculate new step size
        let dt_new = if error < tolerance {
            dt * (tolerance / error).powf(0.2).min(2.0)
        } else {
            dt * (tolerance / error).powf(0.25).max(0.1)
        };

        Rk45Trial {
            position: new_pos,
            velocity: new_vel,
            suggested_dt: dt_new,
            error,
        }
    }

    fn apply_cluster_bc_correction(&self, base_bc: f64, velocity_fps: f64) -> f64 {
        if let Some(ref cluster_bc) = self.cluster_bc {
            cluster_bc.apply_correction_for_drag_model(
                base_bc,
                self.inputs.caliber_inches,
                self.inputs.weight_grains,
                velocity_fps,
                self.inputs.bc_type,
            )
        } else {
            base_bc
        }
    }

    fn calculate_acceleration(
        &self,
        position: &Vector3<f64>,
        velocity: &Vector3<f64>,
        wind_vector: &Vector3<f64>,
        resolved_atmo: (f64, f64, f64), // (base_temp_c, base_press_hpa, base_ratio) hoisted per-solve
    ) -> Vector3<f64> {
        // Resolve the wind at this point. Downrange-segmented wind (when supplied)
        // takes precedence and is sampled by downrange distance (position.x) per
        // step; otherwise altitude-dependent shear (if enabled); otherwise the
        // constant `wind_vector`. Segmented wind is not combined with shear (the
        // CLI/WASM front-ends reject that combination), so the order is safe.
        let actual_wind = if let Some(ref sock) = self.wind_sock {
            sock.vector_for_range_stateless(position.x)
        } else if self.inputs.enable_wind_shear {
            self.get_wind_at_altitude(position.y)
        } else {
            *wind_vector
        };
        let actual_wind =
            crate::derivatives::level_vector_to_shot_frame(actual_wind, self.inputs.shooting_angle);

        let relative_velocity = velocity - actual_wind;
        let velocity_magnitude = relative_velocity.magnitude();

        if velocity_magnitude < 0.001 {
            return self.gravity_acceleration();
        }

        // MBA-1136 (rank 30): recompute the atmosphere at the LOCAL substep altitude instead of
        // holding the frozen station scalars for the whole flight. This mirrors what
        // derivatives.rs / fast_trajectory.rs already do, so all three solver families vary air
        // density AND speed of sound with altitude (matters on elevated / long-range shots; a
        // no-op at the shooter altitude, where the ratio-based density recovers the station value
        // exactly). base_* were resolved once per solve via resolved_atmosphere().
        //
        // `base_temp_c` / `base_press_hpa` are the STATION conditions that seed the local
        // atmosphere calculation below. Magnus dynamic stability consumes the resulting local
        // density rather than freezing a launch-density correction.
        let (base_temp_c, base_press_hpa, station_ratio) = resolved_atmo;

        // MBA-1137: downrange-segmented atmosphere. When an AtmoSock is present, swap the BASE
        // (station-referenced) T/P/H tuple for the active zone selected by downrange distance
        // (position.x), recomputing the per-zone base density ratio via CIPM. That swapped base
        // then flows through the SAME altitude-lapse pipeline, so downrange zone selection and the
        // world-vertical altitude lapse compose — the zone sets the base density/humidity, and the
        // lapse multiplies on top of it (no double-count). When None, this is the resolved
        // single-station base.
        let (drag_base_temp_c, drag_base_press_hpa, drag_base_ratio, drag_humidity_percent) =
            if let Some(ref sock) = self.atmo_sock {
                let (zone_temp_c, zone_press_hpa, zone_humidity) = sock.atmo_for_range(position.x);
                let zone_base_ratio = crate::atmosphere::calculate_air_density_cimp(
                    zone_temp_c,
                    zone_press_hpa,
                    zone_humidity,
                ) / 1.225;
                (zone_temp_c, zone_press_hpa, zone_base_ratio, zone_humidity)
            } else {
                (
                    base_temp_c,
                    base_press_hpa,
                    station_ratio,
                    self.atmosphere.humidity,
                )
            };
        let local_alt = crate::atmosphere::shot_frame_altitude(
            self.atmosphere.altitude,
            position.x,
            position.y,
            self.inputs.shooting_angle,
        );
        let (air_density, speed_of_sound) = crate::atmosphere::get_local_atmosphere_humid(
            local_alt,
            self.atmosphere.altitude,
            drag_base_temp_c,
            drag_base_press_hpa,
            drag_base_ratio,
            drag_humidity_percent,
        );

        // Resolve the Cd and the retardation denominator through the shared resolver. MBA-1423
        // reports a per-point Cd derived from these same two values, so the number a consumer
        // charts cannot drift from the number actually flown here.
        let (cd, retard_denom) = self.drag_terms(velocity_magnitude, speed_of_sound);

        // Convert velocity to fps (still used below outside the BC lookup).
        let velocity_fps = velocity_magnitude * 3.28084;

        // Use proper ballistics retardation formula
        // This matches the proven formula from fast_trajectory.rs
        // The standard retardation factor converts Cd to drag deceleration
        // Note: velocity_fps already calculated above for BC segment lookup
        let cd_to_retard = crate::constants::CD_TO_RETARD;
        let standard_factor = cd * cd_to_retard;
        let density_scale = air_density / 1.225; // Scale relative to standard air (1.225 kg/m³)

        // Drag acceleration in ft/s² then convert to m/s²
        let a_drag_ft_s2 =
            (velocity_fps * velocity_fps) * standard_factor * density_scale / retard_denom;
        let a_drag_m_s2 = a_drag_ft_s2 * 0.3048; // ft/s² to m/s²

        // Apply drag opposite to velocity direction
        let drag_acceleration = -a_drag_m_s2 * (relative_velocity / velocity_magnitude);

        // Total acceleration = drag + gravity. `shooting_angle` rotates gravity into the shot
        // frame for inclined fire; at 0 deg this is the normal vertical-only gravity vector.
        let mut accel = drag_acceleration + self.gravity_acceleration();

        // Coriolis (Earth rotation). McCoy frame: X=downrange, Y=vertical, Z=lateral,
        // azimuth 0 = North. McCoy frame: X=downrange, Y=vertical, Z=lateral.
        if self.inputs.enable_coriolis {
            if let Some(lat_deg) = self.inputs.latitude {
                let omega_earth = 7.2921159e-5_f64; // rad/s
                let lat = lat_deg.to_radians();
                let az = self.inputs.shot_azimuth; // compass bearing (0=N), NOT the aiming offset
                                                   // Earth's angular velocity in level downrange/up/lateral axes.
                                                   // Projecting Omega=(0, Ω cosφ, Ω sinφ) [local E,N,U] by azimuth gives
                                                   // a NEGATIVE lateral component:
                                                   // lateral = downrange × up points East for a North shot, and
                                                   // Omega·East = -Ω cosφ sin(az). The previous code dropped that sign.
                let omega = Vector3::new(
                    omega_earth * lat.cos() * az.cos(),  // X: downrange
                    omega_earth * lat.sin(),             // Y: vertical
                    -omega_earth * lat.cos() * az.sin(), // Z: lateral (MBA-938: corrected sign)
                );
                let omega = crate::derivatives::level_vector_to_shot_frame(
                    omega,
                    self.inputs.shooting_angle,
                );
                // Coriolis acceleration is the physical -2 Ω×v (MBA-938). The old +2 with
                // an "output-preserving relabel" justification produced left-ward drift for
                // a North shot in the Northern hemisphere; first principles (and the +Eötvös
                // lift for East shots) require -2 with the corrected omega above.
                accel += -2.0 * omega.cross(velocity);
            }
        }

        // Magnus force (spinning projectile). SI units in this solver.
        // MBA-1134 (rank 35): the canonical empirical Litz spin-drift post-process
        // (apply_spin_drift) already captures the gyroscopic yaw-of-repose lateral, so the
        // explicit Magnus side force must NOT be added on top of it — otherwise the two lateral
        // models stack and double-count the drift. Suppress Magnus whenever Litz spin drift is
        // active. (The inverse is intentionally NOT done: Litz is not suppressed when Magnus is on.)
        if self.inputs.enable_magnus
            && !self.inputs.use_enhanced_spin_drift
            && self.inputs.bullet_diameter > 0.0
            && self.inputs.twist_rate > 0.0
        {
            let diameter_m = self.inputs.bullet_diameter;
            let (spin_rad_s, spin_param) = crate::spin_drift::calculate_magnus_spin_state(
                self.inputs.muzzle_velocity,
                velocity_magnitude,
                self.inputs.twist_rate,
                diameter_m,
            );
            // Mach and dynamic stability both use the LOCAL atmosphere recomputed above.
            let mach = velocity_magnitude / speed_of_sound;

            // Imperial conversions for the stability / yaw-of-repose helpers.
            let d_in = self.inputs.bullet_diameter / 0.0254;
            let m_gr = self.inputs.bullet_mass / crate::constants::GRAINS_TO_KG;
            let l_in = if self.inputs.bullet_length > 0.0 {
                self.inputs.bullet_length / 0.0254
            } else {
                // MBA-1135: mass-based length estimate (was a mass-blind 4.5-caliber default).
                let est_m = crate::stability::estimate_bullet_length_m(
                    self.inputs.bullet_diameter,
                    self.inputs.bullet_mass,
                );
                if est_m > 0.0 {
                    est_m / 0.0254
                } else {
                    4.5 * d_in
                }
            };
            // Use current-flight Sg with the muzzle-set spin. The helper back-calculates the
            // effective twist from fixed spin and current airspeed, so Sg and yaw of repose grow
            // downrange instead of remaining tied to launch conditions.
            let sg = crate::spin_drift::calculate_dynamic_stability(
                m_gr,
                velocity_magnitude,
                spin_rad_s,
                d_in,
                l_in,
                air_density,
            );

            // Yaw of repose (radians); zero for unstable bullets (Sg <= 1).
            let (yaw_rad, _) = crate::spin_drift::calculate_yaw_of_repose(
                sg,
                velocity_magnitude,
                spin_rad_s,
                0.0, // crosswind handled elsewhere
                0.0, // pitch rate not tracked
                air_density,
                d_in,
                l_in,
                m_gr,
                mach,
                "match",
                false,
            );

            // Proper McCoy Magnus FORCE: F = q S C_Npa (pd/2V) sin(alpha_R).
            let c_np = crate::derivatives::calculate_magnus_moment_coefficient(mach);
            let area = std::f64::consts::PI * (diameter_m / 2.0).powi(2);
            let magnus_force = 0.5
                * air_density
                * velocity_magnitude.powi(2)
                * area
                * c_np
                * spin_param
                * yaw_rad.sin();

            // The yaw of repose is lateral, so its Magnus force follows gravity projected normal
            // to flight (down for right-hand twist). Lateral yaw lift belongs to the separate Litz
            // spin-drift model and must not be synthesized from this Magnus magnitude.
            if magnus_force.abs() > 1e-12 {
                if let Some(dir) = crate::derivatives::yaw_of_repose_magnus_direction(
                    relative_velocity,
                    self.gravity_acceleration(),
                    self.inputs.is_twist_right,
                ) {
                    accel += (magnus_force / self.inputs.bullet_mass) * dir;
                }
            }
        }

        accel
    }

    /// The reference drag coefficient and the retardation denominator this solver divides by,
    /// at one speed.
    ///
    /// This is the single place the two are resolved. `calculate_acceleration` divides one by
    /// the other to produce drag; [`Self::effective_drag_coefficient`] recombines them into the
    /// projectile's own Cd for reporting (MBA-1423). Keeping both callers on one resolver is
    /// deliberate: a second copy of the BC-precedence ladder below is exactly the kind of
    /// duplicate that drifts silently once someone edits one of them.
    ///
    /// The denominator is a BC (lb/in²) for a G-model and the sectional density for a custom
    /// drag table, because that table already supplies the projectile's actual Cd.
    fn drag_terms(&self, velocity_magnitude: f64, speed_of_sound: f64) -> (f64, f64) {
        let cd = self.calculate_drag_coefficient(velocity_magnitude, speed_of_sound);

        let velocity_fps = velocity_magnitude * 3.28084;

        // Match the other solver families' BC precedence: enabled velocity-keyed segments first,
        // then legacy Mach-keyed segments, then the scalar BC. `use_bc_segments` gates velocity
        // tables, while explicit Mach segments remain active when it is false; derivatives.rs and
        // the fast solver preserve that legacy contract for callers that provide a Mach table.
        let (base_bc, bc_from_segments) = if let Some(segments) = self
            .inputs
            .bc_segments_data
            .as_ref()
            .filter(|segments| self.inputs.use_bc_segments && !segments.is_empty())
        {
            // Find matching segment for current velocity.
            (
                crate::bc_estimation::velocity_segment_bc(
                    velocity_fps,
                    segments,
                    self.inputs.bc_value,
                ),
                true,
            )
        } else if let Some(segments) = self
            .inputs
            .bc_segments
            .as_ref()
            .filter(|segments| !segments.is_empty())
        {
            (
                crate::derivatives::interpolated_bc(
                    velocity_magnitude / speed_of_sound,
                    segments,
                    Some(&self.inputs),
                ),
                true,
            )
        } else {
            (self.inputs.bc_value, false)
        };

        // Segment tables already own the velocity-dependent BC shape. Stacking the empirical
        // cluster ladder on top would apply that shape twice (MBA-1175). Cluster correction is
        // therefore only a fallback for a scalar BC, regardless of which explicit segment
        // representation supplied the active value.
        let effective_bc = if bc_from_segments {
            base_bc
        } else {
            self.apply_cluster_bc_correction(base_bc, velocity_fps)
        };
        // The scalar BC is validated at the solve boundary. Retain a small denominator floor for
        // explicit segment tables, whose interpolated values are independent caller data.
        let effective_bc = effective_bc.max(1e-6);

        // When a custom drag table is active, calculate_drag_coefficient returned the
        // projectile's ACTUAL Cd, so the retardation denominator must be the sectional
        // density (lb/in²), not a BC: Cd_own / SD == Cd_ref / BC
        // (see BallisticInputs::custom_drag_denominator).
        let retard_denom = if self.inputs.custom_drag_table.is_some() {
            self.inputs.custom_drag_denominator(effective_bc)
        } else {
            effective_bc
        };

        (cd, retard_denom)
    }

    /// The projectile's own drag coefficient at one speed — the curve a shooter would compare
    /// against a measured CDM trace, not the reference-table Cd (MBA-1423).
    ///
    /// A G-model Cd describes the *standard* projectile, so reporting it directly would hand a
    /// consumer the same G7 curve no matter which bullet was flown. Scaling by the form factor
    /// `SD / BC` converts it to this projectile's actual Cd, from the identity that both
    /// descriptions must produce the same retardation:
    ///
    /// ```text
    /// Cd_own / SD == Cd_ref / BC   =>   Cd_own == Cd_ref * SD / BC
    /// ```
    ///
    /// Because the denominator comes from the same internal resolver the integrator divides by,
    /// a velocity- or Mach-segmented BC
    /// carries into the result and its band steps appear in the reported curve. The same
    /// expression is already correct for a custom drag table, where the denominator *is* the
    /// sectional density and the scale factor collapses to 1.
    ///
    /// Returns `None` when the projectile's mass or diameter is unavailable, since sectional
    /// density — and therefore the projectile's own Cd — is undefined without both.
    pub fn effective_drag_coefficient(
        &self,
        velocity_magnitude: f64,
        speed_of_sound: f64,
    ) -> Option<f64> {
        if !velocity_magnitude.is_finite() || speed_of_sound <= 1e-9 {
            return None;
        }
        let sectional_density = self.inputs.sectional_density_lb_in2()?;
        let (cd, retard_denom) = self.drag_terms(velocity_magnitude, speed_of_sound);
        if retard_denom <= 0.0 {
            return None;
        }
        let effective = cd * sectional_density / retard_denom;
        effective.is_finite().then_some(effective)
    }

    /// Fill every point's [`TrajectoryPoint::drag_coefficient`] in one pass (MBA-1423).
    ///
    /// Deliberately a pass over the finished vector rather than a value set where points are
    /// built: this solver constructs them in several places — terminal interpolation plus one
    /// loop per integrator family — and a pass covers a site added later for free.
    ///
    /// `speed_of_sound` must be the station value the result reports Mach against, so that a
    /// consumer plotting Cd against Mach gets a self-consistent pair. Feeding each step's local
    /// speed of sound instead would attribute a Cd to a Mach the document never shows.
    fn annotate_drag_coefficients(&self, points: &mut [TrajectoryPoint], speed_of_sound: f64) {
        for point in points.iter_mut() {
            point.drag_coefficient =
                self.effective_drag_coefficient(point.velocity_magnitude, speed_of_sound);
        }
    }

    fn calculate_drag_coefficient(&self, velocity: f64, speed_of_sound: f64) -> f64 {
        let mach = velocity / speed_of_sound;

        // MBA-940: a user-supplied custom drag table is the final Cd, used as-is — no G-model
        // lookup, no transonic shape correction, no form factor. The supplied curve already
        // encodes the projectile's true drag, so applying those would distort/double-count it.
        if let Some(ref table) = self.inputs.custom_drag_table {
            // MBA-1357: cd_scale is a single whole-curve drag multiplier applied here, at the
            // Cd lookup site. The Mach-keyed DSF table (truing_dsf.rs) is a SEPARATE, drop-only
            // post-processing correction applied to a solved TrajectoryResult's points after
            // integration finishes — it never touches this Cd/drag-force computation.
            return table.interpolate(mach) * self.inputs.cd_scale;
        }

        // A published/measured BC already contains the projectile form factor (BC = SD / i).
        // Multiplying reference Cd by a second name-derived factor double-counts shape.
        crate::drag::get_drag_coefficient(mach, &self.inputs.bc_type)
    }
}

// Monte Carlo parameters
#[derive(Debug, Clone)]
pub struct MonteCarloParams {
    pub num_simulations: usize,
    pub velocity_std_dev: f64,
    pub angle_std_dev: f64,
    pub bc_std_dev: f64,
    pub wind_speed_std_dev: f64,
    pub target_distance: Option<f64>,
    pub base_wind_speed: f64,
    pub base_wind_direction: f64,
    pub azimuth_std_dev: f64, // Horizontal aiming variation in radians
}

impl Default for MonteCarloParams {
    fn default() -> Self {
        Self {
            num_simulations: 1000,
            velocity_std_dev: 1.0,
            angle_std_dev: 0.001,
            bc_std_dev: 0.01,
            wind_speed_std_dev: 1.0,
            target_distance: None,
            base_wind_speed: 0.0,
            base_wind_direction: 0.0,
            azimuth_std_dev: 0.001, // Default horizontal spread ~0.057 degrees
        }
    }
}

// Monte Carlo results
#[derive(Debug, Clone)]
pub struct MonteCarloResults {
    pub ranges: Vec<f64>,
    pub impact_velocities: Vec<f64>,
    /// Deviations from the baseline point of aim at the target plane.
    ///
    /// A sample that falls short of the plane is encoded as
    /// `(0, TARGET_NOT_REACHED_SENTINEL_M, 0)` so it remains aligned with
    /// `ranges` and `impact_velocities` and still counts as a miss.
    pub impact_positions: Vec<Vector3<f64>>,
}

/// Default hit-zone radius (meters) around the point of aim at the target plane — a 30 cm
/// circle. Shared by the CLI, FFI, and WASM so "hit probability" means the same thing everywhere.
pub const DEFAULT_HIT_RADIUS_M: f64 = 0.3;

/// Vertical-position marker for a Monte Carlo sample that never reached the target plane.
///
/// The marker preserves the equal-length result-vector and C-ABI contract. Exclude marked
/// positions from target-plane dispersion statistics, but keep them in the denominator for hit
/// probability because they are definite misses.
pub const TARGET_NOT_REACHED_SENTINEL_M: f64 = -1.0e9;

impl MonteCarloResults {
    /// Whether an encoded impact position represents a finite arrival at the target plane.
    pub fn position_reached_target(position: &Vector3<f64>) -> bool {
        position.iter().all(|component| component.is_finite())
            && position.y != TARGET_NOT_REACHED_SENTINEL_M
    }

    /// Number of recorded simulations that reached the target plane.
    pub fn target_arrival_count(&self) -> usize {
        self.impact_positions
            .iter()
            .filter(|position| Self::position_reached_target(position))
            .count()
    }

    /// Fraction of recorded simulations that fell short of (or otherwise failed to produce a
    /// finite position at) the target plane.
    pub fn target_shortfall_fraction(&self) -> f64 {
        if self.impact_positions.is_empty() {
            return 0.0;
        }
        (self.impact_positions.len() - self.target_arrival_count()) as f64
            / self.impact_positions.len() as f64
    }

    /// Upper-median radial miss among samples that reached the target plane.
    ///
    /// This preserves the CLI's historical radial-to-baseline "CEP (approx)" convention while
    /// preventing the finite target-shortfall marker from becoming the median (MBA-1159).
    /// Returns `None` when no recorded simulation reached the target plane.
    pub fn target_plane_cep(&self) -> Option<f64> {
        let mut radial_misses: Vec<f64> = self
            .impact_positions
            .iter()
            .filter(|position| Self::position_reached_target(position))
            .map(Vector3::norm)
            .filter(|miss| miss.is_finite())
            .collect();
        radial_misses.sort_by(f64::total_cmp);
        if radial_misses.is_empty() {
            None
        } else {
            Some(radial_misses[radial_misses.len() / 2])
        }
    }

    /// Fraction of simulations whose impact at the target plane lands within `hit_radius_m`
    /// of the point of aim. `impact_positions` are deviations from the baseline at the target
    /// plane (the downrange component is 0), so the vector norm is the radial miss distance.
    /// Samples that fall short of the target remain in the denominator and count as misses.
    /// Returns 0.0 when there are no samples.
    ///
    /// Single source of truth for hit probability — previously the CLI used a range-precision
    /// notion and the FFI a position notion with a redundant clause, so they disagreed.
    pub fn hit_probability(&self, hit_radius_m: f64) -> f64 {
        if self.impact_positions.is_empty() {
            return 0.0;
        }
        let hits = self
            .impact_positions
            .iter()
            .filter(|position| {
                Self::position_reached_target(position) && position.norm() < hit_radius_m
            })
            .count();
        hits as f64 / self.impact_positions.len() as f64
    }

    /// Fraction of simulations whose impact at the target plane lands within the axis-aligned
    /// rectangle `width_m` (lateral, Z) x `height_m` (vertical, Y) centered on the point of aim
    /// — i.e. within `width_m / 2` of center laterally AND `height_m / 2` of center vertically.
    ///
    /// This is the WEZ (Weapon Employment Zone, MBA-1317) counterpart of [`hit_probability`]'s
    /// circular radius for rectangular target sizes (e.g. an 18"x30" steel plate). Samples that
    /// fall short of the target plane remain in the denominator and count as misses, matching
    /// `hit_probability`'s convention. Returns 0.0 when there are no samples or when either
    /// dimension is non-positive.
    ///
    /// [`hit_probability`]: Self::hit_probability
    pub fn rect_hit_probability(&self, width_m: f64, height_m: f64) -> f64 {
        let dimensions_invalid = width_m.is_nan()
            || width_m <= 0.0
            || height_m.is_nan()
            || height_m <= 0.0;
        if self.impact_positions.is_empty() || dimensions_invalid {
            return 0.0;
        }
        let half_width = width_m / 2.0;
        let half_height = height_m / 2.0;
        let hits = self
            .impact_positions
            .iter()
            .filter(|position| {
                Self::position_reached_target(position)
                    && position.z.abs() <= half_width
                    && position.y.abs() <= half_height
            })
            .count();
        hits as f64 / self.impact_positions.len() as f64
    }
}

fn wind_from_signed_speed_sample(
    signed_speed: f64,
    sampled_direction: f64,
    vertical_speed: f64,
) -> WindConditions {
    // The base wind's vertical component rides along un-dispersed: vertical wind is a
    // systematic input (MBA-728), not a sampled dispersion source. Dropping it here
    // would make every per-sample solve disagree with the baseline solve by the whole
    // vertical deflection — a phantom bias in the MC statistics.
    if signed_speed < 0.0 {
        WindConditions {
            speed: -signed_speed,
            direction: sampled_direction + std::f64::consts::PI,
            vertical_speed,
        }
    } else {
        WindConditions {
            speed: signed_speed,
            direction: sampled_direction,
            vertical_speed,
        }
    }
}

struct MonteCarloWindSampler {
    speed: rand_distr::Normal<f64>,
    direction: rand_distr::Normal<f64>,
    /// Base wind's vertical component, carried into every sample un-dispersed.
    vertical_speed: f64,
}

impl MonteCarloWindSampler {
    fn new(
        base_wind: &WindConditions,
        wind_speed_std_dev: f64,
        wind_direction_std_dev: f64,
    ) -> Result<Self, BallisticsError> {
        use rand_distr::Normal;

        if !wind_direction_std_dev.is_finite() || wind_direction_std_dev < 0.0 {
            return Err("Wind direction standard deviation must be finite and non-negative".into());
        }

        let speed = Normal::new(base_wind.speed, wind_speed_std_dev)
            .map_err(|e| format!("Invalid wind speed distribution: {e}"))?;
        let direction = Normal::new(base_wind.direction, wind_direction_std_dev)
            .map_err(|e| format!("Invalid wind direction distribution: {e}"))?;
        Ok(Self { speed, direction, vertical_speed: base_wind.vertical_speed })
    }

    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> WindConditions {
        use rand_distr::Distribution;

        wind_from_signed_speed_sample(
            self.speed.sample(rng),
            self.direction.sample(rng),
            self.vertical_speed,
        )
    }
}

// Run Monte Carlo simulation (backwards compatibility)
pub fn run_monte_carlo(
    base_inputs: BallisticInputs,
    params: MonteCarloParams,
) -> Result<MonteCarloResults, BallisticsError> {
    run_monte_carlo_with_direction_std_dev(base_inputs, params, 0.0)
}

/// Run Monte Carlo with an independent wind-direction standard deviation in radians.
///
/// The older [`run_monte_carlo`] entry point remains source compatible and delegates here with
/// zero direction uncertainty.
pub fn run_monte_carlo_with_direction_std_dev(
    base_inputs: BallisticInputs,
    params: MonteCarloParams,
    wind_direction_std_dev: f64,
) -> Result<MonteCarloResults, BallisticsError> {
    let base_wind = WindConditions {
        speed: params.base_wind_speed,
        direction: params.base_wind_direction,
        vertical_speed: 0.0,
    };
    run_monte_carlo_with_wind_and_direction_std_dev(
        base_inputs,
        base_wind,
        params,
        wind_direction_std_dev,
    )
}

// Run Monte Carlo simulation with wind
pub fn run_monte_carlo_with_wind(
    base_inputs: BallisticInputs,
    base_wind: WindConditions,
    params: MonteCarloParams,
) -> Result<MonteCarloResults, BallisticsError> {
    run_monte_carlo_with_wind_and_direction_std_dev(base_inputs, base_wind, params, 0.0)
}

/// Run Monte Carlo with explicit base wind and independent direction uncertainty in radians.
///
/// The older [`run_monte_carlo_with_wind`] entry point delegates here with zero direction
/// uncertainty, preserving its API while removing the former speed-to-angle unit conflation.
pub fn run_monte_carlo_with_wind_and_direction_std_dev(
    base_inputs: BallisticInputs,
    base_wind: WindConditions,
    params: MonteCarloParams,
    wind_direction_std_dev: f64,
) -> Result<MonteCarloResults, BallisticsError> {
    let mut rng = rand::rng();
    run_monte_carlo_with_wind_and_direction_std_dev_using_rng(
        base_inputs,
        base_wind,
        params,
        wind_direction_std_dev,
        &mut rng,
    )
}

/// Run Monte Carlo with an explicit PRNG seed, for deterministic/reproducible output.
///
/// Otherwise identical to [`run_monte_carlo_with_wind_and_direction_std_dev`], which draws
/// from the process-global RNG instead. Intended for callers that need repeatable results
/// across runs — e.g. tests, or a WEZ (Weapon Employment Zone, MBA-1317) sweep that a caller
/// wants to reproduce exactly while tuning target size or wind-call error.
pub fn run_monte_carlo_with_wind_and_direction_std_dev_seeded(
    base_inputs: BallisticInputs,
    base_wind: WindConditions,
    params: MonteCarloParams,
    wind_direction_std_dev: f64,
    seed: u64,
) -> Result<MonteCarloResults, BallisticsError> {
    use rand::{rngs::StdRng, SeedableRng};
    let mut rng = StdRng::seed_from_u64(seed);
    run_monte_carlo_with_wind_and_direction_std_dev_using_rng(
        base_inputs,
        base_wind,
        params,
        wind_direction_std_dev,
        &mut rng,
    )
}

fn run_monte_carlo_with_wind_and_direction_std_dev_using_rng<R: rand::Rng + ?Sized>(
    base_inputs: BallisticInputs,
    base_wind: WindConditions,
    params: MonteCarloParams,
    wind_direction_std_dev: f64,
    rng: &mut R,
) -> Result<MonteCarloResults, BallisticsError> {
    use rand_distr::{Distribution, Normal};

    let mut ranges = Vec::new();
    let mut impact_velocities = Vec::new();
    let mut impact_positions = Vec::new();

    let atmosphere = AtmosphericConditions {
        temperature: base_inputs.temperature,
        pressure: base_inputs.pressure,
        humidity: base_inputs.humidity_percent(),
        altitude: base_inputs.altitude,
    };
    let target_hint = params
        .target_distance
        .unwrap_or(base_inputs.target_distance);
    let solver_max_range = target_hint.max(1000.0) * 2.0;

    // First, calculate baseline trajectory with no variations
    let mut baseline_solver =
        TrajectorySolver::new(base_inputs.clone(), base_wind.clone(), atmosphere.clone());
    baseline_solver.set_max_range(solver_max_range);
    let baseline_result = baseline_solver.solve()?;

    // Determine target distance: use explicit target or baseline max range
    let target_distance = params.target_distance.unwrap_or(baseline_result.max_range);

    // Get baseline position at target distance (interpolated)
    let baseline_at_target = baseline_result
        .position_at_range(target_distance)
        .ok_or("Could not interpolate baseline at target distance")?;

    // Create normal distributions for variations
    // Sample muzzle velocity as a DELTA and apply it after TrajectorySolver::new resolves the
    // powder-temperature model. Sampling an absolute value here let a powder curve overwrite
    // every draw in the constructor, collapsing the requested dispersion to zero (MBA-1176).
    let velocity_delta_dist = Normal::new(0.0, params.velocity_std_dev)
        .map_err(|e| format!("Invalid velocity distribution: {}", e))?;
    let angle_dist = Normal::new(base_inputs.muzzle_angle, params.angle_std_dev)
        .map_err(|e| format!("Invalid angle distribution: {}", e))?;
    let bc_dist = Normal::new(base_inputs.bc_value, params.bc_std_dev)
        .map_err(|e| format!("Invalid BC distribution: {}", e))?;
    // Direction uncertainty is an independent angular quantity in radians. Do not derive it from
    // wind-speed uncertainty: meters/second cannot supply an angular standard deviation.
    let wind_sampler = MonteCarloWindSampler::new(
        &base_wind,
        params.wind_speed_std_dev,
        wind_direction_std_dev,
    )?;
    let azimuth_dist = Normal::new(base_inputs.azimuth_angle, params.azimuth_std_dev)
        .map_err(|e| format!("Invalid azimuth distribution: {}", e))?;

    for _ in 0..params.num_simulations {
        // Create varied inputs
        let mut inputs = base_inputs.clone();
        let muzzle_velocity_delta = velocity_delta_dist.sample(&mut *rng);
        inputs.muzzle_angle = angle_dist.sample(&mut *rng);
        inputs.bc_value = bc_dist.sample(&mut *rng).max(0.01);
        inputs.azimuth_angle = azimuth_dist.sample(&mut *rng); // Add horizontal variation

        // Create varied wind (now based on base wind conditions)
        let wind = wind_sampler.sample(&mut *rng);

        // Run trajectory
        let mut solver = TrajectorySolver::new(inputs, wind, atmosphere.clone());
        solver.inputs.muzzle_velocity =
            (solver.inputs.muzzle_velocity + muzzle_velocity_delta).max(0.0);
        solver.set_max_range(solver_max_range);
        match solver.solve() {
            Ok(result) => {
                // MBA-967: do NOT skip samples that fall short of the target. range/velocity are
                // recorded at GROUND IMPACT for EVERY sample, so "Mean Range" is the ground-impact
                // distribution — independent of target_distance and consistent with `trajectory`.
                // All three result vectors still grow together per sample, so the equal-length FFI
                // ABI (exposed under one count) is preserved.
                let deviation = if result.max_range < target_distance {
                    // This sample never reached the target plane -> definite miss. Keep the
                    // encoded miss finite but far outside any practical target radius.
                    Vector3::new(0.0, TARGET_NOT_REACHED_SENTINEL_M, 0.0)
                } else {
                    let pos_at_target = match result.position_at_range(target_distance) {
                        Some(p) => p,
                        None => continue, // defensive: skip the whole sample (keeps vectors aligned)
                    };
                    // Deviation from baseline at the SAME target distance (McCoy): X = downrange
                    // (0 here), Y = vertical (elevation), Z = lateral (windage). Muzzle-angle
                    // sampling already models vertical pointing dispersion, so do not add a
                    // second independent vertical pointing draw here.
                    Vector3::new(
                        0.0,
                        pos_at_target.y - baseline_at_target.y,
                        pos_at_target.z - baseline_at_target.z,
                    )
                };

                ranges.push(result.max_range);
                impact_velocities.push(result.impact_velocity);
                impact_positions.push(deviation);
            }
            Err(_) => {
                // Skip failed simulations
                continue;
            }
        }
    }

    if ranges.is_empty() {
        return Err("No successful simulations".into());
    }

    Ok(MonteCarloResults {
        ranges,
        impact_velocities,
        impact_positions,
    })
}

// Calculate zero angle for a target
pub fn calculate_zero_angle(
    inputs: BallisticInputs,
    target_distance: f64,
    target_height: f64,
) -> Result<f64, BallisticsError> {
    calculate_zero_angle_with_conditions(
        inputs,
        target_distance,
        target_height,
        WindConditions::default(),
        AtmosphericConditions::default(),
    )
}

pub fn calculate_zero_angle_with_conditions(
    inputs: BallisticInputs,
    target_distance: f64,
    target_height: f64,
    wind: WindConditions,
    atmosphere: AtmosphericConditions,
) -> Result<f64, BallisticsError> {
    let mut solver = TrajectorySolver::new(inputs, wind, atmosphere);
    solver.calculate_and_set_zero_angle(target_distance, target_height, ZeroTargetFrame::SightLine)
}

/// [`calculate_zero_angle_with_conditions`] for a presence-aware caller that has already
/// resolved station temperature/pressure (MBA-1397; e.g. reduced a declared QNH via
/// [`crate::atmosphere::resolve_station_conditions_with_pressure_mode`]). Unlike the base
/// function, `atmosphere`'s temperature/pressure are trusted as-is and never re-interpreted
/// through the legacy default-sentinel heuristic.
pub fn calculate_zero_angle_with_resolved_conditions(
    inputs: BallisticInputs,
    target_distance: f64,
    target_height: f64,
    wind: WindConditions,
    atmosphere: AtmosphericConditions,
) -> Result<f64, BallisticsError> {
    let mut solver = TrajectorySolver::new_with_resolved_station_atmosphere(inputs, wind, atmosphere);
    solver.calculate_and_set_zero_angle(target_distance, target_height, ZeroTargetFrame::SightLine)
}

/// Generous solve envelope for [`calculate_zero_range_from_angle_with_conditions`] /
/// [`calculate_zero_range_from_angle_with_resolved_conditions`] (MBA-1402): a "zero" bore angle
/// is a near-horizontal shot, so its far line-of-sight crossing sits well inside this even at
/// extended small-arms ranges. Double the engine's own default 1000 m solve envelope
/// ([`TrajectorySolver::new`]) rather than something scaled to an (unknown, being solved for)
/// target distance; `solve()` still stops at ground impact well before this in the typical
/// case, so the larger envelope costs effectively nothing.
pub const ZERO_RANGE_FROM_ANGLE_MAX_RANGE_M: f64 = 2000.0;

/// Solve the zero RANGE(S) that a fixed bore angle produces. Runs the trajectory at
/// `zero_angle_rad` and returns BOTH line-of-sight crossings it finds, as
/// [`ZeroCrossings`] — a rifle sighted above the bore generally crosses a level line of
/// sight twice on a rising shot (near, ascending, close to the muzzle; far, descending past
/// the apex).
///
/// **This is NOT the exact inverse of [`calculate_zero_angle_with_conditions`].** A single
/// bore angle generally implies two valid zero distances (the classic 25/300-yard
/// battle-zero relationship is exactly this: one angle, two zeros), and
/// [`calculate_zero_angle_with_conditions`]'s own choice of root depends on which target
/// distance it was asked to solve for. Round-tripping a solved angle back through this
/// function recovers the ORIGINAL distance as one of the two returned crossings, not
/// necessarily as `far_m` specifically — for short/flat zeros the forward solver's answer is
/// often the NEAR crossing.
pub fn calculate_zero_range_from_angle_with_conditions(
    inputs: BallisticInputs,
    zero_angle_rad: f64,
    target_height: f64,
    wind: WindConditions,
    atmosphere: AtmosphericConditions,
) -> Result<ZeroCrossings, BallisticsError> {
    let mut solver = TrajectorySolver::new(inputs, wind, atmosphere);
    solver.set_max_range(ZERO_RANGE_FROM_ANGLE_MAX_RANGE_M);
    solver.find_zero_range(zero_angle_rad, target_height, ZeroTargetFrame::SightLine)
}

/// [`calculate_zero_range_from_angle_with_conditions`] for a presence-aware caller that has
/// already resolved station temperature/pressure — same relationship as
/// [`calculate_zero_angle_with_resolved_conditions`] has to
/// [`calculate_zero_angle_with_conditions`] (MBA-1397 pattern). See that function's doc
/// comment for why this returns both crossings rather than a single "the" zero range.
pub fn calculate_zero_range_from_angle_with_resolved_conditions(
    inputs: BallisticInputs,
    zero_angle_rad: f64,
    target_height: f64,
    wind: WindConditions,
    atmosphere: AtmosphericConditions,
) -> Result<ZeroCrossings, BallisticsError> {
    let mut solver = TrajectorySolver::new_with_resolved_station_atmosphere(inputs, wind, atmosphere);
    solver.set_max_range(ZERO_RANGE_FROM_ANGLE_MAX_RANGE_M);
    solver.find_zero_range(zero_angle_rad, target_height, ZeroTargetFrame::SightLine)
}

/// What a BC estimate is fit against.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BcFitMode {
    /// Data points are `(distance_m, drop_m)` — the classic drop-curve fit.
    Drop,
    /// Data points are `(distance_m, velocity_mps)` — a velocity-retention fit,
    /// which is immune to zero / sight-height / launch-angle error.
    Velocity,
}

/// The result of a single BC fit (one drag model, one fit basis).
#[derive(Debug, Clone, Copy)]
pub struct BcEstimate {
    /// The estimated ballistic coefficient.
    pub bc: f64,
    /// RMS residual across the data points, in fit units (meters of drop, or m/s of speed).
    pub rms_error: f64,
    /// Which standard drag model this BC is referenced to.
    pub drag_model: DragModel,
    /// Whether the fit was against drop or velocity data.
    pub mode: BcFitMode,
    /// True if the best fit landed at the edge of the physical BC search range — i.e. the
    /// data did not pin down an interior optimum (too sparse/short-range, or wrong units /
    /// atmosphere / zero). The reported `bc` is then a floor/ceiling, not a real estimate.
    pub at_bound: bool,
}

/// Interpolate the fitted quantity (drop in meters, or speed in m/s) at a downrange
/// distance from a solved trajectory. `None` if the trajectory never reaches `target_dist`.
///
/// `drop_offset` is subtracted-from convention: for `Drop` the returned value is
/// `drop_offset - y`. With `drop_offset = 0` this is bore-referenced drop (flat fire);
/// with `drop_offset = sight_height` and a zeroed trajectory it is drop below the
/// (horizontal) line of sight — i.e. dope-card drop.
fn fit_value_at(
    points: &[TrajectoryPoint],
    target_dist: f64,
    mode: BcFitMode,
    drop_offset: f64,
) -> Option<f64> {
    let val = |p: &TrajectoryPoint| match mode {
        BcFitMode::Drop => drop_offset - p.position.y,
        BcFitMode::Velocity => p.velocity_magnitude,
    };
    for i in 0..points.len() {
        if points[i].position.x >= target_dist {
            if i == 0 {
                return Some(val(&points[0]));
            }
            let p1 = &points[i - 1];
            let p2 = &points[i];
            let dx = p2.position.x - p1.position.x;
            if dx.abs() < 1e-9 {
                return Some(val(p2));
            }
            let t = (target_dist - p1.position.x) / dx;
            return Some(val(p1) + t * (val(p2) - val(p1)));
        }
    }
    None
}

fn fit_residual_sse(
    trajectory: &[TrajectoryPoint],
    observations: &[(f64, f64)],
    mode: BcFitMode,
    drop_offset: f64,
) -> Option<f64> {
    if observations.is_empty() {
        return None;
    }
    let mut total = 0.0;
    for (target_dist, target_val) in observations {
        // Scores are comparable only when every candidate contains every residual term.
        // Reject a trajectory that terminates before even one observation (MBA-1178).
        let value = fit_value_at(trajectory, *target_dist, mode, drop_offset)?;
        let error = value - target_val;
        total += error * error;
    }
    Some(total)
}

/// Estimate a BC by fitting a simulated trajectory to measured data, for a chosen drag
/// model (G1, G7, …) and fit basis (drop or velocity). Uses a coarse 0.01 sweep over
/// plausible BCs followed by a 0.001 local refine around the coarse best.
///
/// `points` are `(distance_m, value_m_or_mps)` where the second element is drop in meters
/// (`BcFitMode::Drop`) or remaining speed in m/s (`BcFitMode::Velocity`).
///
/// The fit runs under `atmosphere` — BC is only meaningful relative to the air density the
/// data was measured at, so this must match the conditions the drop/velocity came from
/// (pass ICAO standard for a standard-atmosphere dope card).
///
/// `zero_range` selects the drop reference frame (ignored for velocity fits):
/// - `None` → **bore-referenced**: flat 0° fire, drop below the extended bore axis.
/// - `Some(range_m)` → **sight/dope-card-referenced**: the trajectory is zeroed at
///   `range_m` (using `sight_height`), and drop is measured below the horizontal line of
///   sight — i.e. exactly what a dope card zeroed at that range prints.
#[allow(clippy::too_many_arguments)] // Public compatibility API; grouping would be breaking.
pub fn estimate_bc_fit(
    velocity: f64,
    mass: f64,
    diameter: f64,
    points: &[(f64, f64)],
    drag_model: DragModel,
    mode: BcFitMode,
    atmosphere: AtmosphericConditions,
    zero_range: Option<f64>,
    sight_height: f64,
) -> Result<BcEstimate, BallisticsError> {
    if points.is_empty() {
        return Err(BallisticsError::from(
            "No data points provided for BC estimation.".to_string(),
        ));
    }
    let max_dist = points.iter().map(|(d, _)| *d).fold(0.0_f64, f64::max);
    // For a zeroed drop fit, drop is below the horizontal LOS which sits `sight_height`
    // above the bore at the muzzle: drop = sight_height - y. Bore-referenced fits use 0.
    let drop_offset = if zero_range.is_some() { sight_height } else { 0.0 };

    // Sum of squared residuals for a trial BC; None unless the solve reaches ALL data points.
    let sse = |bc_value: f64| -> Option<f64> {
        let mut inputs = BallisticInputs {
            muzzle_velocity: velocity,
            bc_value,
            bc_type: drag_model,
            bullet_mass: mass,
            bullet_diameter: diameter,
            sight_height,
            ..Default::default()
        };
        // Zeroed fit: tilt the bore so the bullet crosses LOS at the zero range, so the
        // downrange drops match a dope card zeroed there. Bore fit leaves muzzle_angle = 0.
        if let Some(zr) = zero_range {
            // MBA-1130: zero to the LINE OF SIGHT (y = sight_height) at the zero range,
            // not the bore line (y = 0). Drop is measured as `drop_offset - y` with
            // drop_offset = sight_height, so a bore-referenced zero left drop != 0 at the
            // zero range and the drop-fit no longer round-tripped to the true BC. This
            // matches how range-table / come-up / dope-card generation zero.
            let za = calculate_zero_angle_with_conditions(
                inputs.clone(),
                zr,
                sight_height,
                WindConditions::default(),
                atmosphere.clone(),
            )
            .ok()?;
            inputs.muzzle_angle = za;
        }
        let mut solver =
            TrajectorySolver::new(inputs, WindConditions::default(), atmosphere.clone());
        solver.set_max_range(max_dist * 1.5);
        let result = solver.solve().ok()?;
        fit_residual_sse(&result.points, points, mode, drop_offset)
    };

    // Physical BC search range, per drag model. Real G7 BCs top out well under 0.5 (0.7 is
    // a generous ceiling); G1 BCs run higher. Keeping G7 out of G1 territory means a fit
    // that runs to the ceiling reports a sane bound, not a nonsensical 1.2.
    let (bc_min, bc_max) = match drag_model {
        DragModel::G7 => (0.05, 0.70),
        _ => (0.10, 1.20),
    };

    // Coarse sweep across the physical range.
    let mut best_bc = f64::NAN;
    let mut best_sse = f64::MAX;
    let mut bc = bc_min;
    while bc <= bc_max + 1e-9 {
        if let Some(s) = sse(bc) {
            if s < best_sse {
                best_sse = s;
                best_bc = bc;
            }
        }
        bc += 0.01;
    }
    if !best_bc.is_finite() {
        return Err(BallisticsError::from(
            "Unable to estimate BC from provided data. Check that the values and units are correct."
                .to_string(),
        ));
    }

    // Local refine at 0.001 resolution around the coarse best (kept within the range).
    let lo = (best_bc - 0.01).max(bc_min);
    let hi = (best_bc + 0.01).min(bc_max);
    let mut bc = lo;
    while bc <= hi + 1e-9 {
        if let Some(s) = sse(bc) {
            if s < best_sse {
                best_sse = s;
                best_bc = bc;
            }
        }
        bc += 0.001;
    }

    // A solution sitting on the search boundary means the data didn't determine an interior
    // optimum — the fit ran to the floor/ceiling. Flag it so callers don't trust the number.
    let at_bound = best_bc <= bc_min + 0.011 || best_bc >= bc_max - 0.011;
    // fit_residual_sse rejects partial trajectories, so best_sse contains exactly one residual
    // per input point and this denominator is also the honest matched-point count.
    let rms_error = (best_sse / points.len() as f64).sqrt();
    Ok(BcEstimate {
        bc: best_bc,
        rms_error,
        drag_model,
        mode,
        at_bound,
    })
}

/// Estimate a G1 BC from a drop curve. Back-compatible wrapper over [`estimate_bc_fit`];
/// `points` are `(distance_m, drop_m)`.
pub fn estimate_bc_from_trajectory(
    velocity: f64,
    mass: f64,
    diameter: f64,
    points: &[(f64, f64)], // (distance, drop) pairs
) -> Result<f64, BallisticsError> {
    estimate_bc_fit(
        velocity,
        mass,
        diameter,
        points,
        DragModel::G1,
        BcFitMode::Drop,
        AtmosphericConditions::default(),
        None,
        0.05,
    )
    .map(|e| e.bc)
}

// Add rand dependencies for Monte Carlo
use rand;
use rand_distr;

#[cfg(test)]
mod mba737_powder_resolution_tests {
    use super::*;

    #[test]
    fn linear_model_cold_powder_subtracts() {
        // 0.5486 m/s per degC (1 fps/degF), shot day 10 C below the 21.1 C reference.
        let v = resolve_powder_adjusted_velocity(823.0, 11.1, true, 0.5486, 21.1, None, None);
        assert!((v - (823.0 + 0.5486 * (11.1 - 21.1))).abs() < 1e-12);
        assert!(v < 823.0);
    }

    #[test]
    fn linear_model_hot_powder_adds() {
        let v = resolve_powder_adjusted_velocity(823.0, 31.1, true, 0.5486, 21.1, None, None);
        assert!((v - (823.0 + 0.5486 * 10.0)).abs() < 1e-12);
    }

    #[test]
    fn disabled_flag_is_passthrough() {
        let v = resolve_powder_adjusted_velocity(823.0, 40.0, false, 0.5486, 21.1, None, None);
        assert_eq!(v, 823.0);
    }

    #[test]
    fn curve_overrides_linear_and_interpolates_at_powder_temp() {
        let curve = [(4.4, 798.6), (21.1, 823.0), (37.8, 841.2)];
        // Explicit powder temp decouples from ambient: interpolate at 4.4 C, not 30 C.
        let v = resolve_powder_adjusted_velocity(823.0, 30.0, true, 99.0, 21.1, Some(&curve), Some(4.4));
        assert!((v - 798.6).abs() < 1e-9);
    }

    #[test]
    fn curve_falls_back_to_ambient_and_clamps() {
        let curve = [(4.4, 798.6), (37.8, 841.2)];
        // Ambient far below the coldest measured point: clamp, no extrapolation.
        let v = resolve_powder_adjusted_velocity(823.0, -40.0, true, 1.0, 21.1, Some(&curve), None);
        assert!((v - 798.6).abs() < 1e-9);
        let v_hot = resolve_powder_adjusted_velocity(823.0, 60.0, true, 1.0, 21.1, Some(&curve), None);
        assert!((v_hot - 841.2).abs() < 1e-9);
    }

    #[test]
    fn empty_curve_suppresses_linear_fallback() {
        // Historical `else if` semantics: Some-but-empty curve means NO adjustment.
        let v = resolve_powder_adjusted_velocity(823.0, 40.0, true, 0.5486, 21.1, Some(&[]), None);
        assert_eq!(v, 823.0);
    }

    #[test]
    fn sweep_huge_range_errors_instead_of_overflowing() {
        // Review finding: the old count math saturated the usize cast and the +1
        // panicked (debug) or wrapped to zero rows (release), bypassing the cap.
        assert!(parse_powder_sweep("0:1e20:1").is_err());
        assert!(parse_powder_sweep("0:1e308:1e-3").is_err());
    }

    #[test]
    fn sweep_fractional_step_keeps_end_row() {
        // 0.3/0.1 floats to 2.9999999999999996; the END row must survive the floor.
        let rows = parse_powder_sweep("0:0.3:0.1").unwrap();
        assert_eq!(rows.len(), 4);
        assert!((rows[3] - 0.3).abs() < 1e-9);
    }

    #[test]
    fn solver_and_helper_agree_on_linear_model() {
        // The solve() seam must fly exactly what the helper reports (MBA-737 contract).
        let inputs = BallisticInputs {
            use_powder_sensitivity: true,
            powder_temp_sensitivity: 0.5486,
            powder_temp: 21.1,
            temperature: 4.4,
            ..Default::default()
        };
        let expected = resolve_powder_adjusted_velocity(
            inputs.muzzle_velocity,
            inputs.temperature,
            true,
            0.5486,
            21.1,
            None,
            None,
        );
        let solver = TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        assert!((solver.inputs.muzzle_velocity - expected).abs() < 1e-12);
    }
}

#[cfg(test)]
mod mba1302_solver_seam_tests {
    use super::*;
    use crate::wind::WindSegment;

    #[test]
    fn authoritative_station_atmosphere_preserves_explicit_standard_values_at_altitude() {
        let atmosphere = AtmosphericConditions {
            temperature: 15.0,
            pressure: 1013.25,
            humidity: 50.0,
            altitude: 2_000.0,
        };
        let legacy = TrajectorySolver::new(
            BallisticInputs::default(),
            WindConditions::default(),
            atmosphere.clone(),
        );
        let authoritative = TrajectorySolver::new_with_resolved_station_atmosphere(
            BallisticInputs::default(),
            WindConditions::default(),
            atmosphere,
        );

        let (legacy_density, _, legacy_temp_c, legacy_pressure_hpa) = legacy.resolved_atmosphere();
        let (authoritative_density, _, authoritative_temp_c, authoritative_pressure_hpa) =
            authoritative.resolved_atmosphere();
        let (icao_temp_k, icao_pressure_pa) =
            crate::atmosphere::calculate_icao_standard_atmosphere(2_000.0);
        let (expected_authoritative_density, _) =
            crate::atmosphere::calculate_atmosphere(2_000.0, Some(15.0), Some(1013.25), 50.0);

        assert!((legacy_temp_c - (icao_temp_k - 273.15)).abs() < 1e-12);
        assert!((legacy_pressure_hpa - icao_pressure_pa / 100.0).abs() < 1e-12);
        assert_eq!(authoritative_temp_c.to_bits(), 15.0_f64.to_bits());
        assert_eq!(authoritative_pressure_hpa.to_bits(), 1013.25_f64.to_bits());
        assert_eq!(
            authoritative_density.to_bits(),
            expected_authoritative_density.to_bits()
        );
        assert!(
            (authoritative_density - legacy_density).abs() > 0.1,
            "explicit standard values at altitude must differ from ICAO-at-altitude: explicit={authoritative_density}, ICAO={legacy_density}"
        );
    }

    /// MBA-1397: the CLI/WASM `--pressure-type` mechanism precomputes station conditions via
    /// `atmosphere::resolve_station_conditions_with_pressure_mode` and constructs the solver
    /// with `new_with_resolved_station_atmosphere` UNCONDITIONALLY (not only for `Qnh`), to
    /// avoid re-deriving through the legacy sentinel a second time. This must be bit-for-bit
    /// equivalent to the historical `TrajectorySolver::new` (LegacyDefaultSentinels) path for
    /// every NON-ambiguous input -- i.e. for `PressureReferenceMode::Absolute`, precomputing
    /// and using the Authoritative constructor must reproduce `TrajectorySolver::new` exactly,
    /// which is the bit-level guarantee the CLI/WASM byte-identical-output tests rely on.
    #[test]
    fn precomputed_absolute_resolution_via_authoritative_matches_legacy_new() {
        for (temperature, pressure, altitude) in [
            (15.0, 1013.25, 0.0),   // sea-level default
            (15.0, 1013.25, 2000.0), // sentinel: omitted-pressure-at-altitude
            (-5.0, 850.0, 2000.0),  // explicit non-default station values
            (22.0, 950.0, 500.0),
        ] {
            let atmosphere = AtmosphericConditions {
                temperature,
                pressure,
                humidity: 50.0,
                altitude,
            };
            let legacy = TrajectorySolver::new(
                BallisticInputs::default(),
                WindConditions::default(),
                atmosphere.clone(),
            );

            let (resolved_temp_c, resolved_pressure_hpa) =
                crate::atmosphere::resolve_station_conditions_with_pressure_mode(
                    temperature,
                    pressure,
                    altitude,
                    crate::atmosphere::PressureReferenceMode::Absolute,
                );
            let precomputed_atmosphere = AtmosphericConditions {
                temperature: resolved_temp_c,
                pressure: resolved_pressure_hpa,
                humidity: 50.0,
                altitude,
            };
            let precomputed = TrajectorySolver::new_with_resolved_station_atmosphere(
                BallisticInputs::default(),
                WindConditions::default(),
                precomputed_atmosphere,
            );

            let (legacy_density, legacy_sos, legacy_temp_c, legacy_pressure_hpa) =
                legacy.resolved_atmosphere();
            let (pre_density, pre_sos, pre_temp_c, pre_pressure_hpa) =
                precomputed.resolved_atmosphere();

            assert_eq!(
                legacy_temp_c.to_bits(),
                pre_temp_c.to_bits(),
                "temperature=({temperature}, {pressure}, {altitude})"
            );
            assert_eq!(
                legacy_pressure_hpa.to_bits(),
                pre_pressure_hpa.to_bits(),
                "pressure=({temperature}, {pressure}, {altitude})"
            );
            assert_eq!(legacy_density.to_bits(), pre_density.to_bits());
            assert_eq!(legacy_sos.to_bits(), pre_sos.to_bits());
        }
    }

    fn configured_euler_zero(vertical_wind_mps: f64, time_step_s: f64) -> TrajectorySolver {
        let inputs = BallisticInputs {
            muzzle_velocity: 800.0,
            bc_value: 0.5,
            bc_type: DragModel::G7,
            bullet_mass: 0.0109,
            bullet_diameter: 0.00782,
            bullet_length: 0.0309,
            sight_height: 0.05,
            ground_threshold: -100.0,
            use_rk4: false,
            use_adaptive_rk45: false,
            ..BallisticInputs::default()
        };
        let mut solver = TrajectorySolver::new_with_resolved_station_atmosphere(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        solver.set_max_range(300.0);
        solver.set_time_step(time_step_s);
        if vertical_wind_mps != 0.0 {
            solver.set_wind_segments(vec![WindSegment {
                speed_kmh: 0.0,
                angle_deg: 0.0,
                until_m: 400.0,
                vertical_mps: vertical_wind_mps,
            }]);
        }
        solver
    }

    #[test]
    fn inclined_shot_zeroes_like_a_level_rifle() {
        // MBA-1412: since MBA-1302, zero_trial_height_at converted trial heights into the
        // gravity frame (shot_frame_altitude adds d*sin(shooting_angle) ~ 9 m at 5.71deg /
        // 100 yd) while find_zero_angle compared against a sight-frame target height, so ANY
        // inclined shot failed to bracket. A zero is a property of a level rifle's sight
        // geometry (same doctrine as MBA-1286's cant handling): the trial must solve level.
        const ZERO_DISTANCE_M: f64 = 91.44; // 100 yd
        const SIGHT_HEIGHT_M: f64 = 0.0381; // 1.5 in

        let inputs = BallisticInputs {
            bc_value: 0.5,
            bullet_mass: 150.0 * 0.06479891 / 1000.0,
            muzzle_velocity: 2700.0 * 0.3048,
            sight_height: SIGHT_HEIGHT_M,
            ..Default::default()
        };

        let mut level = inputs.clone();
        level.shooting_angle = 0.0;
        let level_angle = TrajectorySolver::new(level, Default::default(), Default::default())
            .find_zero_angle(ZERO_DISTANCE_M, SIGHT_HEIGHT_M, ZeroTargetFrame::SightLine)
            .expect("level zero must solve");

        let mut inclined = inputs;
        inclined.shooting_angle = 5.71_f64.to_radians();
        let inclined_angle =
            TrajectorySolver::new(inclined, Default::default(), Default::default())
                .find_zero_angle(ZERO_DISTANCE_M, SIGHT_HEIGHT_M, ZeroTargetFrame::SightLine)
                .expect("MBA-1412: a 5.71 deg incline at a 100 yd zero must be solvable");

        assert!(
            (inclined_angle - level_angle).abs() < 1e-9,
            "zeroing is level-rifle sight geometry; incline must not move the solved zero: \
             level={level_angle}, inclined={inclined_angle}"
        );
    }

    #[test]
    fn configured_zero_keeps_segments_method_and_time_step_then_sets_base_angle() {
        const TARGET_DISTANCE_M: f64 = 150.0;
        const TARGET_HEIGHT_M: f64 = 0.05;

        // A deliberately coarse Euler step makes an accidental reset to the historical 1 ms
        // zeroing step observable, while remaining stable and physically meaningful.
        let mut segmented = configured_euler_zero(-10.0, 0.02);
        let coarse_height = segmented
            .zero_trial_height_at(0.0, TARGET_DISTANCE_M, ZeroTargetFrame::SightLine)
            .expect("coarse configured trial")
            .expect("coarse trial reaches target");
        let mut fine = segmented.clone();
        fine.set_time_step(0.001);
        let fine_height = fine
            .zero_trial_height_at(0.0, TARGET_DISTANCE_M, ZeroTargetFrame::SightLine)
            .expect("fine configured trial")
            .expect("fine trial reaches target");
        assert!(
            (coarse_height - fine_height).abs() > 1e-5,
            "configured Euler step must affect zero trials: coarse={coarse_height}, fine={fine_height}"
        );

        let segmented_angle = segmented
            .calculate_and_set_zero_angle(TARGET_DISTANCE_M, TARGET_HEIGHT_M, ZeroTargetFrame::SightLine)
            .expect("segmented zero");
        assert_eq!(
            segmented.inputs.muzzle_angle.to_bits(),
            segmented_angle.to_bits(),
            "successful zero must install its angle on the configured solver"
        );
        assert_eq!(segmented.time_step.to_bits(), 0.02_f64.to_bits());
        assert_eq!(segmented.max_range.to_bits(), 300.0_f64.to_bits());
        assert!(segmented.wind_sock.is_some());
        assert_eq!(
            segmented.station_atmosphere_resolution,
            StationAtmosphereResolution::Authoritative
        );
        let zero_height = segmented
            .zero_trial_height_at(segmented_angle, TARGET_DISTANCE_M, ZeroTargetFrame::SightLine)
            .expect("verify segmented zero")
            .expect("zeroed trial reaches target");
        assert!(
            (zero_height - TARGET_HEIGHT_M).abs() < 0.0001,
            "configured zero missed target: height={zero_height}"
        );

        let mut calm = configured_euler_zero(0.0, 0.02);
        let calm_angle = calm
            .calculate_and_set_zero_angle(TARGET_DISTANCE_M, TARGET_HEIGHT_M, ZeroTargetFrame::SightLine)
            .expect("calm zero");
        assert!(
            (segmented_angle - calm_angle).abs() > 1e-5,
            "segmented vertical wind must participate in zero trials: segmented={segmented_angle}, calm={calm_angle}"
        );
    }
}

#[cfg(test)]
mod result_sanity_tests {
    use super::*;

    fn default_solver() -> TrajectorySolver {
        TrajectorySolver::new(
            BallisticInputs::default(),
            WindConditions::default(),
            AtmosphericConditions::default(),
        )
    }

    fn minimal_result() -> TrajectoryResult {
        TrajectoryResult {
            max_range: 100.0,
            max_height: 1.0,
            time_of_flight: 0.5,
            impact_velocity: 700.0,
            impact_energy: 2450.0,
            projectile_mass_kg: 0.01,
            line_of_sight_height_m: 1.5,
            station_speed_of_sound_mps: 340.0,
            termination: TrajectoryTermination::MaxRange,
            points: vec![],
            sampled_points: None,
            min_pitch_damping: None,
            transonic_mach: None,
            angular_state: None,
            max_yaw_angle: None,
            max_precession_angle: None,
            aerodynamic_jump: None,
            mach_1_2_distance_m: None,
            mach_1_0_distance_m: None,
            mach_0_9_distance_m: None,
        }
    }

    #[test]
    fn mba1293_negative_scalars_fail_the_result_postcondition() {
        let solver = default_solver();
        solver
            .validate_result_sanity(&minimal_result())
            .expect("a sane result must pass");

        for (name, mutate) in [
            ("max_range", (|r| r.max_range = -50.588) as fn(&mut TrajectoryResult)),
            ("time_of_flight", |r| r.time_of_flight = -1.0),
            ("impact_velocity", |r| r.impact_velocity = -700.0),
            ("impact_energy", |r| r.impact_energy = -1.0),
        ] {
            let mut result = minimal_result();
            mutate(&mut result);
            let error = solver
                .validate_result_sanity(&result)
                .expect_err("negative scalar must fail");
            assert!(
                error.to_string().contains(name),
                "error for {name} did not name the field: {error}"
            );
        }
    }

    #[test]
    fn mba1293_speed_budget_bounds_legitimate_states_and_rejects_divergence() {
        let solver = default_solver();
        let mv = solver.inputs.muzzle_velocity;

        // A state at muzzle speed is always inside the budget.
        let position = Vector3::new(10.0, 0.0, 0.0);
        solver
            .validate_integration_state(&position, &Vector3::new(mv, 0.0, 0.0), 0.01)
            .expect("muzzle-speed state must pass");

        // The MBA-1293 explosion (13x the muzzle speed) must be rejected as divergence.
        let error = solver
            .validate_integration_state(&position, &Vector3::new(-13.0 * mv, 0.0, 0.0), 0.01)
            .expect_err("13x muzzle speed must fail the budget");
        assert!(error.to_string().contains("diverged"), "{error}");

        // The budget grows with gravity's g*t so long lobbed arcs never trip it.
        let after_fall = mv + crate::constants::G_ACCEL_MPS2 * 60.0;
        solver
            .validate_integration_state(&position, &Vector3::new(0.0, -after_fall, 0.0), 60.0)
            .expect("gravity-accelerated speed within g*t must pass");
    }
}

#[cfg(test)]
mod trajectory_point_budget_tests {
    use super::*;
    use crate::MAX_TRAJECTORY_SAMPLES;

    fn solver_with_budget(
        use_rk4: bool,
        use_adaptive_rk45: bool,
        point_budget: usize,
        max_range: f64,
    ) -> TrajectorySolver {
        let inputs = BallisticInputs {
            use_rk4,
            use_adaptive_rk45,
            ground_threshold: f64::NEG_INFINITY,
            ..BallisticInputs::default()
        };
        let mut solver = TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        solver.max_trajectory_points = point_budget;
        solver.set_max_range(max_range);
        solver.set_time_step(0.001);
        solver
    }

    #[test]
    fn mba1283_every_solver_errors_instead_of_exceeding_point_budget() {
        for (mode, use_rk4, use_adaptive_rk45) in [
            ("Euler", false, false),
            ("RK4", true, false),
            ("RK45", true, true),
        ] {
            let error = solver_with_budget(use_rk4, use_adaptive_rk45, 3, 10.0)
                .solve()
                .expect_err("a solve requiring more than three points must fail");
            assert!(
                error.to_string().contains("point limit of 3"),
                "unexpected {mode} point-budget error: {error}"
            );
        }
    }

    #[test]
    fn mba1283_interpolated_endpoint_counts_toward_point_budget() {
        for (mode, use_rk4, use_adaptive_rk45) in [
            ("Euler", false, false),
            ("RK4", true, false),
            ("RK45", true, true),
        ] {
            let result = solver_with_budget(use_rk4, use_adaptive_rk45, 2, 0.1)
                .solve()
                .expect("the initial point plus exact endpoint fit a two-point budget");
            assert_eq!(result.points.len(), 2, "unexpected {mode} point count");

            let error = solver_with_budget(use_rk4, use_adaptive_rk45, 1, 0.1)
                .solve()
                .expect_err("the exact endpoint must not exceed a one-point budget");
            assert!(
                error.to_string().contains("point limit of 1"),
                "unexpected {mode} endpoint-budget error: {error}"
            );
        }
    }

    #[test]
    fn mba1299_every_solver_preflights_the_sample_budget() {
        for (mode, use_rk4, use_adaptive_rk45) in [
            ("Euler", false, false),
            ("RK4", true, false),
            ("RK45", true, true),
        ] {
            let inputs = BallisticInputs {
                use_rk4,
                use_adaptive_rk45,
                enable_trajectory_sampling: true,
                sample_interval: 1.0,
                ground_threshold: f64::NEG_INFINITY,
                ..BallisticInputs::default()
            };
            let mut solver = TrajectorySolver::new(
                inputs,
                WindConditions::default(),
                AtmosphericConditions::default(),
            );
            solver.set_max_range(MAX_TRAJECTORY_SAMPLES as f64);
            // If validation does not reject the sample grid before dispatch, the first attempted
            // integration point produces a distinct point-budget error.
            solver.max_trajectory_points = 0;

            let error = solver
                .solve()
                .expect_err("an over-limit sample grid must fail before integration");
            assert!(
                error
                    .to_string()
                    .contains("trajectory sample limit of 250000 exceeded"),
                "unexpected {mode} sample-budget error: {error}"
            );
        }
    }

    #[test]
    fn mba1299_normal_sampling_does_not_change_solver_results() {
        for (mode, use_rk4, use_adaptive_rk45) in [
            ("Euler", false, false),
            ("RK4", true, false),
            ("RK45", true, true),
        ] {
            let solve = |enable_trajectory_sampling| {
                let inputs = BallisticInputs {
                    use_rk4,
                    use_adaptive_rk45,
                    enable_trajectory_sampling,
                    sample_interval: 0.5,
                    ground_threshold: f64::NEG_INFINITY,
                    ..BallisticInputs::default()
                };
                let mut solver = TrajectorySolver::new(
                    inputs,
                    WindConditions::default(),
                    AtmosphericConditions::default(),
                );
                solver.set_max_range(2.0);
                solver.solve().expect("normal short-range solve")
            };

            let baseline = solve(false);
            let sampled = solve(true);
            for (field, left, right) in [
                ("max_range", baseline.max_range, sampled.max_range),
                ("max_height", baseline.max_height, sampled.max_height),
                (
                    "time_of_flight",
                    baseline.time_of_flight,
                    sampled.time_of_flight,
                ),
                (
                    "impact_velocity",
                    baseline.impact_velocity,
                    sampled.impact_velocity,
                ),
                (
                    "impact_energy",
                    baseline.impact_energy,
                    sampled.impact_energy,
                ),
            ] {
                assert_eq!(
                    left.to_bits(),
                    right.to_bits(),
                    "{mode} sampling changed {field}"
                );
            }
            assert_eq!(baseline.points.len(), sampled.points.len());
            for (index, (left, right)) in baseline
                .points
                .iter()
                .zip(&sampled.points)
                .enumerate()
            {
                assert_eq!(left.time.to_bits(), right.time.to_bits(), "{mode} point {index}");
                assert_eq!(
                    left.position.map(f64::to_bits),
                    right.position.map(f64::to_bits),
                    "{mode} point {index} position"
                );
                assert_eq!(
                    left.velocity_magnitude.to_bits(),
                    right.velocity_magnitude.to_bits(),
                    "{mode} point {index} velocity"
                );
                assert_eq!(
                    left.kinetic_energy.to_bits(),
                    right.kinetic_energy.to_bits(),
                    "{mode} point {index} energy"
                );
            }
            assert!(baseline.sampled_points.is_none());
            let samples = sampled
                .sampled_points
                .expect("sampling-enabled solve should return observations");
            assert_eq!(
                samples
                    .iter()
                    .map(|sample| sample.distance_m)
                    .collect::<Vec<_>>(),
                vec![0.0, 0.5, 1.0, 1.5, 2.0],
                "{mode} normal sampling grid changed"
            );
        }
    }
}

#[cfg(test)]
mod monte_carlo_result_tests {
    use super::*;

    fn make_results(impact_positions: Vec<Vector3<f64>>) -> MonteCarloResults {
        let count = impact_positions.len();
        MonteCarloResults {
            ranges: vec![500.0; count],
            impact_velocities: vec![300.0; count],
            impact_positions,
        }
    }

    #[test]
    fn target_plane_cep_excludes_shortfall_markers() {
        let mut positions: Vec<Vector3<f64>> = (1..=5)
            .map(|radius| Vector3::new(0.0, radius as f64, 0.0))
            .collect();
        positions.extend(
            (0..5).map(|_| Vector3::new(0.0, TARGET_NOT_REACHED_SENTINEL_M, 0.0)),
        );
        let results = make_results(positions);

        assert_eq!(results.target_arrival_count(), 5);
        assert_eq!(results.target_shortfall_fraction(), 0.5);
        assert_eq!(results.target_plane_cep(), Some(3.0));

        let one_shortfall = make_results(vec![
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 2.0, 0.0),
            Vector3::new(0.0, 3.0, 0.0),
            Vector3::new(0.0, 4.0, 0.0),
            Vector3::new(0.0, 5.0, 0.0),
            Vector3::new(0.0, TARGET_NOT_REACHED_SENTINEL_M, 0.0),
        ]);
        assert_eq!(one_shortfall.target_plane_cep(), Some(3.0));
    }

    #[test]
    fn all_shortfalls_have_no_cep_but_still_count_as_misses() {
        let all_shortfalls = make_results(vec![
            Vector3::new(0.0, TARGET_NOT_REACHED_SENTINEL_M, 0.0),
            Vector3::new(0.0, TARGET_NOT_REACHED_SENTINEL_M, 0.0),
        ]);
        assert_eq!(all_shortfalls.target_arrival_count(), 0);
        assert_eq!(all_shortfalls.target_shortfall_fraction(), 1.0);
        assert_eq!(all_shortfalls.target_plane_cep(), None);
        assert_eq!(all_shortfalls.hit_probability(0.3), 0.0);

        let one_hit_one_shortfall = make_results(vec![
            Vector3::new(0.0, 0.1, 0.0),
            Vector3::new(0.0, TARGET_NOT_REACHED_SENTINEL_M, 0.0),
        ]);
        assert_eq!(one_hit_one_shortfall.hit_probability(0.3), 0.5);
    }

    // MBA-1317: WEZ rectangular hit probability.
    #[test]
    fn rect_hit_probability_checks_independent_axis_halves() {
        let results = make_results(vec![
            // Inside a 0.4 (lateral) x 0.6 (vertical) box: half-width 0.2, half-height 0.3.
            Vector3::new(0.0, 0.1, 0.1),
            // On the lateral edge (exactly half-width) -> counts as a hit ("<=").
            Vector3::new(0.0, 0.0, 0.2),
            // Outside laterally (just past half-width).
            Vector3::new(0.0, 0.0, 0.201),
            // Outside vertically (just past half-height).
            Vector3::new(0.0, 0.301, 0.0),
            // Shortfall marker: stays in the denominator, never a hit.
            Vector3::new(0.0, TARGET_NOT_REACHED_SENTINEL_M, 0.0),
        ]);
        // 2 hits out of 5 samples.
        assert!((results.rect_hit_probability(0.4, 0.6) - 0.4).abs() < 1e-12);
    }

    #[test]
    fn rect_hit_probability_matches_circular_hit_probability_for_a_centered_hit() {
        let results = make_results(vec![Vector3::new(0.0, 0.0, 0.0)]);
        assert_eq!(results.rect_hit_probability(0.5, 0.5), 1.0);
        assert_eq!(results.hit_probability(0.3), 1.0);
    }

    #[test]
    fn rect_hit_probability_is_zero_for_empty_or_nonpositive_dimensions() {
        let empty = make_results(vec![]);
        assert_eq!(empty.rect_hit_probability(1.0, 1.0), 0.0);

        let results = make_results(vec![Vector3::new(0.0, 0.0, 0.0)]);
        assert_eq!(results.rect_hit_probability(0.0, 1.0), 0.0);
        assert_eq!(results.rect_hit_probability(1.0, 0.0), 0.0);
        assert_eq!(results.rect_hit_probability(-1.0, 1.0), 0.0);
    }
}

#[cfg(test)]
mod monte_carlo_seeded_tests {
    use super::*;

    #[test]
    fn seeded_runs_are_deterministic_and_match_the_using_rng_path() {
        let inputs = BallisticInputs {
            muzzle_velocity: 800.0,
            ..BallisticInputs::default()
        };
        let params = MonteCarloParams {
            num_simulations: 64,
            target_distance: Some(200.0),
            ..MonteCarloParams::default()
        };

        let a = run_monte_carlo_with_wind_and_direction_std_dev_seeded(
            inputs.clone(),
            WindConditions::default(),
            params.clone(),
            0.01,
            42,
        )
        .expect("seeded run a");
        let b = run_monte_carlo_with_wind_and_direction_std_dev_seeded(
            inputs,
            WindConditions::default(),
            params,
            0.01,
            42,
        )
        .expect("seeded run b");

        assert_eq!(a.ranges.len(), b.ranges.len());
        for (ra, rb) in a.ranges.iter().zip(b.ranges.iter()) {
            assert_eq!(ra.to_bits(), rb.to_bits());
        }
        for (pa, pb) in a.impact_positions.iter().zip(b.impact_positions.iter()) {
            assert_eq!(pa.x.to_bits(), pb.x.to_bits());
            assert_eq!(pa.y.to_bits(), pb.y.to_bits());
            assert_eq!(pa.z.to_bits(), pb.z.to_bits());
        }
    }

    #[test]
    fn different_seeds_generally_produce_different_draws() {
        let inputs = BallisticInputs {
            muzzle_velocity: 800.0,
            ..BallisticInputs::default()
        };
        let params = MonteCarloParams {
            num_simulations: 32,
            velocity_std_dev: 5.0,
            target_distance: Some(200.0),
            ..MonteCarloParams::default()
        };

        let a = run_monte_carlo_with_wind_and_direction_std_dev_seeded(
            inputs.clone(),
            WindConditions::default(),
            params.clone(),
            0.0,
            1,
        )
        .expect("seeded run a");
        let b = run_monte_carlo_with_wind_and_direction_std_dev_seeded(
            inputs,
            WindConditions::default(),
            params,
            0.0,
            2,
        )
        .expect("seeded run b");

        assert_ne!(a.impact_velocities, b.impact_velocities);
    }
}

#[cfg(test)]
mod monte_carlo_powder_curve_tests {
    use super::*;
    use rand::{rngs::StdRng, SeedableRng};

    #[test]
    fn powder_curve_preserves_sampled_muzzle_velocity_dispersion() {
        let inputs = BallisticInputs {
            muzzle_velocity: 700.0,
            powder_temp_curve: Some(vec![(15.0, 800.0)]),
            powder_curve_temp_c: Some(15.0),
            ..BallisticInputs::default()
        };
        let params = MonteCarloParams {
            num_simulations: 16,
            velocity_std_dev: 20.0,
            angle_std_dev: 1e-12,
            bc_std_dev: 1e-12,
            wind_speed_std_dev: 1e-12,
            target_distance: Some(100.0),
            azimuth_std_dev: 1e-12,
            ..MonteCarloParams::default()
        };

        let mut rng = StdRng::seed_from_u64(0x5EED_1176);
        let results = run_monte_carlo_with_wind_and_direction_std_dev_using_rng(
            inputs,
            WindConditions::default(),
            params,
            0.0,
            &mut rng,
        )
        .expect("Monte Carlo solve");
        let min_velocity = results
            .impact_velocities
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        let max_velocity = results
            .impact_velocities
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);

        assert!(
            max_velocity - min_velocity > 1.0,
            "20 m/s muzzle spread collapsed after curve resolution: impact-velocity span={} m/s",
            max_velocity - min_velocity
        );
    }
}

#[cfg(test)]
mod monte_carlo_wind_sampling_tests {
    use super::*;
    use rand::{rngs::StdRng, SeedableRng};

    #[test]
    fn wind_speed_sigma_does_not_change_seeded_direction_draws() {
        let base_wind = WindConditions {
            speed: 100.0,
            direction: 0.37,
            vertical_speed: 0.0,
        };
        let narrow_speed = MonteCarloWindSampler::new(&base_wind, 0.5, 0.2).unwrap();
        let wide_speed = MonteCarloWindSampler::new(&base_wind, 4.0, 0.2).unwrap();
        let mut narrow_rng = StdRng::seed_from_u64(0x5EED_1223);
        let mut wide_rng = StdRng::seed_from_u64(0x5EED_1223);
        let mut speed_changed = false;

        for _ in 0..32 {
            let narrow = narrow_speed.sample(&mut narrow_rng);
            let wide = wide_speed.sample(&mut wide_rng);
            assert!(narrow.speed > 0.0 && wide.speed > 0.0);
            assert_eq!(narrow.direction.to_bits(), wide.direction.to_bits());
            speed_changed |= narrow.speed.to_bits() != wide.speed.to_bits();
        }
        assert!(
            speed_changed,
            "different speed sigmas must still vary speed draws"
        );
    }

    #[test]
    fn zero_direction_sigma_has_no_angular_jitter() {
        let base_wind = WindConditions {
            speed: 100.0,
            direction: 0.37,
            vertical_speed: 0.0,
        };
        let sampler = MonteCarloWindSampler::new(&base_wind, 4.0, 0.0).unwrap();
        let mut rng = StdRng::seed_from_u64(0x5EED_1223);
        let mut speed_changed = false;

        for _ in 0..32 {
            let wind = sampler.sample(&mut rng);
            speed_changed |= wind.speed.to_bits() != base_wind.speed.to_bits();
            assert_eq!(wind.direction.to_bits(), base_wind.direction.to_bits());
        }
        assert!(speed_changed, "speed uncertainty should remain active");
    }

    #[test]
    fn direction_sigma_controls_seeded_angular_spread_in_radians() {
        let base_wind = WindConditions {
            speed: 100.0,
            direction: 0.37,
            vertical_speed: 0.0,
        };
        let narrow = MonteCarloWindSampler::new(&base_wind, 4.0, 0.1).unwrap();
        let wide = MonteCarloWindSampler::new(&base_wind, 4.0, 0.2).unwrap();
        let mut narrow_rng = StdRng::seed_from_u64(0x5EED_1223);
        let mut wide_rng = StdRng::seed_from_u64(0x5EED_1223);
        let mut nonzero_direction_draw = false;

        for _ in 0..32 {
            let narrow_wind = narrow.sample(&mut narrow_rng);
            let wide_wind = wide.sample(&mut wide_rng);
            assert_eq!(narrow_wind.speed.to_bits(), wide_wind.speed.to_bits());

            let narrow_delta = narrow_wind.direction - base_wind.direction;
            let wide_delta = wide_wind.direction - base_wind.direction;
            assert!((wide_delta - 2.0 * narrow_delta).abs() < 1e-12);
            nonzero_direction_draw |= narrow_delta.abs() > 1e-6;
        }
        assert!(
            nonzero_direction_draw,
            "positive radians sigma must vary direction"
        );
    }

    #[test]
    fn direction_sigma_rejects_negative_or_nonfinite_values() {
        let base_wind = WindConditions::default();
        for sigma in [-0.1, f64::NAN, f64::INFINITY] {
            assert!(MonteCarloWindSampler::new(&base_wind, 1.0, sigma).is_err());
        }
    }

    #[test]
    fn base_vertical_wind_rides_into_every_mc_sample() {
        // MBA-728: vertical wind is a systematic input, not a dispersion source —
        // every sampled wind must carry the base vertical un-dispersed. (Before
        // this fix, samples dropped it, biasing the whole MC cloud vs the baseline.)
        use rand::SeedableRng;
        let base_wind = WindConditions { vertical_speed: 4.2, ..Default::default() };
        let sampler = MonteCarloWindSampler::new(&base_wind, 1.0, 0.2).unwrap();
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);
        for _ in 0..32 {
            let w = sampler.sample(&mut rng);
            assert_eq!(w.vertical_speed, 4.2);
        }
    }

    #[test]
    fn negative_speed_sample_reverses_wind_direction() {
        let direction = 0.25;
        let signed_speed = -2.5;
        let wind = wind_from_signed_speed_sample(signed_speed, direction, 0.0);
        let positive_wind = wind_from_signed_speed_sample(2.5, direction, 0.0);

        assert_eq!(wind.speed, 2.5);
        assert!(
            (wind.direction - (direction + std::f64::consts::PI)).abs() < f64::EPSILON,
            "negative speed must reverse direction by pi: got {}",
            wind.direction
        );
        assert_eq!(positive_wind.speed, 2.5);
        assert_eq!(positive_wind.direction, direction);

        let normalized_x = -wind.speed * wind.direction.cos();
        let normalized_z = -wind.speed * wind.direction.sin();
        let signed_x = -signed_speed * direction.cos();
        let signed_z = -signed_speed * direction.sin();
        assert!((normalized_x - signed_x).abs() < 1e-12);
        assert!((normalized_z - signed_z).abs() < 1e-12);
    }
}

#[cfg(test)]
mod bc_fit_objective_tests {
    use super::*;

    fn velocity_point(range_m: f64, velocity_mps: f64) -> TrajectoryPoint {
        TrajectoryPoint {
            time: 0.0,
            position: Vector3::new(range_m, 0.0, 0.0),
            velocity_magnitude: velocity_mps,
            kinetic_energy: 0.0,
            drag_coefficient: None,
        }
    }

    #[test]
    fn candidate_that_misses_an_observation_has_no_score() {
        let trajectory = vec![velocity_point(0.0, 800.0), velocity_point(100.0, 700.0)];
        let observations = vec![(50.0, 750.0), (150.0, 600.0)];

        assert!(
            fit_residual_sse(&trajectory, &observations, BcFitMode::Velocity, 0.0).is_none(),
            "a candidate that reaches only one of two observations must not compete on partial SSE"
        );

        let complete_observations = vec![(50.0, 740.0), (100.0, 680.0)];
        assert_eq!(
            fit_residual_sse(
                &trajectory,
                &complete_observations,
                BcFitMode::Velocity,
                0.0,
            ),
            Some(500.0)
        );
    }
}

#[cfg(test)]
mod cluster_bc_reference_space_tests {
    use super::*;

    fn acceleration_at_1100_fps(inputs: BallisticInputs) -> Vector3<f64> {
        let solver = TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        let position = Vector3::zeros();
        let velocity = Vector3::new(1100.0 / 3.28084, 0.0, 0.0);
        let (density, _, temp_c, pressure_hpa) = solver.resolved_atmosphere();
        solver.calculate_acceleration(
            &position,
            &velocity,
            &Vector3::zeros(),
            (temp_c, pressure_hpa, density / 1.225),
        )
    }

    #[test]
    fn solver_passes_g7_reference_model_to_cluster_classifier() {
        let inputs = BallisticInputs {
            bc_value: 0.190,
            bc_type: DragModel::G7,
            bullet_mass: 77.0 * crate::constants::GRAINS_TO_KG,
            bullet_diameter: 0.224 * 0.0254,
            use_cluster_bc: true,
            ..BallisticInputs::default()
        };

        let solver = TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        let corrected = solver.apply_cluster_bc_correction(0.190, 2800.0);

        assert!(
            (corrected / 0.190 - 1.004).abs() < 1e-12,
            "solver selected the wrong G7 cluster multiplier: {}",
            corrected / 0.190
        );
    }

    #[test]
    fn velocity_bc_segments_are_not_cluster_corrected_twice() {
        let segmented_clustered = BallisticInputs {
            bc_value: 0.5,
            bc_type: DragModel::G7,
            use_bc_segments: true,
            bc_segments_data: Some(vec![
                crate::BCSegmentData {
                    velocity_min: 0.0,
                    velocity_max: 1_600.0,
                    bc_value: 0.4,
                },
                crate::BCSegmentData {
                    velocity_min: 1_600.0,
                    velocity_max: 5_000.0,
                    bc_value: 0.45,
                },
            ]),
            use_cluster_bc: true,
            ..BallisticInputs::default()
        };
        let mut segmented_only = segmented_clustered.clone();
        segmented_only.use_cluster_bc = false;
        let mut constant_clustered = segmented_clustered.clone();
        constant_clustered.bc_value = 0.4;
        constant_clustered.bc_segments_data = None;

        let stacked = acceleration_at_1100_fps(segmented_clustered);
        let segment_only = acceleration_at_1100_fps(segmented_only);
        let cluster_only = acceleration_at_1100_fps(constant_clustered);

        assert!(
            (stacked.x - segment_only.x).abs() < 1e-12,
            "segment BC already owns the velocity shape: stacked ax={} segment-only ax={}",
            stacked.x,
            segment_only.x
        );
        assert!(
            (cluster_only.x - segment_only.x).abs() > 1e-6,
            "cluster correction must remain active for a constant BC"
        );
    }

    #[test]
    fn mach_bc_segments_are_not_cluster_corrected_twice() {
        let mach_segmented_clustered = BallisticInputs {
            bc_value: 0.5,
            bc_type: DragModel::G7,
            use_bc_segments: false,
            bc_segments: Some(vec![(0.5, 0.3), (1.5, 0.5)]),
            use_cluster_bc: true,
            ..BallisticInputs::default()
        };
        let mut mach_segmented_only = mach_segmented_clustered.clone();
        mach_segmented_only.use_cluster_bc = false;

        let stacked = acceleration_at_1100_fps(mach_segmented_clustered);
        let segment_only = acceleration_at_1100_fps(mach_segmented_only);

        assert!(
            (stacked.x - segment_only.x).abs() < 1e-12,
            "Mach segment BC already owns the velocity shape: stacked ax={} segment-only ax={}",
            stacked.x,
            segment_only.x
        );
    }
}

#[cfg(test)]
mod velocity_bc_flag_tests {
    use super::*;

    fn acceleration_at_600_mps(inputs: BallisticInputs) -> Vector3<f64> {
        let solver = TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        let (density, _, temp_c, pressure_hpa) = solver.resolved_atmosphere();
        solver.calculate_acceleration(
            &Vector3::zeros(),
            &Vector3::new(600.0, 0.0, 0.0),
            &Vector3::zeros(),
            (temp_c, pressure_hpa, density / 1.225),
        )
    }

    #[test]
    fn velocity_bc_data_requires_opt_in_in_trajectory_solver() {
        let scalar_inputs = BallisticInputs {
            bc_value: 0.5,
            bc_type: DragModel::G7,
            ..BallisticInputs::default()
        };
        let mut disabled_inputs = scalar_inputs.clone();
        disabled_inputs.bc_segments_data = Some(vec![crate::BCSegmentData {
            velocity_min: 0.0,
            velocity_max: 4_000.0,
            bc_value: 0.46,
        }]);
        disabled_inputs.use_bc_segments = false;
        let mut enabled_inputs = disabled_inputs.clone();
        enabled_inputs.use_bc_segments = true;
        let mut mach_only_inputs = scalar_inputs.clone();
        mach_only_inputs.bc_segments = Some(vec![(0.0, 0.4), (3.0, 0.4)]);
        let mut disabled_with_both = mach_only_inputs.clone();
        disabled_with_both.bc_segments_data = disabled_inputs.bc_segments_data.clone();

        let scalar = acceleration_at_600_mps(scalar_inputs);
        let disabled = acceleration_at_600_mps(disabled_inputs);
        let enabled = acceleration_at_600_mps(enabled_inputs);
        let mach_only = acceleration_at_600_mps(mach_only_inputs);
        let disabled_with_both = acceleration_at_600_mps(disabled_with_both);

        assert_eq!(
            disabled.x.to_bits(),
            scalar.x.to_bits(),
            "a populated velocity table must not change drag while use_bc_segments is false"
        );
        assert!(
            enabled.x < disabled.x - 1.0,
            "enabling the lower BC table must increase drag: disabled ax={} enabled ax={}",
            disabled.x,
            enabled.x
        );
        assert_eq!(
            disabled_with_both.x.to_bits(),
            mach_only.x.to_bits(),
            "disabling velocity data must fall through to an explicit Mach table"
        );
    }
}

#[cfg(test)]
mod mach_bc_segment_tests {
    use super::*;

    #[test]
    fn trajectory_solver_interpolates_explicit_mach_bc_segments() {
        let segmented_inputs = BallisticInputs {
            bc_value: 0.8,
            use_bc_segments: false,
            bc_segments: Some(vec![(1.0, 0.2), (2.0, 0.4)]),
            bc_segments_data: None,
            ..BallisticInputs::default()
        };

        let mut expected_inputs = segmented_inputs.clone();
        expected_inputs.bc_value = 0.3;
        expected_inputs.bc_segments = None;

        let atmosphere = AtmosphericConditions::default();
        let segmented_solver = TrajectorySolver::new(
            segmented_inputs,
            WindConditions::default(),
            atmosphere.clone(),
        );
        let expected_solver = TrajectorySolver::new(
            expected_inputs,
            WindConditions::default(),
            atmosphere,
        );
        let position = Vector3::zeros();
        let (density, _, temp_c, pressure_hpa) = segmented_solver.resolved_atmosphere();
        let (_, local_speed_of_sound) = crate::atmosphere::get_local_atmosphere_humid(
            segmented_solver.atmosphere.altitude,
            segmented_solver.atmosphere.altitude,
            temp_c,
            pressure_hpa,
            density / 1.225,
            segmented_solver.atmosphere.humidity,
        );
        let velocity = Vector3::new(1.5 * local_speed_of_sound, 0.0, 0.0);
        let resolved_atmo = (temp_c, pressure_hpa, density / 1.225);

        let segmented_acceleration = segmented_solver.calculate_acceleration(
            &position,
            &velocity,
            &Vector3::zeros(),
            resolved_atmo,
        );
        let expected_acceleration = expected_solver.calculate_acceleration(
            &position,
            &velocity,
            &Vector3::zeros(),
            resolved_atmo,
        );

        assert!(
            (segmented_acceleration.x - expected_acceleration.x).abs() < 1e-12,
            "Mach 1.5 must interpolate BC 0.3: segmented ax={} expected ax={}",
            segmented_acceleration.x,
            expected_acceleration.x
        );
    }
}

#[cfg(test)]
mod custom_drag_table_validation_tests {
    use super::*;

    #[test]
    fn solve_accepts_zero_bc_when_custom_table_present() {
        let inputs = BallisticInputs {
            bc_value: 0.0, // ignored when a table is set
            bullet_mass: 0.0106,
            bullet_diameter: 0.00782,
            muzzle_velocity: 850.0,
            custom_drag_table: Some(crate::drag::DragTable::new(
                vec![0.5, 1.0, 2.0, 3.0],
                vec![0.23, 0.40, 0.30, 0.26],
            )),
            ..BallisticInputs::default()
        };
        let solver = TrajectorySolver::new(inputs, WindConditions::default(), AtmosphericConditions::default());
        // Must not error on the bc_value gate.
        assert!(solver.solve().is_ok());
    }

    #[test]
    fn solve_still_requires_bc_without_table() {
        let inputs = BallisticInputs {
            bc_value: 0.0,
            bullet_mass: 0.0106,
            bullet_diameter: 0.00782,
            muzzle_velocity: 850.0,
            ..BallisticInputs::default()
        };
        let solver = TrajectorySolver::new(inputs, WindConditions::default(), AtmosphericConditions::default());
        assert!(solver.solve().is_err());
    }
}

/// MBA-1356: cd_scale — a whole-curve multiplier on the custom-deck Cd.
#[cfg(test)]
mod cd_scale_tests {
    use super::*;

    fn deck() -> crate::drag::DragTable {
        crate::drag::DragTable::new(vec![0.5, 1.0, 2.0, 3.0], vec![0.23, 0.40, 0.30, 0.26])
    }

    fn deck_inputs(cd_scale: f64) -> BallisticInputs {
        BallisticInputs {
            bullet_mass: 0.0106,
            bullet_diameter: 0.00782,
            muzzle_velocity: 850.0,
            custom_drag_table: Some(deck()),
            cd_scale,
            ..BallisticInputs::default()
        }
    }

    #[test]
    fn default_cd_scale_is_one() {
        assert_eq!(BallisticInputs::default().cd_scale, 1.0);
    }

    /// (a) Default invariance: omitting `cd_scale` (picked up from `..Default::default()`)
    /// must be byte-identical to explicitly setting `1.0`, on the custom-deck Cd itself.
    #[test]
    fn cd_scale_absent_is_byte_identical_to_explicit_one() {
        let omitted = BallisticInputs {
            bullet_mass: 0.0106,
            bullet_diameter: 0.00782,
            muzzle_velocity: 850.0,
            custom_drag_table: Some(deck()),
            ..BallisticInputs::default()
        };
        let explicit = BallisticInputs {
            cd_scale: 1.0,
            ..omitted.clone()
        };

        let solver_omitted =
            TrajectorySolver::new(omitted, WindConditions::default(), AtmosphericConditions::default());
        let solver_explicit =
            TrajectorySolver::new(explicit, WindConditions::default(), AtmosphericConditions::default());

        let cd_omitted = solver_omitted.calculate_drag_coefficient(700.0, 340.0);
        let cd_explicit = solver_explicit.calculate_drag_coefficient(700.0, 340.0);
        assert_eq!(
            cd_omitted.to_bits(),
            cd_explicit.to_bits(),
            "default cd_scale must be bit-identical to an explicit 1.0"
        );

        // And a full custom-deck solve (the existing pre-MBA-1356 test surface) must still
        // succeed unchanged with the field simply absent from the literal.
        let result = solver_omitted.solve();
        assert!(result.is_ok(), "existing custom-deck solves must pass unchanged");
    }

    /// The custom-deck interpolation site multiplies the deck's Cd by `cd_scale` exactly.
    #[test]
    fn cd_scale_multiplies_the_interpolated_cd_exactly() {
        let velocity = 700.0;
        let speed_of_sound = 340.0;
        let mach = velocity / speed_of_sound;
        let expected_unscaled = deck().interpolate(mach);

        for &scale in &[0.90, 1.0, 1.10, 1.5] {
            let solver = TrajectorySolver::new(
                deck_inputs(scale),
                WindConditions::default(),
                AtmosphericConditions::default(),
            );
            let cd = solver.calculate_drag_coefficient(velocity, speed_of_sound);
            assert!(
                (cd - expected_unscaled * scale).abs() < 1e-12,
                "scale={scale}: cd={cd} expected={}",
                expected_unscaled * scale
            );
        }
    }

    /// (b) Scale-direction test, cli_api (RK4/RK45) solver path: a higher cd_scale means more
    /// drag, so the solve loses more velocity (and correspondingly drops more) over the same
    /// downrange distance; a lower cd_scale loses less.
    #[test]
    fn cd_scale_direction_on_cli_api_solver() {
        let solve = |scale: f64| {
            TrajectorySolver::new(
                deck_inputs(scale),
                WindConditions::default(),
                AtmosphericConditions::default(),
            )
            .solve()
            .expect("custom-deck solve should succeed")
        };

        let baseline = solve(1.0);
        let scaled_up = solve(1.10);
        let scaled_down = solve(0.90);

        assert!(
            scaled_up.impact_velocity < baseline.impact_velocity,
            "cd_scale=1.10 must increase drag -> lower impact velocity: base={} up={}",
            baseline.impact_velocity,
            scaled_up.impact_velocity
        );
        assert!(
            scaled_down.impact_velocity > baseline.impact_velocity,
            "cd_scale=0.90 must decrease drag -> higher impact velocity: base={} down={}",
            baseline.impact_velocity,
            scaled_down.impact_velocity
        );
    }

    /// (d) `validate_for_solve` rejects a non-finite or non-positive `cd_scale`.
    #[test]
    fn validate_for_solve_rejects_invalid_cd_scale() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let solver = TrajectorySolver::new(
                deck_inputs(bad),
                WindConditions::default(),
                AtmosphericConditions::default(),
            );
            assert!(
                solver.solve().is_err(),
                "cd_scale={bad} must be rejected by validate_for_solve"
            );
        }
    }

    /// `require_positive("cd_scale", ...)` is unconditional in `validate_for_solve` (see the
    /// comment there) — it must reject an invalid scale even with NO custom drag table present,
    /// not just on the deck path the other tests in this module exercise. Otherwise a caller
    /// that sets an invalid `cd_scale` without ever touching `custom_drag_table` would slip
    /// through unvalidated (the field would simply go unread on this path, but the validation
    /// gate itself must not silently skip it).
    #[test]
    fn validate_for_solve_rejects_invalid_cd_scale_without_a_custom_drag_table() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let inputs = BallisticInputs {
                bc_value: 0.5,
                bc_type: crate::DragModel::G1,
                bullet_mass: 0.0106,
                bullet_diameter: 0.00782,
                muzzle_velocity: 850.0,
                cd_scale: bad,
                ..BallisticInputs::default()
            };
            assert!(inputs.custom_drag_table.is_none(), "precondition: no custom deck");
            let solver = TrajectorySolver::new(
                inputs,
                WindConditions::default(),
                AtmosphericConditions::default(),
            );
            assert!(
                solver.solve().is_err(),
                "cd_scale={bad} must be rejected by validate_for_solve even without a custom \
                 drag table"
            );
        }
    }

    /// cd_scale must be inert on the standard G-model/BC path (no custom_drag_table): the
    /// scale is read only inside the `custom_drag_table` branch, so a scale far from 1.0 must
    /// not perturb a plain G1/G7 solve at all.
    #[test]
    fn cd_scale_is_inert_without_a_custom_drag_table() {
        let make = |cd_scale: f64| BallisticInputs {
            bc_value: 0.5,
            bc_type: crate::DragModel::G1,
            bullet_mass: 0.0106,
            bullet_diameter: 0.00782,
            muzzle_velocity: 850.0,
            cd_scale,
            ..BallisticInputs::default()
        };
        let solver_neutral = TrajectorySolver::new(
            make(1.0),
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        let solver_far = TrajectorySolver::new(
            make(1.5),
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        let cd_neutral = solver_neutral.calculate_drag_coefficient(700.0, 340.0);
        let cd_far = solver_far.calculate_drag_coefficient(700.0, 340.0);
        assert_eq!(
            cd_neutral.to_bits(),
            cd_far.to_bits(),
            "cd_scale must not affect the G-model/BC drag path"
        );
    }

    /// (c) Cross-solver consistency: cli_api (RK4/RK45), derivatives-driven, and fast_trajectory
    /// all shift together (same relative direction) under cd_scale != 1.0 on the same deck.
    #[test]
    fn cd_scale_shifts_all_three_solver_paths_in_the_same_direction() {
        // --- cli_api ---
        let cli_solve = |scale: f64| {
            TrajectorySolver::new(
                deck_inputs(scale),
                WindConditions::default(),
                AtmosphericConditions::default(),
            )
            .solve()
            .expect("cli_api custom-deck solve should succeed")
        };
        let cli_baseline = cli_solve(1.0);
        let cli_scaled = cli_solve(1.10);
        assert!(
            cli_scaled.impact_velocity < cli_baseline.impact_velocity,
            "cli_api: cd_scale=1.10 must lower impact velocity"
        );

        // --- derivatives (the RK4/RK45 generic integrator's kernel) ---
        let derivatives_accel_x = |scale: f64| {
            let inputs = deck_inputs(scale);
            crate::derivatives::compute_derivatives(
                nalgebra::Vector3::zeros(),
                nalgebra::Vector3::new(700.0, 0.0, 0.0),
                &inputs,
                nalgebra::Vector3::zeros(),
                (1.225, 340.0, 0.0, 0.0),
                inputs.bc_value,
                None,
                0.0,
                None,
            )[3]
        };
        let deriv_baseline = derivatives_accel_x(1.0);
        let deriv_scaled = derivatives_accel_x(1.10);
        assert!(
            deriv_scaled < deriv_baseline,
            "derivatives: cd_scale=1.10 must make x-acceleration more negative (more drag): \
             base={deriv_baseline} scaled={deriv_scaled}"
        );

        // --- fast_trajectory (fast_integrate) ---
        let fast_final_speed = |scale: f64| {
            let inputs = deck_inputs(scale);
            let wind_sock = crate::wind::WindSock::new(vec![]);
            let params = crate::fast_trajectory::FastIntegrationParams {
                horiz: 500.0,
                vert: 0.0,
                initial_state: [0.0, 0.0, 0.0, 850.0, 0.0, 0.0],
                t_span: (0.0, 5.0),
                atmo_params: (0.0, 15.0, 1013.25, 1.0),
                atmo_sock: None,
            };
            let solution = crate::fast_trajectory::fast_integrate(&inputs, &wind_sock, params);
            assert!(solution.success, "fast_integrate must succeed for scale={scale}");
            let last = solution.t.len() - 1;
            let (vx, vy, vz) = (
                solution.y[3][last],
                solution.y[4][last],
                solution.y[5][last],
            );
            (vx * vx + vy * vy + vz * vz).sqrt()
        };
        let fast_baseline = fast_final_speed(1.0);
        let fast_scaled = fast_final_speed(1.10);
        assert!(
            fast_scaled < fast_baseline,
            "fast_trajectory: cd_scale=1.10 must lower final speed: base={fast_baseline} scaled={fast_scaled}"
        );
    }
}

#[cfg(test)]
mod humid_local_mach_tests {
    use super::*;

    fn solver_with_station_humidity(humidity_percent: f64) -> TrajectorySolver {
        let inputs = BallisticInputs {
            custom_drag_table: Some(crate::drag::DragTable::new(vec![0.5, 1.5], vec![0.1, 1.1])),
            ..BallisticInputs::default()
        };
        TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions {
                temperature: 30.0,
                pressure: 1013.25,
                humidity: humidity_percent,
                altitude: 0.0,
            },
        )
    }

    fn acceleration(solver: &TrajectorySolver, base_ratio: f64) -> Vector3<f64> {
        solver.calculate_acceleration(
            &Vector3::zeros(),
            &Vector3::new(350.0, 0.0, 0.0),
            &Vector3::zeros(),
            (30.0, 1013.25, base_ratio),
        )
    }

    #[test]
    fn local_mach_uses_station_humidity_when_density_is_held_constant() {
        let dry = acceleration(&solver_with_station_humidity(0.0), 1.0);
        let humid = acceleration(&solver_with_station_humidity(100.0), 1.0);

        assert!(
            humid.x > dry.x,
            "humid sound speed should lower Mach and drag on the rising test curve: dry ax={} humid ax={}",
            dry.x,
            humid.x
        );
    }

    #[test]
    fn active_atmosphere_zone_uses_zone_humidity_instead_of_station_humidity() {
        let zone_humidity = 80.0;
        let zone_ratio =
            crate::atmosphere::calculate_air_density_cimp(30.0, 1013.25, zone_humidity) / 1.225;
        let station_solver = solver_with_station_humidity(zone_humidity);
        let mut zoned_solver = solver_with_station_humidity(0.0);
        zoned_solver.set_atmo_segments(vec![(30.0, 1013.25, zone_humidity, 1_000.0)]);

        let station = acceleration(&station_solver, zone_ratio);
        let zoned = acceleration(&zoned_solver, zone_ratio);

        assert!(
            (zoned - station).norm() < 1e-12,
            "active zone T/P/RH should override the station atmosphere: station={station:?} zoned={zoned:?}"
        );
    }
}

#[cfg(test)]
mod inclined_atmosphere_frame_tests {
    use super::*;

    fn expected_shot_frame_vector(level: Vector3<f64>, angle: f64) -> Vector3<f64> {
        let (sin_angle, cos_angle) = angle.sin_cos();
        Vector3::new(
            level.x * cos_angle + level.y * sin_angle,
            -level.x * sin_angle + level.y * cos_angle,
            level.z,
        )
    }

    #[test]
    fn inclined_positions_at_same_world_altitude_have_same_solver_acceleration() {
        let angle = std::f64::consts::FRAC_PI_6;
        let inputs = BallisticInputs {
            shooting_angle: angle,
            ..BallisticInputs::default()
        };
        let atmosphere = AtmosphericConditions {
            altitude: 100.0,
            ..AtmosphericConditions::default()
        };
        let solver = TrajectorySolver::new(inputs, WindConditions::default(), atmosphere);
        let (density, _, temp_c, pressure_hpa) = solver.resolved_atmosphere();
        let resolved_atmo = (temp_c, pressure_hpa, density / 1.225);
        let velocity = Vector3::new(600.0, 0.0, 0.0);
        let along_slant = Vector3::new(1_000.0, 0.0, 0.0);
        let across_slant = Vector3::new(0.0, 500.0 / angle.cos(), 0.0);

        let a = solver.calculate_acceleration(
            &along_slant,
            &velocity,
            &Vector3::zeros(),
            resolved_atmo,
        );
        let b = solver.calculate_acceleration(
            &across_slant,
            &velocity,
            &Vector3::zeros(),
            resolved_atmo,
        );

        assert!(
            (a - b).norm() < 1e-10,
            "solver acceleration differs at equal world altitude: {a:?} vs {b:?}"
        );
    }

    #[test]
    fn inclined_headwind_is_rotated_into_solver_frame() {
        let angle = std::f64::consts::FRAC_PI_6;
        let inputs = BallisticInputs {
            shooting_angle: angle,
            ..BallisticInputs::default()
        };
        let solver = TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        let level_headwind = Vector3::new(-100.0, 0.0, 0.0);
        let velocity = expected_shot_frame_vector(level_headwind, angle);
        let (density, _, temp_c, pressure_hpa) = solver.resolved_atmosphere();
        let actual = solver.calculate_acceleration(
            &Vector3::zeros(),
            &velocity,
            &level_headwind,
            (temp_c, pressure_hpa, density / 1.225),
        );

        assert!(
            (actual - solver.gravity_acceleration()).norm() < 1e-12,
            "co-moving horizontal wind must leave only shot-frame gravity: {actual:?}"
        );
    }

    #[test]
    fn inclined_coriolis_is_rotated_into_solver_frame() {
        let angle = std::f64::consts::FRAC_PI_6;
        let latitude_deg = 45.0_f64;
        let shot_azimuth = 0.4_f64;
        let velocity = Vector3::new(600.0, 20.0, 5.0);
        let base_inputs = BallisticInputs {
            shooting_angle: angle,
            latitude: Some(latitude_deg),
            shot_azimuth,
            ..BallisticInputs::default()
        };
        let acceleration = |enable_coriolis| {
            let mut inputs = base_inputs.clone();
            inputs.enable_coriolis = enable_coriolis;
            let solver = TrajectorySolver::new(
                inputs,
                WindConditions::default(),
                AtmosphericConditions::default(),
            );
            let (density, _, temp_c, pressure_hpa) = solver.resolved_atmosphere();
            solver.calculate_acceleration(
                &Vector3::zeros(),
                &velocity,
                &Vector3::zeros(),
                (temp_c, pressure_hpa, density / 1.225),
            )
        };

        let omega_earth = 7.2921159e-5_f64;
        let latitude = latitude_deg.to_radians();
        let level_omega = Vector3::new(
            omega_earth * latitude.cos() * shot_azimuth.cos(),
            omega_earth * latitude.sin(),
            -omega_earth * latitude.cos() * shot_azimuth.sin(),
        );
        let expected = -2.0 * expected_shot_frame_vector(level_omega, angle).cross(&velocity);
        let actual = acceleration(true) - acceleration(false);

        assert!(
            (actual - expected).norm() < 1e-12,
            "inclined Coriolis mismatch: actual={actual:?}, expected={expected:?}"
        );
    }
}

#[cfg(test)]
mod terminal_range_interpolation_tests {
    use super::*;

    #[test]
    fn terminal_finalizer_selects_the_earliest_crossed_boundary() {
        let inputs = BallisticInputs {
            ground_threshold: 0.0,
            ..BallisticInputs::default()
        };
        let mut solver = TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        solver.set_max_range(120.0);

        let previous_speed = 700.0;
        let mut points = vec![TrajectoryPoint {
            time: 99.0,
            position: Vector3::new(90.0, 1.0, -1.0),
            velocity_magnitude: previous_speed,
            kinetic_energy: 0.5 * solver.inputs.bullet_mass * previous_speed.powi(2),
            drag_coefficient: None,
        }];
        let mut max_height = 1.0;
        let termination = solver
            .append_terminal_endpoint(
                &mut points,
                Vector3::new(130.0, -3.0, 3.0),
                Vector3::new(600.0, 0.0, 0.0),
                101.0,
                &mut max_height,
            )
            .expect("the final step brackets supported boundaries");

        assert_eq!(termination, TrajectoryTermination::GroundThreshold);
        assert_eq!(points.len(), 2);
        let terminal = points.last().expect("terminal point");
        assert_eq!(terminal.time, 99.5);
        assert_eq!(terminal.position, Vector3::new(100.0, 0.0, 0.0));
        assert_eq!(terminal.velocity_magnitude, 675.0);
        assert_eq!(
            terminal.kinetic_energy,
            0.5 * solver.inputs.bullet_mass * 675.0_f64.powi(2)
        );

        // At x=100 the range and ground crossings tie; physical ground impact wins explicitly.
        solver.set_max_range(100.0);
        let mut tied_points = vec![points[0].clone()];
        assert_eq!(
            solver
                .append_terminal_endpoint(
                    &mut tied_points,
                    Vector3::new(130.0, -3.0, 3.0),
                    Vector3::new(600.0, 0.0, 0.0),
                    101.0,
                    &mut max_height,
                )
                .expect("tied boundaries remain a valid terminal"),
            TrajectoryTermination::GroundThreshold
        );
    }

    #[test]
    fn sub_ulp_terminal_crossing_replaces_instead_of_duplicating_range() {
        let ground_threshold = f64::from_bits(1.0_f64.to_bits() - 1);
        let inputs = BallisticInputs {
            ground_threshold,
            ..BallisticInputs::default()
        };
        let mut solver = TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        solver.set_max_range(1_000.0);

        let speed = 700.0;
        let mut points = vec![TrajectoryPoint {
            time: 0.0,
            position: Vector3::new(100.0, 1.0, 0.0),
            velocity_magnitude: speed,
            kinetic_energy: 0.5 * solver.inputs.bullet_mass * speed.powi(2),
            drag_coefficient: None,
        }];
        let mut max_height = 1.0;
        let termination = solver
            .append_terminal_endpoint(
                &mut points,
                Vector3::new(101.0, 0.0, 0.0),
                Vector3::new(699.0, 0.0, 0.0),
                1.0,
                &mut max_height,
            )
            .expect("sub-ULP ground crossing remains representable as one terminal state");

        assert_eq!(termination, TrajectoryTermination::GroundThreshold);
        assert_eq!(points.len(), 1);
        assert_eq!(points[0].position.x, 100.0);
        assert_eq!(points[0].position.y.to_bits(), ground_threshold.to_bits());
        assert!(points[0].time > 0.0);
    }

    #[test]
    fn every_solver_appends_an_exact_max_range_endpoint() {
        let target_range = 0.1;
        let modes = [
            ("Euler", false, false),
            ("RK4", true, false),
            ("RK45", true, true),
        ];

        for (name, use_rk4, use_adaptive_rk45) in modes {
            let inputs = BallisticInputs {
                use_rk4,
                use_adaptive_rk45,
                ground_threshold: f64::NEG_INFINITY,
                enable_trajectory_sampling: true,
                sample_interval: target_range,
                ..BallisticInputs::default()
            };
            let mut solver = TrajectorySolver::new(
                inputs,
                WindConditions::default(),
                AtmosphericConditions::default(),
            );
            solver.set_max_range(target_range);

            let result = solver.solve().expect("short-range solve should succeed");
            let terminal = result.points.last().expect("terminal point is missing");
            let muzzle = result.points.first().expect("muzzle point is missing");

            assert_eq!(result.termination, TrajectoryTermination::MaxRange);
            assert_eq!(
                terminal.position.x.to_bits(),
                target_range.to_bits(),
                "{name} did not terminate exactly at max_range"
            );
            assert_eq!(result.max_range.to_bits(), target_range.to_bits());
            assert!(
                result.time_of_flight > 0.0 && result.time_of_flight < solver.time_step,
                "{name} terminal time was not interpolated within the crossing step: {}",
                result.time_of_flight
            );
            assert_eq!(result.time_of_flight.to_bits(), terminal.time.to_bits());
            assert_eq!(
                result.impact_velocity.to_bits(),
                terminal.velocity_magnitude.to_bits()
            );
            assert_eq!(
                result.impact_energy.to_bits(),
                terminal.kinetic_energy.to_bits()
            );
            let expected_energy = 0.5 * solver.inputs.bullet_mass * result.impact_velocity.powi(2);
            assert!((result.impact_energy - expected_energy).abs() < 1e-12);
            assert!(terminal.velocity_magnitude < muzzle.velocity_magnitude);
            assert!(terminal.kinetic_energy < muzzle.kinetic_energy);

            let terminal_sample = result
                .sampled_points
                .as_ref()
                .and_then(|samples| samples.last())
                .expect("terminal trajectory sample is missing");
            assert_eq!(
                terminal_sample.distance_m.to_bits(),
                target_range.to_bits(),
                "{name} sampling did not include max_range"
            );
            assert_eq!(
                terminal_sample.time_s.to_bits(),
                result.time_of_flight.to_bits()
            );
            assert_eq!(
                terminal_sample.velocity_mps.to_bits(),
                result.impact_velocity.to_bits()
            );
            assert!((terminal_sample.energy_j - result.impact_energy).abs() < 1e-12);
        }
    }
}

#[cfg(test)]
mod precession_inertia_wiring_tests {
    use super::*;

    #[test]
    fn solver_uses_projectile_specific_moments_of_inertia() {
        let mass_kg = 55.0 * crate::constants::GRAINS_TO_KG;
        let caliber_m = 0.224 * 0.0254;
        let length_m = 0.75 * 0.0254;
        let inputs = BallisticInputs {
            bullet_mass: mass_kg,
            bullet_diameter: caliber_m,
            bullet_length: length_m,
            muzzle_velocity: 800.0,
            twist_rate: 7.0,
            enable_precession_nutation: true,
            use_rk4: false,
            use_adaptive_rk45: false,
            ..BallisticInputs::default()
        };
        let mut solver = TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        solver.set_max_range(0.1);

        let (air_density, speed_of_sound, _, _) = solver.resolved_atmosphere();
        let velocity_mps = solver.inputs.muzzle_velocity;
        let velocity_fps = velocity_mps * 3.28084;
        let twist_rate_ft = solver.inputs.twist_rate / 12.0;
        let spin_rate_rad_s = (velocity_fps / twist_rate_ft) * 2.0 * std::f64::consts::PI;
        let initial_state = AngularState {
            pitch_angle: 0.001,
            yaw_angle: 0.001,
            pitch_rate: 0.0,
            yaw_rate: 0.0,
            precession_angle: 0.0,
            nutation_phase: 0.0,
        };
        let params = PrecessionNutationParams {
            mass_kg,
            caliber_m,
            length_m,
            spin_rate_rad_s,
            spin_inertia: crate::spin_decay::calculate_moment_of_inertia(
                mass_kg, caliber_m, length_m, "ogive",
            ),
            transverse_inertia: crate::pitch_damping::calculate_transverse_moment_of_inertia(
                mass_kg, caliber_m, length_m, "ogive",
            ),
            velocity_mps,
            air_density_kg_m3: air_density,
            mach: velocity_mps / speed_of_sound,
            pitch_damping_coeff: PitchDampingCoefficients::default().subsonic,
            nutation_damping_factor: 0.05,
        };
        let expected = calculate_combined_angular_motion(
            &params,
            &initial_state,
            0.0,
            solver.time_step,
            0.001,
        );
        let actual = solver
            .solve()
            .expect("one-step solve should succeed")
            .angular_state
            .expect("precession/nutation was enabled");

        assert!(
            (actual.precession_angle - expected.precession_angle).abs() < 1e-15,
            "precession phase used the wrong inertia: actual={}, expected={}",
            actual.precession_angle,
            expected.precession_angle
        );
        assert!(
            (actual.nutation_phase - expected.nutation_phase).abs() < 1e-15,
            "nutation phase used the wrong inertia: actual={}, expected={}",
            actual.nutation_phase,
            expected.nutation_phase
        );
    }
}

#[cfg(test)]
mod form_factor_drag_tests {
    use super::*;

    fn acceleration_with_form_factor_flag(enabled: bool) -> Vector3<f64> {
        let inputs = BallisticInputs {
            bc_value: 0.462,
            bc_type: DragModel::G1,
            bullet_model: Some("168gr SMK Match".to_string()),
            use_form_factor: enabled,
            ..BallisticInputs::default()
        };
        let solver = TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        let (density, _, temp_c, pressure_hpa) = solver.resolved_atmosphere();
        solver.calculate_acceleration(
            &Vector3::zeros(),
            &Vector3::new(600.0, 0.0, 0.0),
            &Vector3::zeros(),
            (temp_c, pressure_hpa, density / 1.225),
        )
    }

    #[test]
    fn measured_bc_drag_does_not_apply_name_based_form_factor_again() {
        let baseline = acceleration_with_form_factor_flag(false);
        let flagged = acceleration_with_form_factor_flag(true);

        assert!(
            (flagged - baseline).norm() < 1e-12,
            "published BC already encodes form factor: baseline={baseline:?} flagged={flagged:?}"
        );
    }
}

#[cfg(test)]
mod rk45_adaptivity_tests {
    use super::*;

    #[test]
    fn cli_rk45_error_norm_scales_components_independently() {
        let position = Vector3::new(1.0e9, 0.0, 0.0);
        let velocity = Vector3::new(800.0, 0.0, 0.0);
        let fifth_position = position;
        let fifth_velocity = velocity;
        let fourth_position = position;
        let fourth_velocity = Vector3::new(800.0, 1.0e-3, 0.0);

        let error = cli_rk45_error_norm(
            &position,
            &velocity,
            &fifth_position,
            &fifth_velocity,
            &fourth_position,
            &fourth_velocity,
        );
        let expected = 1.0e-3 / 6.0_f64.sqrt();

        assert!(
            (error - expected).abs() <= 1e-15,
            "large downrange position masked a velocity-component error: {error}"
        );
    }

    fn discontinuous_wind_solver() -> TrajectorySolver {
        let inputs = BallisticInputs::default();
        let mut solver = TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        solver.set_wind_segments(vec![
            crate::wind::WindSegment::new(0.0, 90.0, 4.0),
            crate::wind::WindSegment::new(1_000.0, 90.0, 10_000.0),
        ]);
        solver
    }

    #[test]
    fn rk45_retries_discontinuous_trial_before_advancing() {
        let solver = discontinuous_wind_solver();
        let position = Vector3::new(0.0, solver.inputs.muzzle_height, 0.0);
        let velocity = Vector3::new(solver.inputs.muzzle_velocity, 0.0, 0.0);
        let (density, _, temp_c, pressure_hpa) = solver.resolved_atmosphere();
        let resolved_atmo = (temp_c, pressure_hpa, density / 1.225);
        let dt = 0.01;

        let rejected_trial = solver.rk45_step(
            &position,
            &velocity,
            dt,
            &Vector3::zeros(),
            RK45_TOLERANCE,
            resolved_atmo,
        );
        assert!(
            rejected_trial.error > RK45_TOLERANCE,
            "discontinuous full step must exceed tolerance, got {}",
            rejected_trial.error
        );

        let accepted = solver.adaptive_rk45_step(
            &position,
            &velocity,
            dt,
            &Vector3::zeros(),
            resolved_atmo,
        );
        assert!(accepted.used_dt < dt, "oversized trial was not retried");
        assert!(
            accepted.error <= RK45_TOLERANCE || accepted.used_dt <= RK45_MIN_DT,
            "accepted error {} exceeds tolerance at dt {}",
            accepted.error,
            accepted.used_dt
        );

        let accepted_trial = solver.rk45_step(
            &position,
            &velocity,
            accepted.used_dt,
            &Vector3::zeros(),
            RK45_TOLERANCE,
            resolved_atmo,
        );
        assert_eq!(accepted.position, accepted_trial.position);
        assert_eq!(accepted.velocity, accepted_trial.velocity);
        assert!((RK45_MIN_DT..=RK45_MAX_DT).contains(&accepted.next_dt));
    }
}

#[cfg(test)]
mod ground_termination_tests {
    use super::*;
    use crate::trajectory_observation::TrajectoryObservationFlag;

    #[test]
    fn every_solver_reports_one_exact_early_ground_endpoint() {
        for (name, use_rk4, use_adaptive_rk45) in [
            ("Euler", false, false),
            ("RK4", true, false),
            ("RK45", true, true),
        ] {
            let inputs = BallisticInputs {
                muzzle_height: 1.0,
                muzzle_angle: -0.2,
                ground_threshold: 0.0,
                use_rk4,
                use_adaptive_rk45,
                ..BallisticInputs::default()
            };
            let mut solver = TrajectorySolver::new(
                inputs,
                WindConditions::default(),
                AtmosphericConditions::default(),
            );
            solver.set_max_range(1_000.0);

            let result = solver.solve().expect("early-ground solve should succeed");
            let terminal = result.points.last().expect("terminal point is missing");

            assert_eq!(result.termination, TrajectoryTermination::GroundThreshold);
            assert_eq!(terminal.position.y.to_bits(), 0.0_f64.to_bits());
            assert!(
                terminal.position.x < 1_000.0,
                "{name} incorrectly reached max range"
            );
            assert_eq!(result.max_range.to_bits(), terminal.position.x.to_bits());
            assert_eq!(
                result
                    .points
                    .iter()
                    .filter(|point| point.position.y == 0.0)
                    .count(),
                1,
                "{name} did not retain exactly one ground endpoint"
            );

            let observations = result
                .sample_observations(1.0, 100)
                .expect("checked early-ground sampling should succeed");
            assert!(observations[..observations.len() - 1]
                .iter()
                .all(|observation| observation.distance_m < terminal.position.x));
            let terminal_observation = observations.last().expect("terminal observation");
            assert_eq!(
                terminal_observation.distance_m.to_bits(),
                terminal.position.x.to_bits()
            );
            assert!(terminal_observation
                .flags
                .contains(&TrajectoryObservationFlag::Terminal));
            assert!(terminal_observation
                .flags
                .contains(&TrajectoryObservationFlag::GroundThreshold));
            assert_eq!(
                observations
                    .iter()
                    .filter(|observation| observation
                        .flags
                        .contains(&TrajectoryObservationFlag::Terminal))
                    .count(),
                1,
                "{name} repeated the terminal observation"
            );
        }
    }

    // Regression lock for the unified ground termination: solve_euler/solve_rk4/solve_rk45 all
    // loop while `position.y > ground_threshold` (default -100.0), so they agree with RK45. A
    // lofted shot that returns to launch level before reaching max_range must keep descending to
    // the -100 m floor instead of stopping at y = 0 — and RK4-fixed and RK45 must behave the same.
    #[test]
    fn rk4_and_rk45_descend_to_ground_threshold() {
        for adaptive in [false, true] {
            let inputs = BallisticInputs {
                muzzle_angle: 0.1, // ~5.7 deg: arcs up, then descends past launch level
                use_rk4: true,
                use_adaptive_rk45: adaptive,
                ..BallisticInputs::default()
            };
            assert_eq!(
                inputs.ground_threshold, -100.0,
                "default ground_threshold is -100 m"
            );

            let mut solver = TrajectorySolver::new(
                inputs,
                WindConditions::default(),
                AtmosphericConditions::default(),
            );
            // Huge max range: termination must be driven by ground_threshold, not the range cap.
            solver.set_max_range(1.0e7);

            let result = solver.solve().expect("solve should succeed");
            let final_y = result
                .points
                .last()
                .expect("trajectory has points")
                .position
                .y;
            assert!(
                final_y < -1.0,
                "adaptive_rk45={adaptive}: final y = {final_y} m; a lofted shot should descend \
                 past launch level toward the ground_threshold floor, not stop at y = 0"
            );
        }
    }
}

#[cfg(test)]
mod magnus_stability_tests {
    use super::*;

    #[test]
    fn yaw_of_repose_magnus_force_is_vertical_and_twist_signed() {
        let acceleration = |enable_magnus, is_twist_right| {
            let inputs = BallisticInputs {
                muzzle_velocity: 822.96,
                bullet_mass: 168.0 * crate::constants::GRAINS_TO_KG,
                bullet_diameter: 0.308 * 0.0254,
                bullet_length: 1.215 * 0.0254,
                twist_rate: 10.0,
                is_twist_right,
                enable_magnus,
                ..BallisticInputs::default()
            };
            let solver = TrajectorySolver::new(
                inputs,
                WindConditions::default(),
                AtmosphericConditions::default(),
            );
            let (density, _, temp_c, pressure_hpa) = solver.resolved_atmosphere();
            solver.calculate_acceleration(
                &Vector3::zeros(),
                &Vector3::new(822.96, 0.0, 0.0),
                &Vector3::zeros(),
                (temp_c, pressure_hpa, density / 1.225),
            )
        };

        let baseline = acceleration(false, true);
        let right_twist = acceleration(true, true) - baseline;
        let left_twist = acceleration(true, false) - baseline;

        assert!(
            right_twist.y < 0.0,
            "right-hand Magnus must point down, got {right_twist:?}"
        );
        assert!(
            left_twist.y > 0.0,
            "left-hand Magnus must point up, got {left_twist:?}"
        );
        assert!((right_twist.y + left_twist.y).abs() < 1e-12);
        assert!(right_twist.x.abs() < 1e-12 && right_twist.z.abs() < 1e-12);
        assert!(left_twist.x.abs() < 1e-12 && left_twist.z.abs() < 1e-12);
    }

    #[test]
    fn magnus_uses_velocity_corrected_muzzle_stability_gate() {
        let muzzle_velocity = 1_400.0 / 3.28084;
        let inputs = BallisticInputs {
            muzzle_velocity,
            bullet_mass: 168.0 * crate::constants::GRAINS_TO_KG,
            bullet_diameter: 0.308 * 0.0254,
            bullet_length: 1.215 * 0.0254,
            twist_rate: 15.0,
            enable_magnus: true,
            ..BallisticInputs::default()
        };
        let solver = TrajectorySolver::new(
            inputs.clone(),
            WindConditions::default(),
            AtmosphericConditions::default(),
        );

        let bare_sg = crate::spin_drift::miller_stability(0.308, 168.0, 15.0, 1.215);
        let canonical_sg = solver.effective_spin_drift_sg();
        assert!(bare_sg > 1.0, "test requires bare Sg above the Magnus gate");
        assert!(
            canonical_sg < 1.0,
            "velocity-corrected Sg must be below the gate, got {canonical_sg}"
        );

        let (density, _, temp_c, pressure_hpa) = solver.resolved_atmosphere();
        let acceleration = solver.calculate_acceleration(
            &Vector3::zeros(),
            &Vector3::new(muzzle_velocity, 0.0, 0.0),
            &Vector3::zeros(),
            (temp_c, pressure_hpa, density / 1.225),
        );
        let mut baseline_inputs = inputs;
        baseline_inputs.enable_magnus = false;
        let baseline_solver = TrajectorySolver::new(
            baseline_inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        let baseline = baseline_solver.calculate_acceleration(
            &Vector3::zeros(),
            &Vector3::new(muzzle_velocity, 0.0, 0.0),
            &Vector3::zeros(),
            (temp_c, pressure_hpa, density / 1.225),
        );

        assert_eq!(
            acceleration, baseline,
            "canonical Sg below 1 must suppress every Magnus acceleration component"
        );
    }

    #[test]
    fn magnus_force_grows_as_fixed_spin_projectile_slows() {
        let inputs = BallisticInputs {
            muzzle_velocity: 800.0,
            bullet_mass: 168.0 * crate::constants::GRAINS_TO_KG,
            bullet_diameter: 0.308 * 0.0254,
            bullet_length: 1.215 * 0.0254,
            twist_rate: 12.0,
            enable_magnus: true,
            ..BallisticInputs::default()
        };

        let magnus_acceleration = |speed_mps| {
            let evaluate = |enable_magnus| {
                let mut run_inputs = inputs.clone();
                run_inputs.enable_magnus = enable_magnus;
                let solver = TrajectorySolver::new(
                    run_inputs,
                    WindConditions::default(),
                    AtmosphericConditions::default(),
                );
                let (density, _, temp_c, pressure_hpa) = solver.resolved_atmosphere();
                solver
                    .calculate_acceleration(
                        &Vector3::zeros(),
                        &Vector3::new(speed_mps, 0.0, 0.0),
                        &Vector3::zeros(),
                        (temp_c, pressure_hpa, density / 1.225),
                    )
                    .y
            };
            (evaluate(true) - evaluate(false)).abs()
        };

        let fast = magnus_acceleration(200.0);
        let slow = magnus_acceleration(100.0);
        let ratio = slow / fast;
        let expected_ratio = 2.0_f64.powf(5.0 / 3.0);

        assert!(fast > 0.0 && slow > 0.0, "fast={fast}, slow={slow}");
        assert!(
            (ratio - expected_ratio).abs() < 1e-3,
            "fixed-spin Magnus acceleration must grow downrange; slow/fast={ratio}, \
             expected={expected_ratio}"
        );
    }
}

#[cfg(test)]
mod coriolis_direction_tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    #[test]
    fn supersonic_crossing_flags_a_positive_range_sample() {
        // A supersonic shot that slows past Mach 1 must flag a sampled point as a Mach
        // transition. The underlying transonic_distances were a Vec::new() TODO, so this
        // flag was NEVER set regardless of trajectory — this is the regression guard.
        use crate::trajectory_sampling::TrajectoryFlag;

        for (solver_name, use_rk4, use_adaptive_rk45) in [
            ("Euler", false, false),
            ("RK4", true, false),
            ("RK45", true, true),
        ] {
            let inputs = BallisticInputs {
                muzzle_velocity: 850.0,
                bc_value: 0.2,
                bc_type: DragModel::G7,
                muzzle_angle: 0.03,
                enable_trajectory_sampling: true,
                sample_interval: 50.0,
                use_rk4,
                use_adaptive_rk45,
                ..BallisticInputs::default()
            };
            let mut solver = TrajectorySolver::new(
                inputs,
                WindConditions::default(),
                AtmosphericConditions::default(),
            );
            solver.set_max_range(2000.0);
            let samples = solver
                .solve()
                .expect("supersonic solve should succeed")
                .sampled_points
                .expect("sampling was enabled");
            let flagged_distances: Vec<_> = samples
                .iter()
                .filter(|sample| sample.flags.contains(&TrajectoryFlag::MachTransition))
                .map(|sample| sample.distance_m)
                .collect();

            assert!(
                !flagged_distances.is_empty()
                    && flagged_distances.iter().all(|distance| *distance > 0.0),
                "{solver_name} must flag genuine crossings only at positive range: {flagged_distances:?}"
            );
        }
    }

    #[test]
    fn subsonic_launch_does_not_flag_a_muzzle_transition() {
        use crate::trajectory_sampling::TrajectoryFlag;

        for (solver_name, use_rk4, use_adaptive_rk45) in [
            ("Euler", false, false),
            ("RK4", true, false),
            ("RK45", true, true),
        ] {
            let inputs = BallisticInputs {
                muzzle_velocity: 250.0,
                muzzle_angle: 0.02,
                enable_trajectory_sampling: true,
                sample_interval: 25.0,
                use_rk4,
                use_adaptive_rk45,
                ..BallisticInputs::default()
            };
            let mut solver = TrajectorySolver::new(
                inputs,
                WindConditions::default(),
                AtmosphericConditions::default(),
            );
            solver.set_max_range(300.0);
            let samples = solver
                .solve()
                .expect("subsonic solve should succeed")
                .sampled_points
                .expect("sampling was enabled");

            assert!(
                samples
                    .iter()
                    .all(|sample| !sample.flags.contains(&TrajectoryFlag::MachTransition)),
                "{solver_name} marked a Mach transition for a launch already below Mach 1"
            );
        }
    }

    #[test]
    fn mach_transition_tracker_requires_a_downward_crossing() {
        fn record(mach_values: &[f64]) -> Vec<f64> {
            let mut tracker = MachTransitionTracker::default();
            let mut distances = Vec::new();
            for (index, mach) in mach_values.iter().copied().enumerate() {
                tracker.record_downward_crossings(mach, index as f64 * 10.0, &mut distances);
            }
            distances
        }

        assert!(record(&[0.9, 0.8, 0.7]).is_empty());
        assert_eq!(record(&[1.1, 1.05, 0.99]), vec![20.0]);
        assert_eq!(record(&[1.2, 1.19, 1.0, 0.99]), vec![10.0, 30.0]);
        assert_eq!(record(&[0.9, 1.3, 1.1, 0.9, 1.3, 0.8]), vec![20.0, 30.0]);
        assert!(record(&[1.3, f64::NAN, 1.1]).is_empty());
    }

    #[test]
    fn mach_transition_tracker_labels_0_9_without_touching_the_flat_vec() {
        // MBA-1405: the tracker gains a third (0.9) threshold, but it must NEVER be appended
        // to the flat `distances` Vec — only its own labeled field. Reuses the exact input
        // sequences pinned by `mach_transition_tracker_requires_a_downward_crossing` above (the
        // pre-change flat-Vec behavior for those sequences) and additionally asserts the
        // labeled fields.
        fn record(mach_values: &[f64]) -> (Vec<f64>, MachTransitionTracker) {
            let mut tracker = MachTransitionTracker::default();
            let mut distances = Vec::new();
            for (index, mach) in mach_values.iter().copied().enumerate() {
                tracker.record_downward_crossings(mach, index as f64 * 10.0, &mut distances);
            }
            (distances, tracker)
        }

        // Already at/below 1.2 and 1.0 at the muzzle (never crosses those downward from
        // above), but DOES cross 0.9 downward at x=10 (previous_mach 0.9 >= 0.9, mach 0.8 < 0.9).
        let (distances, tracker) = record(&[0.9, 0.8, 0.7]);
        assert!(distances.is_empty()); // flat Vec unchanged (pinned pre-change value)
        assert_eq!(tracker.mach_1_2_distance_m, None);
        assert_eq!(tracker.mach_1_0_distance_m, None);
        assert_eq!(tracker.mach_0_9_distance_m, Some(10.0));

        // Crosses 1.0 only (never reaches 1.2 from above, never reaches 0.9).
        let (distances, tracker) = record(&[1.1, 1.05, 0.99]);
        assert_eq!(distances, vec![20.0]);
        assert_eq!(tracker.mach_1_2_distance_m, None);
        assert_eq!(tracker.mach_1_0_distance_m, Some(20.0));
        assert_eq!(tracker.mach_0_9_distance_m, None);

        // Crosses 1.2 then 1.0, never reaches 0.9 (lowest sample is 0.99).
        let (distances, tracker) = record(&[1.2, 1.19, 1.0, 0.99]);
        assert_eq!(distances, vec![10.0, 30.0]); // flat Vec unchanged
        assert_eq!(tracker.mach_1_2_distance_m, Some(10.0));
        assert_eq!(tracker.mach_1_0_distance_m, Some(30.0));
        assert_eq!(tracker.mach_0_9_distance_m, None);

        // Crosses 1.2, then 1.0, then 0.9 — the flat Vec must be EXACTLY {20.0, 30.0} as
        // before (0.9 crossing at x=50.0 must NOT appear in it).
        let (distances, tracker) = record(&[0.9, 1.3, 1.1, 0.9, 1.3, 0.8]);
        assert_eq!(distances, vec![20.0, 30.0]); // flat Vec unchanged (pinned pre-change value)
        assert_eq!(tracker.mach_1_2_distance_m, Some(20.0));
        assert_eq!(tracker.mach_1_0_distance_m, Some(30.0));
        assert_eq!(tracker.mach_0_9_distance_m, Some(50.0));
        assert!(
            tracker.mach_1_2_distance_m < tracker.mach_1_0_distance_m
                && tracker.mach_1_0_distance_m < tracker.mach_0_9_distance_m,
            "labeled crossings must be strictly increasing downrange"
        );

        // NAN resets tracking exactly as before; no labels set either.
        let (distances, tracker) = record(&[1.3, f64::NAN, 1.1]);
        assert!(distances.is_empty());
        assert_eq!(tracker.mach_1_2_distance_m, None);
        assert_eq!(tracker.mach_1_0_distance_m, None);
        assert_eq!(tracker.mach_0_9_distance_m, None);
    }

    #[test]
    fn humidity_percent_converts_and_clamps() {
        // MBA-722: BallisticInputs.humidity is a 0-1 fraction; the helper yields 0-100 percent.
        let mut i = BallisticInputs {
            humidity: 0.5,
            ..BallisticInputs::default()
        };
        assert!((i.humidity_percent() - 50.0).abs() < 1e-9, "0.5 -> 50%");
        i.humidity = 0.0;
        assert_eq!(i.humidity_percent(), 0.0);
        i.humidity = 1.0;
        assert_eq!(i.humidity_percent(), 100.0);
        i.humidity = 1.5; // out of range -> clamped, never > 100
        assert_eq!(i.humidity_percent(), 100.0);
    }

    /// Vertical position (m) at a given downrange `range_m`, for a shot fired along
    /// compass bearing `shot_azimuth` (radians, 0=N) with Coriolis enabled.
    fn vertical_at(shot_azimuth: f64, range_m: f64) -> f64 {
        let inputs = BallisticInputs {
            muzzle_velocity: 800.0,
            bc_value: 0.5,
            bc_type: DragModel::G7,
            muzzle_angle: 0.02, // ~20 mrad so it carries well past range_m
            enable_coriolis: true,
            latitude: Some(45.0),
            shot_azimuth,
            ground_threshold: f64::NEG_INFINITY, // never terminate early
            ..BallisticInputs::default()
        };
        let mut solver = TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        solver.set_max_range(range_m + 50.0);
        let r = solver.solve().expect("solve");
        let pts = &r.points;
        for i in 1..pts.len() {
            if pts[i].position.x >= range_m {
                let p1 = &pts[i - 1];
                let p2 = &pts[i];
                let t = (range_m - p1.position.x) / (p2.position.x - p1.position.x);
                return p1.position.y + t * (p2.position.y - p1.position.y);
            }
        }
        panic!("range {range_m} not reached");
    }

    /// Regression for the shot-direction Coriolis bug: the Eötvös vertical term
    /// `a_up = +2Ω cosφ v_east` lifts an EAST shot and depresses a WEST shot, so at a
    /// common range east must sit HIGHER than west, with north in between. Before the
    /// fix, `--shot-direction` never reached the solver and E/W/N were identical.
    #[test]
    fn eotvos_east_higher_than_west() {
        let range = 600.0;
        let east = vertical_at(FRAC_PI_2, range); // 90° E
        let west = vertical_at(3.0 * FRAC_PI_2, range); // 270° W
        let north = vertical_at(0.0, range); // 0° N
        assert!(
            east > west,
            "east ({east:.5}) must be higher than west ({west:.5}) at {range} m (Eötvös)"
        );
        assert!(
            east > north && north > west,
            "north ({north:.5}) must lie between east ({east:.5}) and west ({west:.5})"
        );
        assert!(
            (east - west) > 1e-3,
            "E-W vertical separation ({:.6} m) should be physically meaningful, not FP noise",
            east - west
        );
    }

    /// MBA-1405: labeled Mach-crossing distances on `TrajectoryResult`, pinned exactly against
    /// a capture of the (pre-change) flat `transonic_distances` Vec taken from this same
    /// muzzle_velocity=850/bc=0.2 G7/angle=0.03/max_range=2000 scenario, for all three solver
    /// paths. Because `MachTransitionTracker::record_downward_crossings` writes the labeled
    /// field and pushes to the flat Vec from the SAME `downrange_m` value at the SAME crossing,
    /// these pinned labels ARE the flat-Vec content for the 1.2/1.0 thresholds — proving the
    /// flat Vec is untouched by the 0.9 addition.
    #[test]
    fn labeled_mach_crossings_match_pinned_pre_change_flat_vec_across_solvers() {
        // (solver_name, use_rk4, use_adaptive_rk45, expected mach_1_2, expected mach_1_0)
        let cases = [
            ("Euler", false, false, 670.9878683238721_f64, 805.5274119916264_f64),
            ("RK4", true, false, 671.7257336844475_f64, 805.933409072171_f64),
            ("RK45", true, true, 672.4905711917901_f64, 806.5709746782849_f64),
        ];

        for (solver_name, use_rk4, use_adaptive_rk45, expected_1_2, expected_1_0) in cases {
            let inputs = BallisticInputs {
                muzzle_velocity: 850.0,
                bc_value: 0.2,
                bc_type: DragModel::G7,
                muzzle_angle: 0.03,
                use_rk4,
                use_adaptive_rk45,
                ..BallisticInputs::default()
            };
            let mut solver = TrajectorySolver::new(
                inputs,
                WindConditions::default(),
                AtmosphericConditions::default(),
            );
            solver.set_max_range(2000.0);
            let result = solver.solve().expect("solve should succeed");

            assert_eq!(
                result.mach_1_2_distance_m,
                Some(expected_1_2),
                "{solver_name}: mach_1_2_distance_m must match the pinned pre-change flat-Vec value"
            );
            assert_eq!(
                result.mach_1_0_distance_m,
                Some(expected_1_0),
                "{solver_name}: mach_1_0_distance_m must match the pinned pre-change flat-Vec value"
            );

            let mach_1_2 = result.mach_1_2_distance_m.expect("crosses 1.2");
            let mach_1_0 = result.mach_1_0_distance_m.expect("crosses 1.0");
            let mach_0_9 = result
                .mach_0_9_distance_m
                .expect("this trajectory also goes past 0.9 within 2000 m");
            assert!(
                mach_1_2 < mach_1_0 && mach_1_0 < mach_0_9,
                "{solver_name}: labeled crossings must be strictly increasing downrange \
                 (1.2={mach_1_2}, 1.0={mach_1_0}, 0.9={mach_0_9})"
            );
        }
    }

    /// A trajectory that terminates while still well above Mach 1.2 must leave all three
    /// labeled crossings `None` (there is nothing to cross), across all three solver paths.
    #[test]
    fn labeled_mach_crossings_are_none_for_a_fully_supersonic_trajectory() {
        for (solver_name, use_rk4, use_adaptive_rk45) in [
            ("Euler", false, false),
            ("RK4", true, false),
            ("RK45", true, true),
        ] {
            let inputs = BallisticInputs {
                muzzle_velocity: 850.0,
                bc_value: 0.2,
                bc_type: DragModel::G7,
                muzzle_angle: 0.03,
                use_rk4,
                use_adaptive_rk45,
                ..BallisticInputs::default()
            };
            let mut solver = TrajectorySolver::new(
                inputs,
                WindConditions::default(),
                AtmosphericConditions::default(),
            );
            // Well short of the ~671 m downward 1.2 crossing measured for this load.
            solver.set_max_range(200.0);
            let result = solver.solve().expect("solve should succeed");

            assert_eq!(
                result.mach_1_2_distance_m, None,
                "{solver_name}: must not report a 1.2 crossing that never happens"
            );
            assert_eq!(
                result.mach_1_0_distance_m, None,
                "{solver_name}: must not report a 1.0 crossing that never happens"
            );
            assert_eq!(
                result.mach_0_9_distance_m, None,
                "{solver_name}: must not report a 0.9 crossing that never happens"
            );
        }
    }
}

#[cfg(test)]
mod cant_tests {
    use super::*;

    fn base_inputs() -> BallisticInputs {
        BallisticInputs {
            muzzle_velocity: 800.0,
            bc_value: 0.5,
            bc_type: DragModel::G7,
            bullet_mass: 0.0109,
            bullet_diameter: 0.00782,
            bullet_length: 0.0309,
            sight_height: 0.05,
            twist_rate: 10.0,
            use_rk4: true,
            ..BallisticInputs::default()
        }
    }

    fn solve_with(inputs: BallisticInputs, max_range: f64) -> TrajectoryResult {
        let mut s = TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        s.set_max_range(max_range);
        s.solve().expect("solve")
    }

    /// Interpolate (y, z) at downrange x.
    fn yz_at(result: &TrajectoryResult, x: f64) -> (f64, f64) {
        let pts = &result.points;
        for i in 1..pts.len() {
            if pts[i].position.x >= x {
                let (p1, p2) = (&pts[i - 1], &pts[i]);
                let dx = p2.position.x - p1.position.x;
                let t = if dx.abs() < 1e-12 { 0.0 } else { (x - p1.position.x) / dx };
                return (
                    p1.position.y + t * (p2.position.y - p1.position.y),
                    p1.position.z + t * (p2.position.z - p1.position.z),
                );
            }
        }
        panic!("trajectory never reached {x} m");
    }

    #[test]
    fn cant_sign_clockwise_up_offset_goes_right_and_low() {
        // Upward zero offset + clockwise cant => POI right (+z) and low vs un-canted.
        let mut level = base_inputs();
        level.muzzle_angle = 0.003; // ~10 MOA up
        let mut canted = level.clone();
        canted.cant_angle = 10f64.to_radians();

        let (y0, z0) = yz_at(&solve_with(level, 400.0), 300.0);
        let (y1, z1) = yz_at(&solve_with(canted, 400.0), 300.0);
        assert!(z1 > z0 + 0.01, "clockwise cant must move POI right: z0={z0} z1={z1}");
        assert!(y1 < y0 - 0.001, "clockwise cant must move POI low: y0={y0} y1={y1}");
    }

    #[test]
    fn pure_cant_shows_bore_offset_near_range() {
        // No aim offset: the only lateral source near the muzzle is the swung bore,
        // z0 = -sight_height*sin(cant) (left of the aim plane for clockwise cant).
        let mut i = base_inputs();
        i.muzzle_angle = 0.0;
        i.cant_angle = 10f64.to_radians();
        let sh = i.sight_height;
        let r = solve_with(i, 60.0);
        let first = &r.points[1]; // just past the muzzle
        let expected = -sh * 10f64.to_radians().sin();
        assert!(
            (first.position.z - expected).abs() < 0.005,
            "near-muzzle lateral {} should be ~bore offset {expected}",
            first.position.z
        );
    }

    #[test]
    fn zero_angle_is_independent_of_cant() {
        let a = base_inputs();
        let mut b = base_inputs();
        b.cant_angle = 15f64.to_radians();
        let za = calculate_zero_angle(a.clone(), 100.0, 0.0).expect("zero a");
        let zb = calculate_zero_angle(b.clone(), 100.0, 0.0).expect("zero b");
        assert_eq!(za.to_bits(), zb.to_bits(), "zeroing must ignore cant: {za} vs {zb}");
        // silence unused warnings
        let _ = (a.cant_angle, b.cant_angle);
    }

    #[test]
    fn nonfinite_cant_is_rejected() {
        let mut i = base_inputs();
        i.cant_angle = f64::NAN;
        let s = TrajectorySolver::new(i, WindConditions::default(), AtmosphericConditions::default());
        assert!(s.solve().is_err());
    }

    #[test]
    fn incline_and_cant_compose_without_breaking() {
        // 15-degree incline + 10-degree cant: finite result, cant still pushes right.
        let mut flat = base_inputs();
        flat.muzzle_angle = 0.003;
        flat.shooting_angle = 15f64.to_radians();
        let mut canted = flat.clone();
        canted.cant_angle = 10f64.to_radians();
        let (_, z_flat) = yz_at(&solve_with(flat, 400.0), 300.0);
        let (_, z_cant) = yz_at(&solve_with(canted, 400.0), 300.0);
        assert!(z_cant > z_flat, "cant must still deflect right on an incline");
    }
}

#[cfg(test)]
mod vertical_wind_tests {
    use super::*;

    fn base_inputs() -> BallisticInputs {
        BallisticInputs {
            muzzle_velocity: 800.0,
            bc_value: 0.5,
            bc_type: DragModel::G7,
            bullet_mass: 0.0109,
            bullet_diameter: 0.00782,
            bullet_length: 0.0309,
            sight_height: 0.05,
            twist_rate: 10.0,
            use_rk4: true,
            ..BallisticInputs::default()
        }
    }

    /// Interpolate trajectory height (McCoy Y) at downrange distance `x`.
    fn y_at(result: &TrajectoryResult, x: f64) -> f64 {
        let pts = &result.points;
        for i in 1..pts.len() {
            if pts[i].position.x >= x {
                let (p1, p2) = (&pts[i - 1], &pts[i]);
                let dx = p2.position.x - p1.position.x;
                let t = if dx.abs() < 1e-12 { 0.0 } else { (x - p1.position.x) / dx };
                return p1.position.y + t * (p2.position.y - p1.position.y);
            }
        }
        panic!("trajectory never reached {x} m");
    }

    fn solve_with(inputs: BallisticInputs, wind: WindConditions, max_range: f64) -> TrajectoryResult {
        let mut s = TrajectorySolver::new(inputs, wind, AtmosphericConditions::default());
        s.set_max_range(max_range);
        s.solve().expect("solve")
    }

    #[test]
    fn updraft_raises_poi_downrange() {
        // No shear, no segments: this exercises the constant-wind sites in
        // solve_euler/solve_rk4/solve_rk45 directly (MBA-728).
        let calm_inputs = base_inputs();
        let calm_wind = WindConditions::default();
        let updraft = WindConditions {
            vertical_speed: 5.0,
            ..Default::default()
        };

        let calm = solve_with(calm_inputs.clone(), calm_wind, 500.0);
        let updraft_result = solve_with(calm_inputs, updraft, 500.0);

        let y_calm = y_at(&calm, 400.0);
        let y_updraft = y_at(&updraft_result, 400.0);
        assert!(
            y_updraft > y_calm,
            "5 m/s updraft must raise POI at 400m: calm={y_calm}, updraft={y_updraft}"
        );
    }

    #[test]
    fn zero_vertical_is_default_and_finite_required() {
        assert_eq!(WindConditions::default().vertical_speed, 0.0);

        let inputs = base_inputs();
        let wind = WindConditions {
            vertical_speed: f64::NAN,
            ..Default::default()
        };
        let s = TrajectorySolver::new(inputs, wind, AtmosphericConditions::default());
        assert!(
            s.solve().is_err(),
            "NaN wind.vertical_speed must be rejected by validate_for_solve"
        );
    }
}

/// MBA-1365: declare whether a BC is ICAO- or Army-Standard-Metro-referenced.
#[cfg(test)]
mod bc_reference_standard_tests {
    use super::*;

    fn base_inputs() -> BallisticInputs {
        BallisticInputs {
            muzzle_velocity: 800.0,
            bc_value: 0.5,
            bc_type: DragModel::G7,
            bullet_mass: 0.0109,
            bullet_diameter: 0.00782,
            bullet_length: 0.0309,
            sight_height: 0.05,
            twist_rate: 10.0,
            use_rk4: true,
            ..BallisticInputs::default()
        }
    }

    /// Interpolate trajectory height (McCoy Y) and speed at downrange distance `x`.
    fn y_and_speed_at(result: &TrajectoryResult, x: f64) -> (f64, f64) {
        let pts = &result.points;
        for i in 1..pts.len() {
            if pts[i].position.x >= x {
                let (p1, p2) = (&pts[i - 1], &pts[i]);
                let dx = p2.position.x - p1.position.x;
                let t = if dx.abs() < 1e-12 {
                    0.0
                } else {
                    (x - p1.position.x) / dx
                };
                return (
                    p1.position.y + t * (p2.position.y - p1.position.y),
                    p1.velocity_magnitude + t * (p2.velocity_magnitude - p1.velocity_magnitude),
                );
            }
        }
        panic!("trajectory never reached {x} m");
    }

    // ---- The constant itself -----------------------------------------------------------

    #[test]
    fn asm_to_icao_ratio_matches_documented_value() {
        assert!(
            (crate::constants::ASM_TO_ICAO_BC - 0.98237).abs() < 1e-5,
            "ASM_TO_ICAO_BC = {} must equal 0.98237 to 5 decimal places",
            crate::constants::ASM_TO_ICAO_BC
        );
        // Derived from the two named densities, not a bare literal. (The ratio being < 1.0
        // — a denser reference atmosphere implying a numerically smaller ICAO-equivalent
        // BC — is already established by the 0.98237 check above; a bare `assert!` on a
        // const-only comparison here would be compiled away, so it isn't repeated.)
        assert_eq!(
            crate::constants::ASM_TO_ICAO_BC,
            crate::constants::ASM_DENSITY_LB_FT3 / crate::constants::ICAO_DENSITY_LB_FT3
        );
    }

    // ---- Default / byte-identical --------------------------------------------------------

    #[test]
    fn default_bc_reference_standard_is_icao() {
        assert_eq!(
            BallisticInputs::default().bc_reference_standard,
            BcReferenceStandard::Icao
        );
    }

    /// RED (pre-MBA-1365) would not compile at all — there was no field to normalize.
    /// GREEN: with the field defaulted to `Icao`, `TrajectorySolver::new` must never
    /// touch `bc_value` — bit-for-bit, not just approximately — so every existing
    /// caller that never sets this field is byte-identical to today.
    #[test]
    fn icao_reference_leaves_bc_value_bit_identical() {
        let raw_bc: f64 = 0.4372911; // an arbitrary, non-round value
        let inputs = BallisticInputs {
            bc_value: raw_bc,
            bc_reference_standard: BcReferenceStandard::Icao,
            ..base_inputs()
        };
        let solver = TrajectorySolver::new(inputs, WindConditions::default(), AtmosphericConditions::default());
        assert_eq!(solver.inputs.bc_value.to_bits(), raw_bc.to_bits());
    }

    #[test]
    fn default_inputs_solve_is_unaffected_by_the_new_field_existing() {
        // A solve built entirely from BallisticInputs::default() (which now carries
        // bc_reference_standard: Icao) must match a solve of an equivalent struct that
        // never mentions the field at all in its literal (relying on ..default()).
        let a = TrajectorySolver::new(base_inputs(), WindConditions::default(), AtmosphericConditions::default())
            .solve()
            .expect("solve a");
        let b = TrajectorySolver::new(
            BallisticInputs { ..base_inputs() },
            WindConditions::default(),
            AtmosphericConditions::default(),
        )
        .solve()
        .expect("solve b");
        assert_eq!(a.impact_velocity.to_bits(), b.impact_velocity.to_bits());
        assert_eq!(a.max_range.to_bits(), b.max_range.to_bits());
    }

    // ---- ArmyStandardMetro normalization --------------------------------------------------

    #[test]
    fn army_standard_metro_scales_bc_value_by_exactly_the_derived_ratio() {
        let raw_bc = 0.5;
        let inputs = BallisticInputs {
            bc_value: raw_bc,
            bc_reference_standard: BcReferenceStandard::ArmyStandardMetro,
            ..base_inputs()
        };
        let solver = TrajectorySolver::new(inputs, WindConditions::default(), AtmosphericConditions::default());
        assert_eq!(
            solver.inputs.bc_value,
            raw_bc * crate::constants::ASM_TO_ICAO_BC
        );
    }

    #[test]
    fn army_standard_metro_scales_mach_keyed_bc_segments() {
        let inputs = BallisticInputs {
            bc_reference_standard: BcReferenceStandard::ArmyStandardMetro,
            bc_segments: Some(vec![(0.5, 0.40), (1.5, 0.30)]),
            ..base_inputs()
        };
        let solver = TrajectorySolver::new(inputs, WindConditions::default(), AtmosphericConditions::default());
        let segments = solver.inputs.bc_segments.as_ref().expect("segments");
        assert_eq!(segments[0], (0.5, 0.40 * crate::constants::ASM_TO_ICAO_BC));
        assert_eq!(segments[1], (1.5, 0.30 * crate::constants::ASM_TO_ICAO_BC));
    }

    #[test]
    fn army_standard_metro_scales_velocity_keyed_bc_segments_data() {
        let inputs = BallisticInputs {
            bc_reference_standard: BcReferenceStandard::ArmyStandardMetro,
            bc_segments_data: Some(vec![
                crate::BCSegmentData {
                    velocity_min: 0.0,
                    velocity_max: 500.0,
                    bc_value: 0.40,
                },
                crate::BCSegmentData {
                    velocity_min: 500.0,
                    velocity_max: 900.0,
                    bc_value: 0.45,
                },
            ]),
            ..base_inputs()
        };
        let solver = TrajectorySolver::new(inputs, WindConditions::default(), AtmosphericConditions::default());
        let segments = solver.inputs.bc_segments_data.as_ref().expect("segments");
        assert_eq!(segments[0].bc_value, 0.40 * crate::constants::ASM_TO_ICAO_BC);
        assert_eq!(segments[1].bc_value, 0.45 * crate::constants::ASM_TO_ICAO_BC);
        // Non-BC fields must be untouched.
        assert_eq!(segments[0].velocity_min, 0.0);
        assert_eq!(segments[1].velocity_max, 900.0);
    }

    /// Empirically verify the drag DIRECTION (not just the arithmetic): declaring the
    /// SAME raw BC number as Army-Standard-Metro must produce MORE drop and LOWER
    /// remaining velocity downrange than declaring it ICAO, because the normalized
    /// (smaller) BC feeds MORE drag into the ICAO-calibrated retardation formula.
    #[test]
    fn army_standard_metro_moves_impact_in_the_more_drag_direction() {
        let solve_at = |standard: BcReferenceStandard| {
            let inputs = BallisticInputs {
                bc_value: 0.475,
                bc_reference_standard: standard,
                ..base_inputs()
            };
            let mut solver =
                TrajectorySolver::new(inputs, WindConditions::default(), AtmosphericConditions::default());
            solver.set_max_range(500.0);
            solver.solve().expect("solve")
        };

        let icao = solve_at(BcReferenceStandard::Icao);
        let asm = solve_at(BcReferenceStandard::ArmyStandardMetro);

        let (y_icao, v_icao) = y_and_speed_at(&icao, 400.0);
        let (y_asm, v_asm) = y_and_speed_at(&asm, 400.0);

        assert!(
            y_asm < y_icao,
            "ArmyStandardMetro must drop MORE (lower y) at 400m than Icao for the same raw \
             bc_value: icao_y={y_icao}, asm_y={y_asm}"
        );
        assert!(
            v_asm < v_icao,
            "ArmyStandardMetro must retain LESS velocity at 400m than Icao for the same raw \
             bc_value: icao_v={v_icao}, asm_v={v_asm}"
        );
    }

    // ---- Downstream inheritance: Monte Carlo ----------------------------------------------

    /// Monte Carlo must inherit the exact same normalization with zero duplicated math:
    /// with every std-dev pinned to 0 (so every sample is deterministically the base
    /// input, unperturbed), a single-sample MC run must match a plain direct solve of the
    /// identical base inputs bit-for-bit.
    #[test]
    fn monte_carlo_inherits_the_normalized_bc_reference() {
        let base_inputs_asm = BallisticInputs {
            bc_value: 0.475,
            bc_reference_standard: BcReferenceStandard::ArmyStandardMetro,
            ..base_inputs()
        };
        let wind = WindConditions::default();

        // Match run_monte_carlo_with_wind_and_direction_std_dev_using_rng's own
        // solver_max_range derivation exactly (target_hint.max(1000.0) * 2.0, with
        // target_hint = base_inputs.target_distance since MonteCarloParams.target_distance
        // is None below) — otherwise the two solves integrate to different caps and a
        // level, non-ground-impacting shot legitimately reports a different max_range.
        let mut direct_solver =
            TrajectorySolver::new(base_inputs_asm.clone(), wind.clone(), AtmosphericConditions::default());
        direct_solver.set_max_range(base_inputs_asm.target_distance.max(1000.0) * 2.0);
        let direct = direct_solver.solve().expect("direct solve");

        let mc_params = MonteCarloParams {
            num_simulations: 1,
            velocity_std_dev: 0.0,
            angle_std_dev: 0.0,
            bc_std_dev: 0.0,
            wind_speed_std_dev: 0.0,
            target_distance: None,
            base_wind_speed: 0.0,
            base_wind_direction: 0.0,
            azimuth_std_dev: 0.0,
        };
        let mc = run_monte_carlo_with_wind_and_direction_std_dev_seeded(
            base_inputs_asm,
            wind,
            mc_params,
            0.0,
            42,
        )
        .expect("monte carlo");

        assert_eq!(mc.ranges.len(), 1);
        assert_eq!(
            mc.ranges[0].to_bits(),
            direct.max_range.to_bits(),
            "a zero-dispersion single MC sample must match a plain solve of the same \
             ASM-referenced inputs bit-for-bit"
        );
        assert_eq!(
            mc.impact_velocities[0].to_bits(),
            direct.impact_velocity.to_bits()
        );
    }

    // ---- Downstream inheritance: estimate-bc fitting --------------------------------------

    /// `estimate_bc_fit` always builds its internal search trials from
    /// `BallisticInputs { ..Default::default() }`, which defaults to `Icao` — so a fitted
    /// BC is always ICAO-referenced regardless of any other configuration in the process.
    /// Prove it by generating synthetic drop data from a KNOWN Icao-referenced bc_value and
    /// confirming the fit recovers that same value (not an ASM-shifted one).
    #[test]
    fn estimate_bc_fit_recovers_an_icao_referenced_bc() {
        let known_bc = 0.475;
        let velocity = 800.0;
        let mass = 0.0109;
        let diameter = 0.00782;
        let atmosphere = AtmosphericConditions::default();

        let synth_inputs = BallisticInputs {
            muzzle_velocity: velocity,
            bc_value: known_bc,
            bc_type: DragModel::G7,
            bullet_mass: mass,
            bullet_diameter: diameter,
            bullet_length: 0.0309,
            sight_height: 0.05,
            twist_rate: 10.0,
            use_rk4: true,
            bc_reference_standard: BcReferenceStandard::Icao,
            ..BallisticInputs::default()
        };
        let mut solver = TrajectorySolver::new(synth_inputs, WindConditions::default(), atmosphere.clone());
        solver.set_max_range(500.0);
        let trajectory = solver.solve().expect("synthetic solve");

        let points: Vec<(f64, f64)> = [100.0, 200.0, 300.0, 400.0]
            .iter()
            .map(|&d| {
                let (y, _) = {
                    let pts = &trajectory.points;
                    let mut found = None;
                    for i in 1..pts.len() {
                        if pts[i].position.x >= d {
                            let (p1, p2) = (&pts[i - 1], &pts[i]);
                            let dx = p2.position.x - p1.position.x;
                            let t = if dx.abs() < 1e-12 {
                                0.0
                            } else {
                                (d - p1.position.x) / dx
                            };
                            found = Some((
                                p1.position.y + t * (p2.position.y - p1.position.y),
                                0.0,
                            ));
                            break;
                        }
                    }
                    found.expect("trajectory reached observation distance")
                };
                (d, -y) // drop is positive-down
            })
            .collect();

        let estimate = estimate_bc_fit(
            velocity,
            mass,
            diameter,
            &points,
            DragModel::G7,
            BcFitMode::Drop,
            atmosphere,
            None,
            0.05,
        )
        .expect("fit should converge");

        assert!(
            (estimate.bc - known_bc).abs() < 0.02,
            "fit should recover the known ICAO-referenced bc={known_bc}, got {}",
            estimate.bc
        );
    }

    // ---- Custom drag table: documented inert behavior -------------------------------------

    #[test]
    fn custom_drag_table_makes_bc_reference_standard_numerically_inert() {
        let table = crate::drag::DragTable::try_new(vec![0.5, 1.0, 2.0, 3.0], vec![0.3, 0.4, 0.3, 0.2])
            .expect("valid table");

        let solve_with = |standard: BcReferenceStandard| {
            let inputs = BallisticInputs {
                bc_value: 0.5, // irrelevant once a custom drag table is active
                bc_reference_standard: standard,
                custom_drag_table: Some(table.clone()),
                ..base_inputs()
            };
            let mut solver =
                TrajectorySolver::new(inputs, WindConditions::default(), AtmosphericConditions::default());
            solver.set_max_range(500.0);
            solver.solve().expect("solve")
        };

        let icao = solve_with(BcReferenceStandard::Icao);
        let asm = solve_with(BcReferenceStandard::ArmyStandardMetro);

        assert_eq!(
            icao.impact_velocity.to_bits(),
            asm.impact_velocity.to_bits(),
            "a custom drag table must make bc_reference_standard fully inert"
        );
        assert_eq!(icao.max_range.to_bits(), asm.max_range.to_bits());
    }

    #[test]
    fn custom_drag_table_inert_warning_fires_only_for_army_standard_metro_with_a_table() {
        let table = crate::drag::DragTable::try_new(vec![0.5, 1.0, 2.0], vec![0.3, 0.4, 0.3])
            .expect("valid table");

        // No table: never warns, regardless of the declared standard.
        let no_table_icao = base_inputs();
        assert!(no_table_icao.bc_reference_standard_inert_warning().is_none());
        let no_table_asm = BallisticInputs {
            bc_reference_standard: BcReferenceStandard::ArmyStandardMetro,
            ..base_inputs()
        };
        assert!(no_table_asm.bc_reference_standard_inert_warning().is_none());

        // Table + Icao: no conversion was ever going to apply, so no warning either.
        let table_icao = BallisticInputs {
            custom_drag_table: Some(table.clone()),
            ..base_inputs()
        };
        assert!(table_icao.bc_reference_standard_inert_warning().is_none());

        // Table + ArmyStandardMetro: the one case that is actually inert -> warn.
        let table_asm = BallisticInputs {
            custom_drag_table: Some(table),
            bc_reference_standard: BcReferenceStandard::ArmyStandardMetro,
            ..base_inputs()
        };
        let warning = table_asm
            .bc_reference_standard_inert_warning()
            .expect("must warn");
        assert!(warning.contains("--bc-reference"));
        assert!(warning.contains("--drag-table"));
    }
}

/// MBA-1423: the per-point effective drag coefficient.
///
/// These pin the DEFINITION (the projectile's own Cd, form-factor scaled) rather than the
/// plumbing, because the plumbing being present is not the risk — reporting a plausible number
/// that is not the one flown is.
#[cfg(test)]
mod effective_drag_coefficient_tests {
    use super::*;

    fn inputs_175gr_g7() -> BallisticInputs {
        let mut inputs = BallisticInputs {
            bc_value: 0.243,
            bc_type: DragModel::G7,
            muzzle_velocity: 823.0,
            ..Default::default()
        };
        // Set the SI fields: TrajectorySolver::new re-derives caliber_inches/weight_grains
        // from bullet_diameter/bullet_mass, so setting only the imperial pair would be
        // silently overwritten by the SI defaults (a 154gr .30, not a 175gr .308).
        inputs.bullet_mass = 175.0 * crate::constants::GRAINS_TO_KG;
        inputs.bullet_diameter = 0.308 * 0.0254;
        inputs.weight_grains = 175.0;
        inputs.caliber_inches = 0.308;
        inputs
    }

    fn solver(inputs: BallisticInputs) -> TrajectorySolver {
        TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        )
    }

    /// The reported Cd is the projectile's own, i.e. the reference Cd scaled by the form factor
    /// SD/BC — not the reference table value itself, which would be identical for every bullet
    /// sharing a drag model and so useless for charting a specific load.
    #[test]
    fn reports_the_projectiles_own_cd_not_the_reference_tables() {
        let inputs = inputs_175gr_g7();
        let sd = inputs.sectional_density_lb_in2().expect("SD");
        let solver = solver(inputs);

        let sos = 340.0;
        let velocity = 800.0;
        let mach = velocity / sos;

        let reference = crate::drag::get_drag_coefficient(mach, &DragModel::G7);
        let reported = solver
            .effective_drag_coefficient(velocity, sos)
            .expect("mass and diameter are set");

        let expected = reference * sd / 0.243;
        assert!(
            (reported - expected).abs() < 1e-12,
            "reported {reported} != Cd_ref * SD / BC {expected}"
        );
        // The form factor here is > 1, so the two genuinely differ: a test that passed with
        // `reported == reference` would not be testing anything.
        assert!(
            (reported - reference).abs() > 1e-6,
            "form factor collapsed to 1; this fixture no longer distinguishes the two values"
        );
    }

    /// A custom drag table already supplies the projectile's actual Cd and divides by sectional
    /// density, so the form-factor scale must collapse to exactly 1 and pass the curve through.
    #[test]
    fn a_custom_drag_table_passes_through_unscaled() {
        let mut inputs = inputs_175gr_g7();
        inputs.custom_drag_table = Some(crate::drag::DragTable::new(
            vec![0.5, 3.0],
            vec![0.15, 0.40],
        ));
        let solver = solver(inputs);

        let sos = 340.0;
        let velocity = 0.9 * sos;
        let table_value = solver
            .inputs
            .custom_drag_table
            .as_ref()
            .expect("table")
            .interpolate(0.9);

        let reported = solver
            .effective_drag_coefficient(velocity, sos)
            .expect("mass and diameter are set");
        assert!(
            (reported - table_value).abs() < 1e-12,
            "custom table Cd {table_value} was rescaled to {reported}"
        );
    }

    /// The band step a segmented BC produces is the whole reason this field exists: it is the
    /// one feature of a real load's drag curve a consumer cannot reconstruct from a published BC.
    #[test]
    fn a_velocity_segmented_bc_steps_the_reported_cd() {
        let mut inputs = inputs_175gr_g7();
        inputs.use_bc_segments = true;
        inputs.bc_segments_data = Some(vec![
            crate::BCSegmentData { velocity_min: 2400.0, velocity_max: 4000.0, bc_value: 0.243 },
            crate::BCSegmentData { velocity_min: 0.0, velocity_max: 2400.0, bc_value: 0.200 },
        ]);
        let solver = solver(inputs);

        let sos = 340.0;
        // Straddle the 2400 fps boundary (fps -> m/s).
        let above = solver.effective_drag_coefficient(2500.0 / 3.28084, sos).expect("cd");
        let below = solver.effective_drag_coefficient(2300.0 / 3.28084, sos).expect("cd");

        // Lower BC below the boundary means MORE drag for the same reference curve.
        assert!(
            below > above,
            "expected the 0.200 band to report a higher Cd than the 0.243 band; got {below} vs {above}"
        );
    }

    /// Sectional density is undefined without both mass and diameter, and so is the projectile's
    /// own Cd. Reporting the reference value there would be a silently wrong number.
    #[test]
    fn is_absent_when_sectional_density_is_unknown() {
        let mut inputs = inputs_175gr_g7();
        inputs.weight_grains = 0.0;
        inputs.bullet_mass = 0.0;
        let solver = solver(inputs);
        assert!(solver.effective_drag_coefficient(800.0, 340.0).is_none());
    }

    /// MBA-1427: the emit rule the browser terminal applies, pinned natively because the WASM
    /// emit site is cfg-gated out of every native build and so is never executed by CI. Three
    /// cases, and the third is the contract: absent — not null — when sectional density is
    /// unknown, even with the flag set.
    #[test]
    fn the_json_emit_rule_is_flag_gated_and_absent_when_cd_is_unknown() {
        let mut point = TrajectoryPoint {
            time: 0.0,
            position: nalgebra::Vector3::new(0.0, 0.0, 0.0),
            velocity_magnitude: 800.0,
            kinetic_energy: 3000.0,
            drag_coefficient: Some(0.31),
        };
        assert_eq!(point.drag_coefficient_json_value(true), Some(0.31));
        assert_eq!(
            point.drag_coefficient_json_value(false),
            None,
            "without the flag the key must not exist, so default JSON stays byte-identical"
        );
        point.drag_coefficient = None;
        assert_eq!(
            point.drag_coefficient_json_value(true),
            None,
            "unknown sectional density must yield an ABSENT key, not null"
        );
    }

    /// Every solver family runs the same post-pass, so none may leave the field empty.
    #[test]
    fn every_point_of_a_solved_trajectory_carries_the_value() {
        let mut solver = solver(inputs_175gr_g7());
        solver.set_max_range(300.0);
        let result = solver.solve().expect("solve");
        assert!(!result.points.is_empty());
        assert!(
            result.points.iter().all(|p| p.drag_coefficient.is_some()),
            "the post-integration pass missed at least one point"
        );
    }
}
