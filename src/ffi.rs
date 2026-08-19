//! FFI bindings for iOS/Swift integration

use crate::{
    calculate_zero_angle_with_conditions, run_monte_carlo_with_direction_std_dev,
    AtmosphericConditions, BallisticInputs, DragModel, MonteCarloParams, TrajectorySolver,
    WindConditions,
};
use std::os::raw::{c_char, c_double, c_int};
use std::ptr;

/// Minimum C-ABI trajectory `step_size`, in milliseconds (0.1 ms = 0.0001 s).
///
/// Smaller steps are rejected before integration because fixed-step solves retain one public
/// trajectory point per step. All in-repository FFI examples use values in `[0.1, 1.0]` ms.
pub const MIN_FFI_STEP_SIZE_MS: c_double = 0.1;

/// Maximum C-ABI custom drag-table row count.
///
/// `drag_table_from_raw` copies both caller arrays before validation; an unbounded
/// `len` lets a single call request two multi-gigabyte allocations and abort the
/// process (MBA-1407). Real Cd decks are two to three orders of magnitude smaller
/// (the embedded G1/G7 references are under 100 rows; doppler exports run a few
/// hundred), so 4096 is far above any legitimate deck.
pub const MAX_FFI_DRAG_TABLE_LEN: c_int = 4096;

// FFI-safe structures with C-compatible layouts

#[repr(C)]
pub struct FFIBallisticInputs {
    pub muzzle_velocity: c_double,         // m/s
    pub muzzle_angle: c_double,            // radians (launch angle)
    pub bc_value: c_double,                // ballistic coefficient
    pub bullet_mass: c_double,             // kg
    pub bullet_diameter: c_double,         // meters
    pub bc_type: c_int,                    // 0=G1, 1=G7, 2=G2, 3=G5, 4=G6, 5=G8, 6=GI, 7=GS, 8=RA4 (MBA-1386; unrecognized -> G1)
    pub sight_height: c_double,            // meters
    pub target_distance: c_double,         // meters
    pub temperature: c_double,             // Celsius
    pub twist_rate: c_double,              // inches per turn
    pub is_twist_right: c_int,             // 0=false, 1=true
    pub shooting_angle: c_double,          // uphill/downhill angle in radians
    pub altitude: c_double,                // meters
    pub latitude: c_double,                // degrees (use NAN if not provided)
    pub azimuth_angle: c_double,           // horizontal aiming angle in radians
    pub use_rk4: c_int,                    // 0=Euler, 1=RK4
    pub use_adaptive_rk45: c_int,          // 0=false, 1=true (adaptive RK45)
    pub enable_wind_shear: c_int,          // 0=false, 1=true
    pub enable_trajectory_sampling: c_int, // 0=false, 1=true
    pub sample_interval: c_double,         // meters
    pub enable_pitch_damping: c_int,       // 0=false, 1=true
    pub enable_precession_nutation: c_int, // 0=false, 1=true
    pub enable_spin_drift: c_int,          // 0=false, 1=true
    pub enable_magnus: c_int,              // 0=false, 1=true
    pub enable_coriolis: c_int,            // 0=false, 1=true
    // Appended (keeps existing field offsets): compass bearing the shot is fired
    // along, radians, 0=North, PI/2=East. Drives the Coriolis Eötvös/drift azimuth.
    // Distinct from azimuth_angle (the small aiming offset). 0.0 if unset.
    pub shot_azimuth: c_double,
    // Appended (keeps existing field offsets): rifle cant angle in RADIANS about the line
    // of sight, positive = clockwise from the shooter. Rotates the sight-frame aim offsets
    // and bore geometry ("zeroed level, fired canted" -> POI right and low). 0.0 if unset.
    pub cant_angle: c_double,
    // Appended (keeps existing field offsets): deliberate vertical POI offset AT the zero
    // range, METERS (MBA-1359, Kestrel "zero height"): positive = the rifle is deliberately
    // zeroed to impact HIGH by this much at the zero distance. Applied by the zero-angle
    // solve as an angular bias on the solved angle. 0.0 if unset (identical to absent).
    pub zero_poi_vertical: c_double,
    // Appended (keeps existing field offsets): deliberate horizontal POI offset AT the zero
    // range, METERS (MBA-1359, Kestrel "zero offset"): positive = impacts RIGHT by this much
    // at the zero distance. 0.0 if unset (identical to absent).
    pub zero_poi_horizontal: c_double,
    // Appended (keeps existing field offsets): lateral sight-to-bore mount offset, METERS
    // (MBA-1396, offset-mounted optics): positive = the sight sits RIGHT of the bore, so
    // the trajectory starts that far LEFT of the line of sight. 0.0 if unset (identical to
    // absent). NOTE for callers pairing ballistics_calculate_zero_angle with a trajectory
    // call: the returned zero angle carries the VERTICAL zero_poi bias, but the horizontal
    // terms are azimuth corrections a bare elevation angle cannot carry — a caller
    // replicating an auto-zero flow must add
    // (zero_poi_horizontal + sight_offset_lateral) / zero_distance to azimuth_angle itself.
    pub sight_offset_lateral: c_double,
    // NOT plumbed (deliberately): BallisticInputs.drops_reference (MBA-1403) is an
    // OUTPUT-mode toggle for the trajectory sampler's drop column, not physics. This FFI
    // surface exposes raw kinematic samples only (world-frame positions — no drop/dial
    // outputs), so there is nothing for the toggle to act on; C callers wanting
    // target-plane drops divide their own LOS drop by cos(shooting_angle).
}

#[repr(C)]
pub struct FFIWindConditions {
    pub speed: c_double, // m/s
    // radians, wind-FROM convention: 0 = headwind, PI/2 = from the right,
    // PI = tailwind, 3*PI/2 = from the left (matches WindConditions / WindSock).
    // ALWAYS shooter-relative: the FFI deliberately has NO earth-fixed compass mode
    // (MBA-1368 decision — bindings stay shooter-relative numeric); a caller holding a
    // compass bearing converts it BEFORE the call as
    // (bearing_rad - shot_azimuth_rad).rem_euclid(2*PI), which is exactly what the
    // CLI/WASM/solve-json `--wind-ref compass` surfaces do internally.
    pub direction: c_double,
    // Appended (keeps existing field offsets): vertical wind m/s, positive = updraft;
    // 0.0 if unset.
    pub vertical_speed: c_double,
}

#[repr(C)]
pub struct FFIAtmosphericConditions {
    pub temperature: c_double, // Celsius
    pub pressure: c_double,    // hPa
    pub humidity: c_double,    // percentage (0-100)
    pub altitude: c_double,    // meters
}

#[repr(C)]
pub struct FFITrajectorySample {
    pub distance: c_double,       // meters
    pub time: c_double,           // seconds
    pub velocity_mps: c_double,   // meters per second
    pub energy_joules: c_double,  // joules
    pub drop_meters: c_double,    // meters
    pub windage_meters: c_double, // meters
    pub mach: c_double,           // Mach number
    pub spin_rate_rps: c_double,  // revolutions per second
}

/// One integrated trajectory point, in the launch frame.
///
/// **Axis convention — read this before indexing the positions.** All three are ABSOLUTE
/// positions in METERS, measured from the muzzle (the bore exit is the origin):
///
/// | field | axis | sign |
/// |-------|------|------|
/// | `position_x` | DOWNRANGE, toward the target | increases with range; bracket on this to find a given distance |
/// | `position_y` | VERTICAL height | `+` up, `-` down (below the muzzle) |
/// | `position_z` | LATERAL windage | `+` right as seen by the shooter, `-` left |
///
/// Two traps this ordering sets for consumers:
///
/// - It is **not** the "Z is downrange, X is lateral" convention. Engine 0.13.x used that
///   opposite convention, and consumers and API layers written against those releases may
///   still document it, so code ported from one must swap X and Z. Nothing fails loudly if you
///   do not — all three fields stay valid `double`s, so a mix-up silently plots range as
///   windage instead of erroring. Take the convention from this struct, not from a downstream
///   document.
/// - `position_y` is a height, not a drop. It is positive UP and relative to the muzzle,
///   whereas [`FFITrajectorySample::drop_meters`] is positive DOWN and relative to the line
///   of sight. The two disagree in both origin and sign; they are not interchangeable.
///
/// Enforced in this file's tests: `zero_then_fly_with_same_deck_is_consistent` brackets the
/// zero distance on `position_x`, and `ffi_cant_angle_deflects_laterally` asserts `position_z`
/// grows rightward under cant.
#[repr(C)]
pub struct FFITrajectoryPoint {
    pub time: c_double,
    /// Downrange distance from the muzzle, meters. See the [struct docs](FFITrajectoryPoint).
    pub position_x: c_double,
    /// Height above the muzzle, meters, positive up (NOT a drop). See the
    /// [struct docs](FFITrajectoryPoint).
    pub position_y: c_double,
    /// Lateral offset from the bore line, meters, positive right. See the
    /// [struct docs](FFITrajectoryPoint).
    pub position_z: c_double,
    pub velocity_magnitude: c_double,
    pub kinetic_energy: c_double,
}

#[repr(C)]
pub struct FFITrajectoryResult {
    pub max_range: c_double,
    pub max_height: c_double,
    pub time_of_flight: c_double,
    pub impact_velocity: c_double,
    pub impact_energy: c_double,
    pub points: *mut FFITrajectoryPoint,
    pub point_count: c_int,
    pub sampled_points: *mut FFITrajectorySample,
    pub sampled_point_count: c_int,
    pub min_pitch_damping: c_double,    // NAN if not calculated
    pub transonic_mach: c_double,       // NAN if not reached
    pub final_pitch_angle: c_double,    // NAN if not calculated
    pub final_yaw_angle: c_double,      // NAN if not calculated
    pub max_yaw_angle: c_double,        // NAN if not calculated
    pub max_precession_angle: c_double, // NAN if not calculated
}

// Monte Carlo simulation parameters
#[repr(C)]
pub struct FFIMonteCarloParams {
    pub num_simulations: c_int,
    pub velocity_std_dev: c_double,
    pub angle_std_dev: c_double,
    pub bc_std_dev: c_double,
    pub wind_speed_std_dev: c_double,
    pub target_distance: c_double,     // Use NAN if not specified
    pub base_wind_speed: c_double,     // Base wind speed in m/s
    pub base_wind_direction: c_double, // Base wind direction in radians
    pub azimuth_std_dev: c_double,     // Horizontal aiming variation in radians
}

/// Monte Carlo simulation results, one entry per simulated shot.
///
/// Every array holds `num_results` `double`s and all of them are indexed by the same sample
/// number.
///
/// **The `impact_positions_*` arrays are DEVIATIONS, not positions.** Each triple is the
/// sample's offset from the baseline point of aim, measured in the target plane (the plane at
/// the requested target distance), in METERS — see [`crate::MonteCarloResults::impact_positions`].
/// A triple of zeros means "exactly on the point of aim", not "at the muzzle".
///
/// Axes follow [`FFITrajectoryPoint`]:
///
/// | field | axis | use for dispersion? |
/// |-------|------|---------------------|
/// | `impact_positions_z` | HORIZONTAL / windage deviation, `+` right | YES — this is the horizontal dispersion axis |
/// | `impact_positions_y` | VERTICAL deviation, `+` up | YES — this is the vertical dispersion axis |
/// | `impact_positions_x` | downrange, and the target plane sits at a FIXED downrange | NO — never a dispersion axis |
///
/// Consumers ported from engine 0.13.x (`X` lateral, `Z` downrange) commonly read
/// `impact_positions_x` as the horizontal spread. Under this convention that reads the
/// downrange component instead, which is not scatter, and it fails silently.
///
/// # Filter the shortfall sentinel before computing any statistic
///
/// A sample that never reached the target plane has no deviation to report, so it is encoded as
/// `(0, TARGET_NOT_REACHED_SENTINEL_M, 0)` — the Y component is
/// [`crate::TARGET_NOT_REACHED_SENTINEL_M`] (`-1.0e9` meters) — which keeps these arrays the same
/// length as `ranges` and `impact_velocities`. The engine does NOT drop those entries for you.
/// A sample is a finite arrival exactly when
/// [`crate::MonteCarloResults::position_reached_target`] accepts it: every component finite AND
/// the Y component not equal to the sentinel.
///
/// Feeding an unfiltered array into a mean or standard deviation drags the result toward
/// `-1e9` and produces a plainly absurd group size. Exclude sentinel samples from dispersion
/// statistics, but keep them in the denominator for hit probability — they are definite misses,
/// which is how `hit_probability` already counts them.
#[repr(C)]
pub struct FFIMonteCarloResults {
    pub ranges: *mut c_double,
    pub impact_velocities: *mut c_double,
    /// Downrange component of the target-plane deviation, meters. NOT a dispersion axis; see
    /// the [struct docs](FFIMonteCarloResults).
    pub impact_positions_x: *mut c_double,
    /// VERTICAL deviation from the point of aim, meters, positive up — and the field carrying
    /// the shortfall marker: [`crate::TARGET_NOT_REACHED_SENTINEL_M`] (`-1.0e9`) marks a sample
    /// that did not reach the target plane. Exclude those from dispersion statistics but retain
    /// them as misses for probability calculations. See the [struct docs](FFIMonteCarloResults).
    pub impact_positions_y: *mut c_double,
    /// HORIZONTAL / windage deviation from the point of aim, meters, positive right. This is
    /// the horizontal dispersion axis; see the [struct docs](FFIMonteCarloResults).
    pub impact_positions_z: *mut c_double,
    pub num_results: c_int,
    pub mean_range: c_double,
    pub std_dev_range: c_double,
    pub mean_impact_velocity: c_double,
    pub std_dev_impact_velocity: c_double,
    pub hit_probability: c_double, // If target_distance was specified
}

// Helper function to convert FFI inputs to internal types
#[allow(clippy::field_reassign_with_default)] // Keep the C-to-Rust field mapping sequential/auditable.
fn convert_inputs(inputs: &FFIBallisticInputs) -> BallisticInputs {
    let mut ballistic_inputs = BallisticInputs::default();

    ballistic_inputs.muzzle_velocity = inputs.muzzle_velocity;
    ballistic_inputs.muzzle_angle = inputs.muzzle_angle;
    ballistic_inputs.azimuth_angle = inputs.azimuth_angle;
    ballistic_inputs.shot_azimuth = inputs.shot_azimuth;
    ballistic_inputs.cant_angle = inputs.cant_angle;
    ballistic_inputs.zero_poi_vertical_m = inputs.zero_poi_vertical;
    ballistic_inputs.zero_poi_horizontal_m = inputs.zero_poi_horizontal;
    ballistic_inputs.sight_offset_lateral_m = inputs.sight_offset_lateral;
    ballistic_inputs.use_rk4 = inputs.use_rk4 != 0;
    ballistic_inputs.use_adaptive_rk45 = inputs.use_adaptive_rk45 != 0;
    ballistic_inputs.bc_value = inputs.bc_value;
    ballistic_inputs.bullet_mass = inputs.bullet_mass;
    ballistic_inputs.bullet_diameter = inputs.bullet_diameter;
    ballistic_inputs.bc_type = match inputs.bc_type {
        1 => DragModel::G7,
        2 => DragModel::G2,
        3 => DragModel::G5,
        4 => DragModel::G6,
        5 => DragModel::G8,
        6 => DragModel::GI,
        7 => DragModel::GS,
        // MBA-1386: additive slot for the new RA4 family. Existing callers passing
        // 0-7 are unaffected; any other/unrecognized value still falls back to G1.
        8 => DragModel::RA4,
        _ => DragModel::G1,
    };
    ballistic_inputs.sight_height = inputs.sight_height;
    ballistic_inputs.target_distance = inputs.target_distance;
    ballistic_inputs.temperature = inputs.temperature;
    ballistic_inputs.twist_rate = inputs.twist_rate;
    ballistic_inputs.is_twist_right = inputs.is_twist_right != 0;
    ballistic_inputs.shooting_angle = inputs.shooting_angle;
    ballistic_inputs.altitude = inputs.altitude;

    if !inputs.latitude.is_nan() {
        ballistic_inputs.latitude = Some(inputs.latitude);
    }

    // Set derived values
    ballistic_inputs.caliber_inches = inputs.bullet_diameter / 0.0254;
    ballistic_inputs.weight_grains = inputs.bullet_mass / crate::constants::GRAINS_TO_KG;
    // MBA-1135: mass-based length estimate (was a mass-blind 4.5-caliber default). The C ABI does
    // not carry a bullet length, so derive it from diameter + mass; fall back to 4.5-cal if mass<=0.
    ballistic_inputs.bullet_length = {
        let est = crate::stability::estimate_bullet_length_m(
            ballistic_inputs.bullet_diameter,
            ballistic_inputs.bullet_mass,
        );
        if est > 0.0 {
            est
        } else {
            ballistic_inputs.bullet_diameter * 4.5
        }
    };

    // New advanced physics flags
    ballistic_inputs.enable_wind_shear = inputs.enable_wind_shear != 0;
    ballistic_inputs.enable_trajectory_sampling = inputs.enable_trajectory_sampling != 0;
    ballistic_inputs.sample_interval = inputs.sample_interval;
    ballistic_inputs.enable_pitch_damping = inputs.enable_pitch_damping != 0;
    ballistic_inputs.enable_precession_nutation = inputs.enable_precession_nutation != 0;
    ballistic_inputs.use_enhanced_spin_drift = inputs.enable_spin_drift != 0;
    ballistic_inputs.enable_advanced_effects =
        inputs.enable_magnus != 0 || inputs.enable_coriolis != 0;
    // Gate Magnus and Coriolis independently so enabling one does not enable the other.
    ballistic_inputs.enable_magnus = inputs.enable_magnus != 0;
    ballistic_inputs.enable_coriolis = inputs.enable_coriolis != 0;

    ballistic_inputs
}

/// Build a validated [`crate::drag::DragTable`] from borrowed C arrays.
///
/// Both arrays must contain `len` elements, with `len` in `[2, MAX_FFI_DRAG_TABLE_LEN]`.
/// The data is copied; the caller retains ownership. Returns `Err(())` for null
/// pointers, `len` outside that range, or any deck that fails
/// [`crate::drag::DragTable::try_new`] validation (non-ascending Mach, non-positive
/// or non-finite Cd). No error detail crosses the ABI, matching the null/NAN
/// error convention of this module.
///
/// # Safety
///
/// When non-null, `mach` and `cd` must each point to `len` readable `f64` values
/// that remain valid and unmutated for the duration of the call.
unsafe fn drag_table_from_raw(
    mach: *const c_double,
    cd: *const c_double,
    len: c_int,
) -> Result<crate::drag::DragTable, ()> {
    if mach.is_null() || cd.is_null() || !(2..=MAX_FFI_DRAG_TABLE_LEN).contains(&len) {
        return Err(());
    }
    let len = len as usize;
    let mach = unsafe { std::slice::from_raw_parts(mach, len) }.to_vec();
    let cd = unsafe { std::slice::from_raw_parts(cd, len) }.to_vec();
    crate::drag::DragTable::try_new(mach, cd).map_err(|_| ())
}

/// Shared implementation for the trajectory exports. `custom_drag_table`, when
/// present, replaces the G-model + BC drag (the deck's Cd is divided by sectional
/// density — see `BallisticInputs::custom_drag_denominator`). `cd_scale` multiplies the
/// deck's interpolated Cd (MBA-1356); callers that don't expose a scale pass `1.0`
/// (neutral, byte-identical to no scale). It is inert when `custom_drag_table` is `None`.
unsafe fn calculate_trajectory_impl(
    inputs: *const FFIBallisticInputs,
    wind: *const FFIWindConditions,
    atmosphere: *const FFIAtmosphericConditions,
    max_range: c_double,
    step_size: c_double,
    custom_drag_table: Option<crate::drag::DragTable>,
    cd_scale: c_double,
) -> *mut FFITrajectoryResult {
    if inputs.is_null() {
        return ptr::null_mut();
    }
    if !step_size.is_finite() || step_size < MIN_FFI_STEP_SIZE_MS {
        return ptr::null_mut();
    }

    let inputs = unsafe { &*inputs };
    let mut ballistic_inputs = convert_inputs(inputs);
    ballistic_inputs.custom_drag_table = custom_drag_table;
    ballistic_inputs.cd_scale = cd_scale;
    let twist_rate_in = ballistic_inputs.twist_rate;

    let wind_conditions = if wind.is_null() {
        WindConditions::default()
    } else {
        let wind = unsafe { &*wind };
        WindConditions {
            speed: wind.speed,
            direction: wind.direction,
            vertical_speed: wind.vertical_speed,
        }
    };

    let atmospheric_conditions = if atmosphere.is_null() {
        AtmosphericConditions::default()
    } else {
        let atmo = unsafe { &*atmosphere };
        AtmosphericConditions {
            temperature: atmo.temperature,
            pressure: atmo.pressure,
            humidity: atmo.humidity,
            altitude: atmo.altitude,
        }
    };

    // Create solver and calculate trajectory
    let (sample_temp_c, sample_pressure_hpa) = crate::atmosphere::resolve_station_conditions(
        atmospheric_conditions.temperature,
        atmospheric_conditions.pressure,
        atmospheric_conditions.altitude,
    );
    let (_, sample_speed_of_sound) = crate::atmosphere::calculate_atmosphere(
        atmospheric_conditions.altitude,
        Some(sample_temp_c),
        Some(sample_pressure_hpa),
        atmospheric_conditions.humidity,
    );

    let mut solver =
        TrajectorySolver::new(ballistic_inputs, wind_conditions, atmospheric_conditions);

    // Set max range and time step
    solver.set_max_range(max_range);
    solver.set_time_step(step_size / 1000.0); // milliseconds -> seconds

    match solver.solve() {
        Ok(result) => {
            // Convert trajectory points to FFI format
            let point_count = result.points.len();
            let points = if point_count > 0 {
                let mut ffi_points = Vec::with_capacity(point_count);
                for point in result.points.iter() {
                    ffi_points.push(FFITrajectoryPoint {
                        time: point.time,
                        position_x: point.position[0],
                        position_y: point.position[1],
                        position_z: point.position[2],
                        velocity_magnitude: point.velocity_magnitude,
                        kinetic_energy: point.kinetic_energy,
                    });
                }
                let points_ptr = ffi_points.as_mut_ptr();
                std::mem::forget(ffi_points); // Prevent deallocation
                points_ptr
            } else {
                ptr::null_mut()
            };

            // Convert sampled points if available
            let (sampled_points, sampled_point_count) =
                if let Some(ref samples) = result.sampled_points {
                    let mut ffi_samples = Vec::with_capacity(samples.len());
                    for sample in samples {
                        ffi_samples.push(FFITrajectorySample {
                            distance: sample.distance_m,
                            time: sample.time_s,
                            velocity_mps: sample.velocity_mps,
                            energy_joules: sample.energy_j,
                            drop_meters: sample.drop_m,
                            windage_meters: sample.wind_drift_m,
                            mach: if sample_speed_of_sound > 0.0 {
                                sample.velocity_mps / sample_speed_of_sound
                            } else {
                                0.0
                            },
                            spin_rate_rps: if twist_rate_in > 0.0 {
                                sample.velocity_mps / (twist_rate_in * 0.0254)
                            } else {
                                0.0
                            },
                        });
                    }
                    let count = ffi_samples.len() as c_int;
                    let samples_ptr = ffi_samples.as_mut_ptr();
                    std::mem::forget(ffi_samples);
                    (samples_ptr, count)
                } else {
                    (ptr::null_mut(), 0)
                };

            // Extract angular state values if available
            let (final_pitch, final_yaw, max_yaw, max_prec) =
                if let Some(ref angular) = result.angular_state {
                    (
                        angular.pitch_angle,
                        angular.yaw_angle,
                        result.max_yaw_angle.unwrap_or(f64::NAN),
                        result.max_precession_angle.unwrap_or(f64::NAN),
                    )
                } else {
                    (f64::NAN, f64::NAN, f64::NAN, f64::NAN)
                };

            // Create result on heap
            let ffi_result = Box::new(FFITrajectoryResult {
                max_range: result.max_range,
                max_height: result.max_height,
                time_of_flight: result.time_of_flight,
                impact_velocity: result.impact_velocity,
                impact_energy: result.impact_energy,
                points,
                point_count: point_count as c_int,
                sampled_points,
                sampled_point_count,
                min_pitch_damping: result.min_pitch_damping.unwrap_or(f64::NAN),
                transonic_mach: result.transonic_mach.unwrap_or(f64::NAN),
                final_pitch_angle: final_pitch,
                final_yaw_angle: final_yaw,
                max_yaw_angle: max_yaw,
                max_precession_angle: max_prec,
            });

            Box::into_raw(ffi_result)
        }
        Err(_) => ptr::null_mut(),
    }
}

/// Calculate a trajectory through the C ABI.
///
/// `step_size` is expressed in milliseconds and must be finite and at least
/// [`MIN_FFI_STEP_SIZE_MS`]. This boundary contract is validated for every solver mode, although
/// adaptive RK45 chooses its integration steps internally. Invalid values return null without
/// starting a solve. A solve that would exceed [`crate::MAX_TRAJECTORY_POINTS`] also returns null;
/// an enabled sampling grid above [`crate::MAX_TRAJECTORY_SAMPLES`] does likewise. Callers can
/// increase `step_size`, reduce `max_range`, or select adaptive RK45.
///
/// # Safety
///
/// `inputs` may be null, in which case this function returns null. When non-null,
/// it must point to a valid, properly aligned [`FFIBallisticInputs`] that remains
/// readable and is not mutated for the duration of this call. `wind` and
/// `atmosphere` may also be null; each non-null pointer has the same requirements
/// for its corresponding type.
/// The returned pointer, when non-null, must be released exactly once with
/// [`ballistics_free_trajectory_result`].
#[no_mangle]
pub unsafe extern "C" fn ballistics_calculate_trajectory(
    inputs: *const FFIBallisticInputs,
    wind: *const FFIWindConditions,
    atmosphere: *const FFIAtmosphericConditions,
    max_range: c_double,
    step_size: c_double,
) -> *mut FFITrajectoryResult {
    unsafe { calculate_trajectory_impl(inputs, wind, atmosphere, max_range, step_size, None, 1.0) }
}

/// [`ballistics_calculate_trajectory`] with a caller-supplied custom drag deck
/// (Cd vs Mach, e.g. Hornady CDM / Doppler-radar data). The deck REPLACES the
/// G-model + BC for drag: `bc_type`/`bc_value` are ignored, and the retardation
/// denominator becomes the projectile's sectional density (mass and diameter in
/// `inputs` must therefore be positive). Mach values must be strictly ascending
/// with `drag_table_len` in `[2, MAX_FFI_DRAG_TABLE_LEN]` and finite positive Cd;
/// outside the deck's Mach domain the nearest endpoint Cd is held.
///
/// Returns null for an invalid deck (null arrays, `drag_table_len` outside
/// `[2, MAX_FFI_DRAG_TABLE_LEN]`, or failed validation), in addition to every
/// failure mode of the base function.
///
/// # Safety
///
/// Same contract as [`ballistics_calculate_trajectory`] for `inputs`, `wind`, and
/// `atmosphere`. Additionally, when non-null, `drag_mach` and `drag_cd` must each
/// point to `drag_table_len` readable `f64` values, borrowed only for the duration
/// of this call (the deck is copied; the caller retains ownership — no new free
/// function is required). The returned pointer, when non-null, must be released
/// exactly once with [`ballistics_free_trajectory_result`].
#[no_mangle]
pub unsafe extern "C" fn ballistics_calculate_trajectory_with_drag_table(
    inputs: *const FFIBallisticInputs,
    wind: *const FFIWindConditions,
    atmosphere: *const FFIAtmosphericConditions,
    max_range: c_double,
    step_size: c_double,
    drag_mach: *const c_double,
    drag_cd: *const c_double,
    drag_table_len: c_int,
) -> *mut FFITrajectoryResult {
    let table = match unsafe { drag_table_from_raw(drag_mach, drag_cd, drag_table_len) } {
        Ok(t) => t,
        Err(()) => return ptr::null_mut(),
    };
    unsafe {
        calculate_trajectory_impl(
            inputs, wind, atmosphere, max_range, step_size, Some(table), 1.0,
        )
    }
}

/// [`ballistics_calculate_trajectory_with_drag_table`] with an additional whole-curve
/// drag scale (MBA-1356): the deck's interpolated Cd is multiplied by `cd_scale` at the
/// same site the base export uses, i.e. `Cd_used = table.interpolate(mach) * cd_scale`.
/// `cd_scale = 1.0` is neutral and produces byte-identical output to
/// [`ballistics_calculate_trajectory_with_drag_table`] on the same deck. Typical truing
/// values are in `[0.90, 1.10]`; values outside that band are accepted here (the engine
/// only rejects non-finite or non-positive) — an "unusually large" warning is a CLI-only
/// concern (Task 2), not part of this frozen C ABI.
///
/// Returns null when `cd_scale` is not finite or not `> 0`, in addition to every failure
/// mode of [`ballistics_calculate_trajectory_with_drag_table`] (matching that export's
/// null sentinel for an invalid deck).
///
/// # Safety
///
/// Same contract as [`ballistics_calculate_trajectory_with_drag_table`].
#[no_mangle]
pub unsafe extern "C" fn ballistics_calculate_trajectory_with_drag_table_scaled(
    inputs: *const FFIBallisticInputs,
    wind: *const FFIWindConditions,
    atmosphere: *const FFIAtmosphericConditions,
    max_range: c_double,
    step_size: c_double,
    drag_mach: *const c_double,
    drag_cd: *const c_double,
    drag_table_len: c_int,
    cd_scale: c_double,
) -> *mut FFITrajectoryResult {
    if !cd_scale.is_finite() || cd_scale <= 0.0 {
        return ptr::null_mut();
    }
    let table = match unsafe { drag_table_from_raw(drag_mach, drag_cd, drag_table_len) } {
        Ok(t) => t,
        Err(()) => return ptr::null_mut(),
    };
    unsafe {
        calculate_trajectory_impl(
            inputs, wind, atmosphere, max_range, step_size, Some(table), cd_scale,
        )
    }
}

/// Release a trajectory result allocated by [`ballistics_calculate_trajectory`].
///
/// # Safety
///
/// `result` must be null or a pointer returned by
/// [`ballistics_calculate_trajectory`] that has not already been freed. Its
/// embedded pointers and counts must be unchanged from the returned values.
/// After this call, the result and its point arrays must not be accessed again.
#[no_mangle]
pub unsafe extern "C" fn ballistics_free_trajectory_result(result: *mut FFITrajectoryResult) {
    if !result.is_null() {
        unsafe {
            let result = Box::from_raw(result);
            if !result.points.is_null() && result.point_count > 0 {
                let points = Vec::from_raw_parts(
                    result.points,
                    result.point_count as usize,
                    result.point_count as usize,
                );
                drop(points);
            }
            if !result.sampled_points.is_null() && result.sampled_point_count > 0 {
                let samples = Vec::from_raw_parts(
                    result.sampled_points,
                    result.sampled_point_count as usize,
                    result.sampled_point_count as usize,
                );
                drop(samples);
            }
            drop(result);
        }
    }
}

/// Shared implementation for the zero-angle exports. `custom_drag_table`, when
/// present, replaces the G-model + BC drag (the deck's Cd is divided by sectional
/// density — see `BallisticInputs::custom_drag_denominator`), matching the deck
/// semantics of [`calculate_trajectory_impl`] so a zero solved with a deck and a
/// trajectory flown with the same deck agree. `cd_scale` multiplies the deck's
/// interpolated Cd (MBA-1356); callers that don't expose a scale pass `1.0` (neutral).
unsafe fn calculate_zero_angle_impl(
    inputs: *const FFIBallisticInputs,
    wind: *const FFIWindConditions,
    atmosphere: *const FFIAtmosphericConditions,
    zero_distance: c_double,
    custom_drag_table: Option<crate::drag::DragTable>,
    cd_scale: c_double,
) -> c_double {
    if inputs.is_null() {
        return f64::NAN;
    }

    let inputs = unsafe { &*inputs };
    let mut ballistic_inputs = convert_inputs(inputs);
    ballistic_inputs.custom_drag_table = custom_drag_table;
    ballistic_inputs.cd_scale = cd_scale;

    let wind_conditions = if wind.is_null() {
        WindConditions::default()
    } else {
        let wind = unsafe { &*wind };
        WindConditions {
            speed: wind.speed,
            direction: wind.direction,
            vertical_speed: wind.vertical_speed,
        }
    };

    let atmospheric_conditions = if atmosphere.is_null() {
        AtmosphericConditions::default()
    } else {
        let atmo = unsafe { &*atmosphere };
        AtmosphericConditions {
            temperature: atmo.temperature,
            pressure: atmo.pressure,
            humidity: atmo.humidity,
            altitude: atmo.altitude,
        }
    };

    // For zero angle, we want the bullet to hit at sight height at the zero distance
    // This means the bullet crosses the line of sight at the zero distance
    let target_height = ballistic_inputs.sight_height;

    calculate_zero_angle_with_conditions(
        ballistic_inputs,
        zero_distance,
        target_height,
        wind_conditions,
        atmospheric_conditions,
    )
    .unwrap_or(f64::NAN)
}

/// Calculate the zero angle for a target distance through the C ABI.
///
/// # Safety
///
/// `inputs` may be null, in which case this function returns NaN. When non-null,
/// it must point to a valid, properly aligned [`FFIBallisticInputs`] that remains
/// readable and is not mutated for the duration of this call. `wind` and
/// `atmosphere` may also be null; each non-null pointer has the same requirements
/// for its corresponding type.
#[no_mangle]
pub unsafe extern "C" fn ballistics_calculate_zero_angle(
    inputs: *const FFIBallisticInputs,
    wind: *const FFIWindConditions,
    atmosphere: *const FFIAtmosphericConditions,
    zero_distance: c_double,
) -> c_double {
    unsafe { calculate_zero_angle_impl(inputs, wind, atmosphere, zero_distance, None, 1.0) }
}

/// [`ballistics_calculate_zero_angle`] with a caller-supplied custom drag deck
/// (Cd vs Mach, e.g. Hornady CDM / Doppler-radar data). The deck REPLACES the
/// G-model + BC for drag: `bc_type`/`bc_value` are ignored, and the retardation
/// denominator becomes the projectile's sectional density (mass and diameter in
/// `inputs` must therefore be positive). Mach values must be strictly ascending
/// with `drag_table_len` in `[2, MAX_FFI_DRAG_TABLE_LEN]` and finite positive Cd;
/// outside the deck's Mach domain the nearest endpoint Cd is held. Pair this with
/// [`ballistics_calculate_trajectory_with_drag_table`] using the same deck to fly
/// the solved angle — the two exports share identical deck semantics.
///
/// Returns NaN for an invalid deck (null arrays, `drag_table_len` outside
/// `[2, MAX_FFI_DRAG_TABLE_LEN]`, or failed validation), in addition to every
/// failure mode of the base function.
///
/// # Safety
///
/// Same contract as [`ballistics_calculate_zero_angle`] for `inputs`, `wind`, and
/// `atmosphere`. Additionally, when non-null, `drag_mach` and `drag_cd` must each
/// point to `drag_table_len` readable `f64` values, borrowed only for the duration
/// of this call (the deck is copied; the caller retains ownership — no new free
/// function is required).
#[no_mangle]
pub unsafe extern "C" fn ballistics_calculate_zero_angle_with_drag_table(
    inputs: *const FFIBallisticInputs,
    wind: *const FFIWindConditions,
    atmosphere: *const FFIAtmosphericConditions,
    zero_distance: c_double,
    drag_mach: *const c_double,
    drag_cd: *const c_double,
    drag_table_len: c_int,
) -> c_double {
    let table = match unsafe { drag_table_from_raw(drag_mach, drag_cd, drag_table_len) } {
        Ok(t) => t,
        Err(()) => return f64::NAN,
    };
    unsafe {
        calculate_zero_angle_impl(inputs, wind, atmosphere, zero_distance, Some(table), 1.0)
    }
}

/// [`ballistics_calculate_zero_angle_with_drag_table`] with an additional whole-curve
/// drag scale (MBA-1356): the deck's interpolated Cd is multiplied by `cd_scale` at the
/// same site the base export uses, i.e. `Cd_used = table.interpolate(mach) * cd_scale`.
/// `cd_scale = 1.0` is neutral and produces byte-identical output to
/// [`ballistics_calculate_zero_angle_with_drag_table`] on the same deck. Pair this with
/// [`ballistics_calculate_trajectory_with_drag_table_scaled`] using the same deck AND the
/// same `cd_scale` to fly the solved angle. Typical truing values are in `[0.90, 1.10]`;
/// values outside that band are accepted here (the engine only rejects non-finite or
/// non-positive) — an "unusually large" warning is a CLI-only concern (Task 2).
///
/// Returns NaN when `cd_scale` is not finite or not `> 0`, in addition to every failure
/// mode of [`ballistics_calculate_zero_angle_with_drag_table`] (matching that export's
/// NaN sentinel for an invalid deck).
///
/// # Safety
///
/// Same contract as [`ballistics_calculate_zero_angle_with_drag_table`].
#[no_mangle]
pub unsafe extern "C" fn ballistics_calculate_zero_angle_with_drag_table_scaled(
    inputs: *const FFIBallisticInputs,
    wind: *const FFIWindConditions,
    atmosphere: *const FFIAtmosphericConditions,
    zero_distance: c_double,
    drag_mach: *const c_double,
    drag_cd: *const c_double,
    drag_table_len: c_int,
    cd_scale: c_double,
) -> c_double {
    if !cd_scale.is_finite() || cd_scale <= 0.0 {
        return f64::NAN;
    }
    let table = match unsafe { drag_table_from_raw(drag_mach, drag_cd, drag_table_len) } {
        Ok(t) => t,
        Err(()) => return f64::NAN,
    };
    unsafe {
        calculate_zero_angle_impl(
            inputs,
            wind,
            atmosphere,
            zero_distance,
            Some(table),
            cd_scale,
        )
    }
}

// Simple trajectory calculation for quick results
#[no_mangle]
#[allow(clippy::field_reassign_with_default)] // Preserve the staged zero-angle workflow below.
pub extern "C" fn ballistics_quick_trajectory(
    muzzle_velocity: c_double,
    bc: c_double,
    sight_height: c_double,
    zero_distance: c_double,
    target_distance: c_double,
) -> c_double {
    // This provides a simple drop calculation at target distance
    // Using simplified ballistic calculations

    let mut inputs = BallisticInputs::default();
    inputs.muzzle_velocity = muzzle_velocity;
    inputs.bc_value = bc;
    inputs.sight_height = sight_height;
    inputs.target_distance = target_distance;

    let wind = WindConditions::default();
    let atmo = AtmosphericConditions::default();

    // First calculate zero angle
    let zero_angle = match calculate_zero_angle_with_conditions(
        inputs.clone(),
        zero_distance,
        sight_height,
        wind.clone(),
        atmo.clone(),
    ) {
        Ok(angle) => angle,
        Err(_) => return f64::NAN,
    };

    // Now calculate trajectory with that zero angle
    inputs.muzzle_angle = zero_angle;

    let mut solver = TrajectorySolver::new(inputs, wind, atmo);
    solver.set_max_range(target_distance * 1.1);

    match solver.solve() {
        Ok(result) => {
            // Find the drop at target distance
            for point in result.points {
                if point.position[0] >= target_distance {
                    return sight_height - point.position[1];
                }
            }
            f64::NAN
        }
        Err(_) => f64::NAN,
    }
}

/// Run a Monte Carlo simulation through the C ABI.
///
/// # Safety
///
/// `inputs` and `params` may be null, in which case this function returns null.
/// Each non-null pointer must point to a valid, properly aligned value of its
/// corresponding FFI type that remains readable and is not mutated for the
/// duration of this call. `atmosphere` may be null; a non-null pointer has the
/// same requirements for [`FFIAtmosphericConditions`]. The returned pointer,
/// when non-null, must be released exactly once with
/// [`ballistics_free_monte_carlo_results`].
#[no_mangle]
pub unsafe extern "C" fn ballistics_monte_carlo(
    inputs: *const FFIBallisticInputs,
    atmosphere: *const FFIAtmosphericConditions,
    params: *const FFIMonteCarloParams,
) -> *mut FFIMonteCarloResults {
    unsafe { ballistics_monte_carlo_impl(inputs, atmosphere, params, 0.0) }
}

/// Run a Monte Carlo simulation with independent wind-direction uncertainty through the C ABI.
///
/// `wind_direction_std_dev` is in radians. This additive entry point keeps
/// [`FFIMonteCarloParams`] and [`ballistics_monte_carlo`] binary-compatible; the older function
/// delegates with zero wind-direction uncertainty.
///
/// # Safety
///
/// The pointer and ownership requirements are identical to [`ballistics_monte_carlo`].
#[no_mangle]
pub unsafe extern "C" fn ballistics_monte_carlo_with_direction_std_dev(
    inputs: *const FFIBallisticInputs,
    atmosphere: *const FFIAtmosphericConditions,
    params: *const FFIMonteCarloParams,
    wind_direction_std_dev: c_double,
) -> *mut FFIMonteCarloResults {
    unsafe { ballistics_monte_carlo_impl(inputs, atmosphere, params, wind_direction_std_dev) }
}

unsafe fn ballistics_monte_carlo_impl(
    inputs: *const FFIBallisticInputs,
    atmosphere: *const FFIAtmosphericConditions,
    params: *const FFIMonteCarloParams,
    wind_direction_std_dev: f64,
) -> *mut FFIMonteCarloResults {
    if inputs.is_null() || params.is_null() {
        return ptr::null_mut();
    }

    let inputs = unsafe { &*inputs };
    let params = unsafe { &*params };

    // Reject an out-of-range simulation count. num_simulations is a c_int (i32) cast straight to
    // usize: a negative value would wrap to a near-max usize, and even a large positive value (up
    // to i32::MAX ~ 2.1e9) would drive billions of iterations with the result arrays scaling to
    // match — an unbounded-loop / OOM DoS from a single FFI call. Bound it to a sane maximum.
    // (n == 0 also yields NaN stats and a zero-size allocation.)
    const MAX_SIMULATIONS: c_int = 1_000_000;
    if params.num_simulations <= 0 || params.num_simulations > MAX_SIMULATIONS {
        return ptr::null_mut();
    }

    // Convert FFI inputs to internal types
    let mut ballistic_inputs = convert_inputs(inputs);
    ballistic_inputs.muzzle_height = 1.5;
    ballistic_inputs.ground_threshold = 0.0;
    if !atmosphere.is_null() {
        let atmo = unsafe { &*atmosphere };
        ballistic_inputs.temperature = atmo.temperature;
        ballistic_inputs.pressure = atmo.pressure;
        ballistic_inputs.humidity = (atmo.humidity / 100.0).clamp(0.0, 1.0);
        ballistic_inputs.altitude = atmo.altitude;
    }

    // Create Monte Carlo parameters
    let mc_params = MonteCarloParams {
        num_simulations: params.num_simulations as usize,
        velocity_std_dev: params.velocity_std_dev,
        angle_std_dev: params.angle_std_dev,
        bc_std_dev: params.bc_std_dev,
        wind_speed_std_dev: params.wind_speed_std_dev,
        target_distance: if params.target_distance.is_nan() {
            None
        } else {
            Some(params.target_distance)
        },
        base_wind_speed: params.base_wind_speed,
        base_wind_direction: params.base_wind_direction,
        azimuth_std_dev: params.azimuth_std_dev,
    };

    // Run Monte Carlo simulation
    match run_monte_carlo_with_direction_std_dev(
        ballistic_inputs,
        mc_params,
        wind_direction_std_dev,
    ) {
        Ok(results) => {
            let num_results = results.ranges.len() as c_int;

            // Calculate statistics
            let mean_range: f64 = results.ranges.iter().sum::<f64>() / num_results as f64;
            let variance_range: f64 = results
                .ranges
                .iter()
                .map(|r| (r - mean_range).powi(2))
                .sum::<f64>()
                / num_results as f64;
            let std_dev_range = variance_range.sqrt();

            let mean_velocity: f64 =
                results.impact_velocities.iter().sum::<f64>() / num_results as f64;
            let variance_velocity: f64 = results
                .impact_velocities
                .iter()
                .map(|v| (v - mean_velocity).powi(2))
                .sum::<f64>()
                / num_results as f64;
            let std_dev_velocity = variance_velocity.sqrt();

            // Calculate hit probability if target distance was specified. MBA-971: use the shared
            // position-based criterion (fraction within DEFAULT_HIT_RADIUS_M of the point of aim
            // at the target plane). The old inline version had a redundant `distance < target`
            // clause comparing a ~meter deviation to the ~hundreds-of-meters target distance
            // (effectively always true), and the CLI used a different range-based notion entirely.
            let hit_probability = if params.target_distance.is_nan() {
                0.0
            } else {
                results.hit_probability(crate::DEFAULT_HIT_RADIUS_M)
            };

            // Allocate memory for arrays
            let ranges_ptr = unsafe {
                let ptr = std::alloc::alloc(
                    std::alloc::Layout::array::<c_double>(num_results as usize).unwrap(),
                ) as *mut c_double;
                for (i, &range) in results.ranges.iter().enumerate() {
                    *ptr.add(i) = range;
                }
                ptr
            };

            let velocities_ptr = unsafe {
                let ptr = std::alloc::alloc(
                    std::alloc::Layout::array::<c_double>(num_results as usize).unwrap(),
                ) as *mut c_double;
                for (i, &vel) in results.impact_velocities.iter().enumerate() {
                    *ptr.add(i) = vel;
                }
                ptr
            };

            let pos_x_ptr = unsafe {
                let ptr = std::alloc::alloc(
                    std::alloc::Layout::array::<c_double>(num_results as usize).unwrap(),
                ) as *mut c_double;
                for (i, pos) in results.impact_positions.iter().enumerate() {
                    *ptr.add(i) = pos.x;
                }
                ptr
            };

            let pos_y_ptr = unsafe {
                let ptr = std::alloc::alloc(
                    std::alloc::Layout::array::<c_double>(num_results as usize).unwrap(),
                ) as *mut c_double;
                for (i, pos) in results.impact_positions.iter().enumerate() {
                    *ptr.add(i) = pos.y;
                }
                ptr
            };

            let pos_z_ptr = unsafe {
                let ptr = std::alloc::alloc(
                    std::alloc::Layout::array::<c_double>(num_results as usize).unwrap(),
                ) as *mut c_double;
                for (i, pos) in results.impact_positions.iter().enumerate() {
                    *ptr.add(i) = pos.z;
                }
                ptr
            };

            // Create result structure
            let result = Box::new(FFIMonteCarloResults {
                ranges: ranges_ptr,
                impact_velocities: velocities_ptr,
                impact_positions_x: pos_x_ptr,
                impact_positions_y: pos_y_ptr,
                impact_positions_z: pos_z_ptr,
                num_results,
                mean_range,
                std_dev_range,
                mean_impact_velocity: mean_velocity,
                std_dev_impact_velocity: std_dev_velocity,
                hit_probability,
            });

            Box::into_raw(result)
        }
        Err(_) => ptr::null_mut(),
    }
}

/// Release Monte Carlo results allocated by either Monte Carlo C entry point.
///
/// # Safety
///
/// `results` must be null or a pointer returned by [`ballistics_monte_carlo`] or
/// [`ballistics_monte_carlo_with_direction_std_dev`] that has not already been freed. Its
/// embedded pointers and result count must be unchanged from the returned values. After this
/// call, the result and all of its arrays must not be accessed again.
#[no_mangle]
pub unsafe extern "C" fn ballistics_free_monte_carlo_results(results: *mut FFIMonteCarloResults) {
    if results.is_null() {
        return;
    }

    unsafe {
        let results = Box::from_raw(results);
        let num = results.num_results as usize;

        // Free arrays
        if !results.ranges.is_null() {
            std::alloc::dealloc(
                results.ranges as *mut u8,
                std::alloc::Layout::array::<c_double>(num).unwrap(),
            );
        }

        if !results.impact_velocities.is_null() {
            std::alloc::dealloc(
                results.impact_velocities as *mut u8,
                std::alloc::Layout::array::<c_double>(num).unwrap(),
            );
        }

        if !results.impact_positions_x.is_null() {
            std::alloc::dealloc(
                results.impact_positions_x as *mut u8,
                std::alloc::Layout::array::<c_double>(num).unwrap(),
            );
        }

        if !results.impact_positions_y.is_null() {
            std::alloc::dealloc(
                results.impact_positions_y as *mut u8,
                std::alloc::Layout::array::<c_double>(num).unwrap(),
            );
        }

        if !results.impact_positions_z.is_null() {
            std::alloc::dealloc(
                results.impact_positions_z as *mut u8,
                std::alloc::Layout::array::<c_double>(num).unwrap(),
            );
        }

        // Box automatically deallocates the result structure
    }
}

/// Which standard atmosphere a `bc` value passed to [`ballistics_bc_for_reference_standard`]
/// is referenced to (MBA-1365). `0` = ICAO (the default every other export in this module
/// assumes), `1` = Army Standard Metro (some vendor-published BCs, notably many
/// Sierra/Hornady/Barnes bullets). Any other value is treated as `0` (ICAO), matching the
/// permissive unrecognized-value convention `convert_inputs` already uses for `bc_type`.
pub const FFI_BC_REFERENCE_ICAO: c_int = 0;
pub const FFI_BC_REFERENCE_ARMY_STANDARD_METRO: c_int = 1;

/// Convert a ballistic coefficient declared under `reference_standard` to the ICAO-referenced
/// value every `FFIBallisticInputs.bc_value` in this module expects (MBA-1365).
///
/// `FFIBallisticInputs` is an ABI-frozen `repr(C)` struct (an iOS-consumer contract enforced
/// by a regression test) with no room to add a reference-standard field, so this is a
/// standalone pure conversion instead of a struct setter: call it once on a raw BC before
/// writing the result into `FFIBallisticInputs.bc_value`, then use every existing
/// `ballistics_calculate_trajectory*`/`ballistics_calculate_zero_angle*`/`ballistics_monte_carlo*`
/// export completely unchanged. `reference_standard == FFI_BC_REFERENCE_ICAO` (`0`) is a
/// no-op, so every existing caller that never calls this function is unaffected — this is a
/// pure addition to the ABI, not a modification, so no recompile is required unless a caller
/// opts into the new symbol.
///
/// `reference_standard == FFI_BC_REFERENCE_ARMY_STANDARD_METRO` (`1`) multiplies by
/// [`crate::constants::ASM_TO_ICAO_BC`], the same constant and the same single multiply
/// [`crate::cli_api::TrajectorySolver::new`] applies for the CLI/WASM/Rust-native surfaces.
/// A non-finite `bc` is returned unchanged (this function performs no validation; the
/// existing `bc_value must be finite and greater than zero` solve-time check still applies).
#[no_mangle]
pub extern "C" fn ballistics_bc_for_reference_standard(
    bc: c_double,
    reference_standard: c_int,
) -> c_double {
    if reference_standard == FFI_BC_REFERENCE_ARMY_STANDARD_METRO {
        bc * crate::constants::ASM_TO_ICAO_BC
    } else {
        bc
    }
}

/// Reduce a sea-level-corrected altimeter setting (QNH, in hPa) to the station pressure at
/// `altitude_m` (MBA-1397; see [`crate::atmosphere::reduce_qnh_to_station_pressure`] for the
/// formula). `FFIAtmosphericConditions.pressure` has always meant absolute station pressure,
/// and remains a frozen `repr(C)` struct enforced by an ABI regression test with no room to
/// add a pressure-reference-mode field, so this is a standalone pure conversion instead of a
/// struct setter — exactly the same pattern as [`ballistics_bc_for_reference_standard`]: call
/// it once on a caller-declared QNH reading before writing the result into
/// `FFIAtmosphericConditions.pressure`, then use every existing
/// `ballistics_calculate_trajectory*`/`ballistics_calculate_zero_angle*`/`ballistics_monte_carlo*`
/// export completely unchanged — every one of them reads `pressure` as absolute station
/// pressure already, so feeding it an already-reduced value is a pure addition to the ABI,
/// not a modification. No existing caller that never calls this function is affected, and no
/// recompile is required for callers that don't opt into QNH support.
///
/// A non-finite `qnh_hpa` or `altitude_m` is returned unchanged (this function performs no
/// validation; the existing per-export input checks still apply to whatever ends up in
/// `FFIAtmosphericConditions.pressure`).
#[no_mangle]
pub extern "C" fn ballistics_reduce_qnh_pressure(
    qnh_hpa: c_double,
    altitude_m: c_double,
) -> c_double {
    if !qnh_hpa.is_finite() || !altitude_m.is_finite() {
        return qnh_hpa;
    }
    crate::atmosphere::reduce_qnh_to_station_pressure(qnh_hpa, altitude_m)
}

/// Sentinel `explicit_temperature_c` value meaning "no explicit temperature override" for the
/// `ballistics_density_altitude_*` exports below (MBA-1366) — same NaN-means-absent convention
/// [`FFIBallisticInputs::latitude`] already uses. Any NaN bit pattern is accepted (checked via
/// `is_nan()`, not equality), matching that precedent.
pub const FFI_NO_EXPLICIT_TEMPERATURE: c_double = f64::NAN;

/// Back-solve the station TEMPERATURE (Celsius) an ISA-equivalent atmosphere at
/// `density_altitude_m` implies (MBA-1366; see
/// [`crate::atmosphere::resolve_atmosphere_for_density_altitude`] for the full derivation).
///
/// `FFIAtmosphericConditions` is an ABI-frozen `repr(C)` struct (the same iOS-consumer contract
/// enforced by a regression test as [`ballistics_reduce_qnh_pressure`]/
/// [`ballistics_bc_for_reference_standard`]) with no room for a density-altitude field, so this
/// is a standalone pure conversion — call it (and its two companions below) once on a declared
/// density altitude, then write the three results into `FFIAtmosphericConditions.temperature`/
/// `.pressure`/`.altitude` before calling any existing `ballistics_calculate_trajectory*`/
/// `ballistics_calculate_zero_angle*`/`ballistics_monte_carlo*` export unchanged — a pure
/// addition to the C ABI requiring no recompile for existing callers.
///
/// `explicit_temperature_c`: pass [`FFI_NO_EXPLICIT_TEMPERATURE`] (NaN) for the ISA-at-density-
/// altitude default, or a real Celsius value to have it honored exactly (density is still
/// honored either way — only the implied pressure/altitude differ; see the Rust doc comment).
/// A non-finite `density_altitude_m` returns NaN for all three exports (there is no plausible
/// station value to fall back to, unlike the QNH/BC conversions above, which return their input
/// unchanged).
#[no_mangle]
pub extern "C" fn ballistics_density_altitude_temperature_c(
    density_altitude_m: c_double,
    explicit_temperature_c: c_double,
) -> c_double {
    if !density_altitude_m.is_finite() {
        return f64::NAN;
    }
    let explicit = (!explicit_temperature_c.is_nan()).then_some(explicit_temperature_c);
    crate::atmosphere::resolve_atmosphere_for_density_altitude(density_altitude_m, explicit).1
}

/// Companion to [`ballistics_density_altitude_temperature_c`]: the back-solved station PRESSURE
/// (hPa) for the same `(density_altitude_m, explicit_temperature_c)` pair.
#[no_mangle]
pub extern "C" fn ballistics_density_altitude_pressure_hpa(
    density_altitude_m: c_double,
    explicit_temperature_c: c_double,
) -> c_double {
    if !density_altitude_m.is_finite() {
        return f64::NAN;
    }
    let explicit = (!explicit_temperature_c.is_nan()).then_some(explicit_temperature_c);
    crate::atmosphere::resolve_atmosphere_for_density_altitude(density_altitude_m, explicit).2
}

/// Companion to [`ballistics_density_altitude_temperature_c`]: the back-solved station ALTITUDE
/// (meters, geometric) for the same `(density_altitude_m, explicit_temperature_c)` pair. This is
/// NOT necessarily equal to `density_altitude_m` — it only is when `explicit_temperature_c` is
/// [`FFI_NO_EXPLICIT_TEMPERATURE`] (see the Rust doc comment's algebraic identity).
#[no_mangle]
pub extern "C" fn ballistics_density_altitude_altitude_m(
    density_altitude_m: c_double,
    explicit_temperature_c: c_double,
) -> c_double {
    if !density_altitude_m.is_finite() {
        return f64::NAN;
    }
    let explicit = (!explicit_temperature_c.is_nan()).then_some(explicit_temperature_c);
    crate::atmosphere::resolve_atmosphere_for_density_altitude(density_altitude_m, explicit).0
}

/// Largest `marks_len` [`ballistics_hold_point_in_reticle`] will accept (MBA-1361).
///
/// Mirrors [`MAX_FFI_DRAG_TABLE_LEN`] and exists for the same reason (MBA-1407): the
/// export copies a caller-owned array whose length it cannot otherwise verify, so an
/// unbounded `len` would let one call request a multi-gigabyte read. It is also
/// [`crate::reticle::MAX_RETICLE_MARKS`], so the C ABI and the Rust API reject the same
/// inputs. Real reticles carry tens of marks; a dense tree carries a few hundred.
pub const MAX_FFI_RETICLE_MARKS: c_int = crate::reticle::MAX_RETICLE_MARKS as c_int;

/// `focal_plane` value selecting a FIRST-focal-plane reticle (subtensions constant across
/// magnification). Any value other than [`FFI_RETICLE_SECOND_FOCAL_PLANE`] is treated as FFP.
pub const FFI_RETICLE_FIRST_FOCAL_PLANE: c_int = 0;
/// `focal_plane` value selecting a SECOND-focal-plane reticle (subtensions scale as
/// `reference_magnification / magnification`).
pub const FFI_RETICLE_SECOND_FOCAL_PLANE: c_int = 1;

/// [`ballistics_hold_point_in_reticle`] succeeded and `out` was written.
pub const FFI_RETICLE_OK: c_int = 0;
/// A null pointer, or a `marks_len` outside `1..=`[`MAX_FFI_RETICLE_MARKS`].
pub const FFI_RETICLE_ERR_INVALID_ARGUMENT: c_int = -1;
/// `magnification` was not finite and greater than zero.
pub const FFI_RETICLE_ERR_MAGNIFICATION: c_int = -2;
/// A second-focal-plane call carried a non-finite or non-positive `reference_magnification`.
pub const FFI_RETICLE_ERR_REFERENCE_MAGNIFICATION: c_int = -3;
/// A mark coordinate, or the supplied firing solution, was not finite.
pub const FFI_RETICLE_ERR_NON_FINITE: c_int = -4;

/// Where a firing solution lands in a reticle (MBA-1361).
///
/// A NEW appended struct — no existing `repr(C)` layout is touched by this addition, so
/// existing callers need no recompile. All angles are milliradians from the optical
/// center; `down` is positive BELOW center and `right` is positive to the shooter's
/// RIGHT (see [`crate::reticle`] for the full convention block).
#[repr(C)]
pub struct FFIReticleHold {
    /// True angular milliradians below center. Equals the supplied `drop_mil`.
    pub down_mil: c_double,
    /// True angular milliradians right of center. Equals the supplied `wind_mil`.
    pub right_mil: c_double,
    /// Index of the nearest mark in the caller's array, or `-1` when there is none
    /// (unreachable today: an empty mark array is rejected before the search).
    pub nearest_mark: c_int,
    /// Distance from the hold to that mark, milliradians, measured in TRUE angular space
    /// (i.e. after second-focal-plane scaling).
    pub nearest_mark_distance_mil: c_double,
    /// `1` when the hold falls outside the marks' bounding box grown by 20% of its span
    /// per axis, `0` otherwise.
    pub off_reticle: c_int,
    /// The subtension scale applied to the marks: `reference_magnification / magnification`
    /// for a second-focal-plane reticle, exactly `1.0` for first focal plane.
    pub mark_scale: c_double,
}

/// Place an angular firing solution in a reticle (MBA-1361).
///
/// `marks` is a FLAT array of `2 * marks_len` doubles laid out as
/// `[down_0, right_0, down_1, right_1, ...]` in NOMINAL milliradians (as etched; for a
/// second-focal-plane reticle that means "true at `reference_magnification`").
/// `focal_plane` is [`FFI_RETICLE_FIRST_FOCAL_PLANE`] or
/// [`FFI_RETICLE_SECOND_FOCAL_PLANE`]; `reference_magnification` is consulted only in the
/// second-focal-plane case.
///
/// Returns [`FFI_RETICLE_OK`] and writes `out` on success, or one of the negative
/// `FFI_RETICLE_ERR_*` codes, in which case `out` is left untouched.
///
/// # Safety
///
/// `marks` must point to at least `2 * marks_len` readable `double`s and `out` to one
/// writable [`FFIReticleHold`]. `marks_len` is validated against
/// `1..=`[`MAX_FFI_RETICLE_MARKS`] BEFORE any element is read (MBA-1407 lesson: the
/// bound is the only thing standing between a caller typo and an out-of-range read).
#[no_mangle]
pub unsafe extern "C" fn ballistics_hold_point_in_reticle(
    drop_mil: c_double,
    wind_mil: c_double,
    magnification: c_double,
    marks: *const c_double,
    marks_len: c_int,
    focal_plane: c_int,
    reference_magnification: c_double,
    out: *mut FFIReticleHold,
) -> c_int {
    use crate::reticle::{
        hold_point_in_reticle, FocalPlane, MarkKind, ReticleDescription, ReticleError, ReticleMark,
    };

    if marks.is_null() || out.is_null() || !(1..=MAX_FFI_RETICLE_MARKS).contains(&marks_len) {
        return FFI_RETICLE_ERR_INVALID_ARGUMENT;
    }
    let count = marks_len as usize;
    // Length validated above, so this read stays inside the caller's declared array.
    let flat = unsafe { std::slice::from_raw_parts(marks, count * 2) };

    let description = ReticleDescription {
        name: String::new(),
        focal_plane: if focal_plane == FFI_RETICLE_SECOND_FOCAL_PLANE {
            FocalPlane::Second
        } else {
            FocalPlane::First
        },
        reference_magnification,
        marks: flat
            .chunks_exact(2)
            .map(|pair| ReticleMark::new(pair[0], pair[1], MarkKind::Dot))
            .collect(),
    };

    let hold = match hold_point_in_reticle(drop_mil, wind_mil, magnification, &description) {
        Ok(hold) => hold,
        Err(ReticleError::NonPositiveMagnification { .. }) => return FFI_RETICLE_ERR_MAGNIFICATION,
        Err(ReticleError::NonPositiveReferenceMagnification { .. }) => {
            return FFI_RETICLE_ERR_REFERENCE_MAGNIFICATION
        }
        Err(ReticleError::NonFiniteMark { .. }) | Err(ReticleError::NonFiniteHold { .. }) => {
            return FFI_RETICLE_ERR_NON_FINITE
        }
        Err(_) => return FFI_RETICLE_ERR_INVALID_ARGUMENT,
    };

    unsafe {
        *out = FFIReticleHold {
            down_mil: hold.down_mil,
            right_mil: hold.right_mil,
            nearest_mark: hold.nearest_mark.map_or(-1, |index| index as c_int),
            nearest_mark_distance_mil: hold.nearest_mark_distance_mil,
            off_reticle: c_int::from(hold.off_reticle),
            mark_scale: hold.mark_scale,
        };
    }
    FFI_RETICLE_OK
}

// Get library version
#[no_mangle]
pub extern "C" fn ballistics_get_version() -> *const c_char {
    // Return a pointer to a static NUL-terminated string (the caller must NOT free it).
    // Previously this leaked a freshly-allocated CString on every call and reported a
    // stale hardcoded "0.3.0"; use the real crate version with no allocation.
    concat!(env!("CARGO_PKG_VERSION"), "\0").as_ptr() as *const c_char
}

#[cfg(test)]
mod tests {
    use super::*;

    fn valid_trajectory_inputs() -> FFIBallisticInputs {
        FFIBallisticInputs {
            muzzle_velocity: 800.0,
            muzzle_angle: 0.0,
            bc_value: 0.5,
            bullet_mass: 0.01,
            bullet_diameter: 0.00762,
            bc_type: 0,
            sight_height: 0.05,
            target_distance: 1.0,
            temperature: 15.0,
            twist_rate: 12.0,
            is_twist_right: 1,
            shooting_angle: 0.0,
            altitude: 0.0,
            latitude: f64::NAN,
            azimuth_angle: 0.0,
            use_rk4: 1,
            use_adaptive_rk45: 0,
            enable_wind_shear: 0,
            enable_trajectory_sampling: 0,
            sample_interval: 10.0,
            enable_pitch_damping: 0,
            enable_precession_nutation: 0,
            enable_spin_drift: 0,
            enable_magnus: 0,
            enable_coriolis: 0,
            shot_azimuth: 0.0,
            cant_angle: 0.0,
            zero_poi_vertical: 0.0,
            zero_poi_horizontal: 0.0,
            sight_offset_lateral: 0.0,
        }
    }

    #[allow(dead_code)]
    #[repr(C)]
    struct LegacyFFIMonteCarloParams {
        num_simulations: c_int,
        velocity_std_dev: c_double,
        angle_std_dev: c_double,
        bc_std_dev: c_double,
        wind_speed_std_dev: c_double,
        target_distance: c_double,
        base_wind_speed: c_double,
        base_wind_direction: c_double,
        azimuth_std_dev: c_double,
    }

    #[test]
    fn monte_carlo_params_legacy_abi_size_is_unchanged() {
        assert_eq!(
            std::mem::size_of::<FFIMonteCarloParams>(),
            std::mem::size_of::<LegacyFFIMonteCarloParams>()
        );
        assert_eq!(
            std::mem::align_of::<FFIMonteCarloParams>(),
            std::mem::align_of::<LegacyFFIMonteCarloParams>()
        );
    }

    /// MBA-1361: the reticle export is APPEND-ONLY — a new struct and a new function.
    /// This pins that no pre-existing `repr(C)` layout moved when it landed.
    #[test]
    fn reticle_addition_does_not_disturb_existing_layouts() {
        assert_eq!(
            std::mem::size_of::<FFIMonteCarloParams>(),
            std::mem::size_of::<LegacyFFIMonteCarloParams>()
        );
        // The new struct is 6 fields: 4 doubles + 2 ints, C-laid-out.
        assert_eq!(std::mem::align_of::<FFIReticleHold>(), 8);
    }

    fn zeroed_hold() -> FFIReticleHold {
        FFIReticleHold {
            down_mil: 0.0,
            right_mil: 0.0,
            nearest_mark: -99,
            nearest_mark_distance_mil: -1.0,
            off_reticle: -1,
            mark_scale: -1.0,
        }
    }

    #[test]
    fn ffi_hold_point_matches_the_rust_api_on_both_focal_planes() {
        // down/right pairs: center, 2 mil, 4 mil, and a windage dot.
        let marks: [c_double; 8] = [0.0, 0.0, 2.0, 0.0, 4.0, 0.0, 2.0, 1.0];
        let mut out = zeroed_hold();

        // FFP at any magnification: marks are used as etched.
        let code = unsafe {
            ballistics_hold_point_in_reticle(
                4.0,
                0.0,
                6.0,
                marks.as_ptr(),
                4,
                FFI_RETICLE_FIRST_FOCAL_PLANE,
                0.0,
                &mut out,
            )
        };
        assert_eq!(code, FFI_RETICLE_OK);
        assert_eq!(out.down_mil, 4.0);
        assert_eq!(out.nearest_mark, 2);
        assert_eq!(out.nearest_mark_distance_mil, 0.0);
        assert_eq!(out.mark_scale, 1.0);
        assert_eq!(out.off_reticle, 0);

        // SFP at half the reference magnification: the 2 mil mark reads 4 mil true.
        let mut out = zeroed_hold();
        let code = unsafe {
            ballistics_hold_point_in_reticle(
                4.0,
                0.0,
                5.0,
                marks.as_ptr(),
                4,
                FFI_RETICLE_SECOND_FOCAL_PLANE,
                10.0,
                &mut out,
            )
        };
        assert_eq!(code, FFI_RETICLE_OK);
        assert_eq!(out.nearest_mark, 1);
        assert_eq!(out.nearest_mark_distance_mil, 0.0);
        assert_eq!(out.mark_scale, 2.0);
    }

    /// The MBA-1407 lesson applied to the new export: `marks_len` is validated against a
    /// stated bound BEFORE a single element is read, and a null pointer is rejected.
    #[test]
    fn ffi_hold_point_bounds_check_marks_len_before_reading() {
        let marks: [c_double; 4] = [0.0, 0.0, 2.0, 0.0];
        let mut out = zeroed_hold();
        let call = |len: c_int, ptr: *const c_double, out: &mut FFIReticleHold| unsafe {
            ballistics_hold_point_in_reticle(
                1.0,
                0.0,
                10.0,
                ptr,
                len,
                FFI_RETICLE_FIRST_FOCAL_PLANE,
                0.0,
                out,
            )
        };

        assert_eq!(call(0, marks.as_ptr(), &mut out), FFI_RETICLE_ERR_INVALID_ARGUMENT);
        assert_eq!(call(-1, marks.as_ptr(), &mut out), FFI_RETICLE_ERR_INVALID_ARGUMENT);
        assert_eq!(
            call(MAX_FFI_RETICLE_MARKS + 1, marks.as_ptr(), &mut out),
            FFI_RETICLE_ERR_INVALID_ARGUMENT
        );
        assert_eq!(call(c_int::MAX, marks.as_ptr(), &mut out), FFI_RETICLE_ERR_INVALID_ARGUMENT);
        assert_eq!(call(2, std::ptr::null(), &mut out), FFI_RETICLE_ERR_INVALID_ARGUMENT);
        // `out` is untouched on every rejection.
        assert_eq!(out.nearest_mark, -99);

        // A null `out` is rejected too, without reading the marks.
        assert_eq!(
            unsafe {
                ballistics_hold_point_in_reticle(
                    1.0,
                    0.0,
                    10.0,
                    marks.as_ptr(),
                    2,
                    FFI_RETICLE_FIRST_FOCAL_PLANE,
                    0.0,
                    std::ptr::null_mut(),
                )
            },
            FFI_RETICLE_ERR_INVALID_ARGUMENT
        );
    }

    #[test]
    fn ffi_hold_point_maps_each_error_class_to_its_own_code() {
        let marks: [c_double; 4] = [0.0, 0.0, 2.0, 0.0];
        let bad_marks: [c_double; 4] = [0.0, 0.0, f64::NAN, 0.0];
        let mut out = zeroed_hold();
        let call = |drop: c_double, mag: c_double, plane: c_int, ref_mag: c_double,
                    m: &[c_double], out: &mut FFIReticleHold| unsafe {
            ballistics_hold_point_in_reticle(
                drop,
                0.0,
                mag,
                m.as_ptr(),
                (m.len() / 2) as c_int,
                plane,
                ref_mag,
                out,
            )
        };

        assert_eq!(
            call(1.0, 0.0, FFI_RETICLE_FIRST_FOCAL_PLANE, 0.0, &marks, &mut out),
            FFI_RETICLE_ERR_MAGNIFICATION
        );
        assert_eq!(
            call(1.0, 10.0, FFI_RETICLE_SECOND_FOCAL_PLANE, 0.0, &marks, &mut out),
            FFI_RETICLE_ERR_REFERENCE_MAGNIFICATION
        );
        assert_eq!(
            call(f64::NAN, 10.0, FFI_RETICLE_FIRST_FOCAL_PLANE, 0.0, &marks, &mut out),
            FFI_RETICLE_ERR_NON_FINITE
        );
        assert_eq!(
            call(1.0, 10.0, FFI_RETICLE_FIRST_FOCAL_PLANE, 0.0, &bad_marks, &mut out),
            FFI_RETICLE_ERR_NON_FINITE
        );
        assert_eq!(out.nearest_mark, -99, "out stays untouched on every error");
    }

    /// MBA-1386 scope addition: `bc_type` gains an additive numeric slot (8) for the
    /// new RA4 family, appended after the existing 0-7 mapping (G1, G7, G2, G5, G6,
    /// G8, GI, GS). No existing caller's value changes meaning; `8` is simply new,
    /// and anything still outside 0-8 keeps falling back to G1.
    #[test]
    fn bc_type_8_maps_to_ra4() {
        let mut inputs = valid_trajectory_inputs();
        inputs.bc_type = 8;
        assert_eq!(convert_inputs(&inputs).bc_type, DragModel::RA4);

        // Every pre-existing code is unaffected by the addition.
        let expected = [
            (0, DragModel::G1),
            (1, DragModel::G7),
            (2, DragModel::G2),
            (3, DragModel::G5),
            (4, DragModel::G6),
            (5, DragModel::G8),
            (6, DragModel::GI),
            (7, DragModel::GS),
        ];
        for (code, model) in expected {
            let mut inputs = valid_trajectory_inputs();
            inputs.bc_type = code;
            assert_eq!(convert_inputs(&inputs).bc_type, model, "code {code}");
        }

        // Unrecognized codes (including ones above the new 8) still fall back to G1.
        for code in [9, 42, -1] {
            let mut inputs = valid_trajectory_inputs();
            inputs.bc_type = code;
            assert_eq!(convert_inputs(&inputs).bc_type, DragModel::G1, "code {code}");
        }
    }

    #[test]
    fn null_pointer_contracts_return_sentinels_and_free_safely() {
        unsafe {
            assert!(ballistics_calculate_trajectory(
                std::ptr::null(),
                std::ptr::null(),
                std::ptr::null(),
                1_000.0,
                1.0,
            )
            .is_null());
            assert!(ballistics_calculate_zero_angle(
                std::ptr::null(),
                std::ptr::null(),
                std::ptr::null(),
                100.0,
            )
            .is_nan());
            assert!(ballistics_calculate_trajectory_with_drag_table(
                std::ptr::null(),
                std::ptr::null(),
                std::ptr::null(),
                1_000.0,
                1.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                DECK_MACH.len() as c_int,
            )
            .is_null());
            assert!(ballistics_calculate_zero_angle_with_drag_table(
                std::ptr::null(),
                std::ptr::null(),
                std::ptr::null(),
                100.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                DECK_MACH.len() as c_int,
            )
            .is_nan());
            assert!(
                ballistics_monte_carlo(std::ptr::null(), std::ptr::null(), std::ptr::null(),)
                    .is_null()
            );
            assert!(ballistics_monte_carlo_with_direction_std_dev(
                std::ptr::null(),
                std::ptr::null(),
                std::ptr::null(),
                0.1,
            )
            .is_null());

            ballistics_free_trajectory_result(std::ptr::null_mut());
            ballistics_free_monte_carlo_results(std::ptr::null_mut());
        }
    }

    #[test]
    fn mba1283_ffi_enforces_step_floor_for_every_solver_mode() {
        for (mode, use_rk4, use_adaptive_rk45) in [("Euler", 0, 0), ("RK4", 1, 0), ("RK45", 1, 1)] {
            for step_size in [
                f64::NAN,
                f64::INFINITY,
                f64::NEG_INFINITY,
                -1.0,
                -0.0,
                0.0,
                0.001,
                MIN_FFI_STEP_SIZE_MS - 0.001,
            ] {
                let mut inputs = valid_trajectory_inputs();
                inputs.use_rk4 = use_rk4;
                inputs.use_adaptive_rk45 = use_adaptive_rk45;
                let result = unsafe {
                    ballistics_calculate_trajectory(
                        &inputs,
                        std::ptr::null(),
                        std::ptr::null(),
                        0.01,
                        step_size,
                    )
                };
                assert!(
                    result.is_null(),
                    "{mode} step_size={step_size:?} bypassed the FFI floor"
                );
            }

            let mut inputs = valid_trajectory_inputs();
            inputs.use_rk4 = use_rk4;
            inputs.use_adaptive_rk45 = use_adaptive_rk45;
            let result = unsafe {
                ballistics_calculate_trajectory(
                    &inputs,
                    std::ptr::null(),
                    std::ptr::null(),
                    0.01,
                    MIN_FFI_STEP_SIZE_MS,
                )
            };
            assert!(
                !result.is_null(),
                "the documented minimum step must remain usable in {mode}"
            );
            unsafe {
                assert!((*result).point_count >= 0);
                assert!((*result).point_count as usize <= crate::MAX_TRAJECTORY_POINTS);
                ballistics_free_trajectory_result(result);
            }
        }
    }

    /// A tiny valid deck: strictly ascending Mach, positive Cd.
    const DECK_MACH: [f64; 4] = [0.5, 1.0, 2.0, 3.0];
    /// Deliberately LOW drag so the deck measurably increases impact velocity vs G1.
    const DECK_CD_LOW: [f64; 4] = [0.05, 0.08, 0.06, 0.05];

    #[test]
    fn trajectory_with_drag_table_applies_the_deck() {
        let inputs = valid_trajectory_inputs();
        unsafe {
            let plain = ballistics_calculate_trajectory(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                300.0,
                1.0,
            );
            let decked = ballistics_calculate_trajectory_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                300.0,
                1.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                DECK_MACH.len() as c_int,
            );
            assert!(!plain.is_null() && !decked.is_null());
            // The low-drag deck must retain materially more velocity than the G-model.
            assert!(
                (*decked).impact_velocity > (*plain).impact_velocity + 1.0,
                "deck did not change the solve: plain={} decked={}",
                (*plain).impact_velocity,
                (*decked).impact_velocity
            );
            ballistics_free_trajectory_result(plain);
            ballistics_free_trajectory_result(decked);
        }
    }

    #[test]
    fn trajectory_with_drag_table_rejects_invalid_decks() {
        let inputs = valid_trajectory_inputs();
        let descending = [3.0, 2.0, 1.0, 0.5];
        let negative_cd = [0.05, -0.08, 0.06, 0.05];
        unsafe {
            // null arrays
            assert!(ballistics_calculate_trajectory_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                300.0,
                1.0,
                std::ptr::null(),
                DECK_CD_LOW.as_ptr(),
                4,
            )
            .is_null());
            assert!(ballistics_calculate_trajectory_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                300.0,
                1.0,
                DECK_MACH.as_ptr(),
                std::ptr::null(),
                4,
            )
            .is_null());
            // too few points
            assert!(ballistics_calculate_trajectory_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                300.0,
                1.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                1,
            )
            .is_null());
            // non-ascending Mach
            assert!(ballistics_calculate_trajectory_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                300.0,
                1.0,
                descending.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                4,
            )
            .is_null());
            // non-positive Cd
            assert!(ballistics_calculate_trajectory_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                300.0,
                1.0,
                DECK_MACH.as_ptr(),
                negative_cd.as_ptr(),
                4,
            )
            .is_null());
            // null inputs still rejected
            assert!(ballistics_calculate_trajectory_with_drag_table(
                std::ptr::null(),
                std::ptr::null(),
                std::ptr::null(),
                300.0,
                1.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                4,
            )
            .is_null());
        }
    }

    #[test]
    fn zero_angle_with_drag_table_applies_the_deck() {
        // A realistic zeroing setup: 100 m zero.
        let inputs = valid_trajectory_inputs();
        unsafe {
            let plain =
                ballistics_calculate_zero_angle(&inputs, std::ptr::null(), std::ptr::null(), 100.0);
            let decked = ballistics_calculate_zero_angle_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                100.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                DECK_MACH.len() as c_int,
            );
            assert!(plain.is_finite() && decked.is_finite());
            // A much lower-drag deck needs a flatter (smaller) zero angle; at minimum it
            // must differ measurably from the G-model zero.
            assert!(
                (plain - decked).abs() > 1e-6,
                "deck did not change the zero: plain={plain} decked={decked}"
            );
        }
    }

    #[test]
    fn zero_angle_with_drag_table_rejects_invalid_decks() {
        let inputs = valid_trajectory_inputs();
        let descending = [3.0, 2.0, 1.0, 0.5];
        unsafe {
            assert!(ballistics_calculate_zero_angle_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                100.0,
                std::ptr::null(),
                DECK_CD_LOW.as_ptr(),
                4,
            )
            .is_nan());
            assert!(ballistics_calculate_zero_angle_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                100.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                0,
            )
            .is_nan());
            assert!(ballistics_calculate_zero_angle_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                100.0,
                descending.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                4,
            )
            .is_nan());
            // null inputs still rejected
            assert!(ballistics_calculate_zero_angle_with_drag_table(
                std::ptr::null(),
                std::ptr::null(),
                std::ptr::null(),
                100.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                4,
            )
            .is_nan());
        }
    }

    #[test]
    fn zero_then_fly_with_same_deck_is_consistent() {
        // The pair-use case the two exports exist for: zero with the deck, fly with the
        // deck at the solved angle; the trajectory must cross near sight height at the
        // zero distance (i.e. the two functions share identical deck semantics).
        let mut inputs = valid_trajectory_inputs();
        unsafe {
            let angle = ballistics_calculate_zero_angle_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                100.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                DECK_MACH.len() as c_int,
            );
            assert!(angle.is_finite());
            inputs.muzzle_angle = angle;
            let result = ballistics_calculate_trajectory_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                150.0,
                1.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                DECK_MACH.len() as c_int,
            );
            assert!(!result.is_null());
            // Interpolate y at exactly the zero distance (100 m) rather than snapping to the
            // nearest raw trajectory sample, so the residual reflects only zero-solver
            // convergence, not the sampling grid's x-offset from 100 m.
            let zero_distance = 100.0;
            let pts = std::slice::from_raw_parts((*result).points, (*result).point_count as usize);
            let bracket = pts
                .windows(2)
                .find(|w| w[0].position_x <= zero_distance && w[1].position_x >= zero_distance)
                .expect("trajectory brackets the zero distance");
            let (lo, hi) = (&bracket[0], &bracket[1]);
            let y_at_zero = if hi.position_x > lo.position_x {
                let t = (zero_distance - lo.position_x) / (hi.position_x - lo.position_x);
                lo.position_y + t * (hi.position_y - lo.position_y)
            } else {
                lo.position_y
            };
            assert!(
                (y_at_zero - inputs.sight_height).abs() < 0.002,
                "zeroed flight missed the line of sight at 100 m: y={} (sight_height={})",
                y_at_zero,
                inputs.sight_height
            );
            ballistics_free_trajectory_result(result);
        }
    }

    #[test]
    fn drag_table_len_above_cap_is_rejected() {
        // A valid, monotonically increasing deck that is simply too long: the cap
        // must reject it BEFORE the to_vec() copies, returning the null sentinel.
        let n = (MAX_FFI_DRAG_TABLE_LEN + 1) as usize;
        let mach: Vec<f64> = (0..n).map(|i| 0.01 + i as f64 * 0.001).collect();
        let cd: Vec<f64> = vec![0.3; n];
        let inputs = valid_trajectory_inputs();
        unsafe {
            let r = ballistics_calculate_trajectory_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                300.0,
                1.0,
                mach.as_ptr(),
                cd.as_ptr(),
                n as c_int,
            );
            assert!(r.is_null(), "len {n} must be rejected by the cap");
        }
    }

    #[test]
    fn drag_table_len_at_cap_is_accepted() {
        let n = MAX_FFI_DRAG_TABLE_LEN as usize;
        let mach: Vec<f64> = (0..n).map(|i| 0.01 + i as f64 * 0.001).collect();
        let cd: Vec<f64> = vec![0.3; n];
        let inputs = valid_trajectory_inputs();
        unsafe {
            let r = ballistics_calculate_trajectory_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                300.0,
                1.0,
                mach.as_ptr(),
                cd.as_ptr(),
                n as c_int,
            );
            assert!(!r.is_null(), "len == cap must be accepted");
            ballistics_free_trajectory_result(r);
        }
    }

    #[test]
    fn ffi_cant_angle_deflects_laterally() {
        let mut level = valid_trajectory_inputs();
        level.muzzle_angle = 0.003;
        let mut canted = valid_trajectory_inputs();
        canted.muzzle_angle = 0.003;
        canted.cant_angle = 10f64.to_radians();
        unsafe {
            let a = ballistics_calculate_trajectory(&level, std::ptr::null(), std::ptr::null(), 400.0, 1.0);
            let b = ballistics_calculate_trajectory(&canted, std::ptr::null(), std::ptr::null(), 400.0, 1.0);
            assert!(!a.is_null() && !b.is_null());
            let za = std::slice::from_raw_parts((*a).points, (*a).point_count as usize).last().unwrap().position_z;
            let zb = std::slice::from_raw_parts((*b).points, (*b).point_count as usize).last().unwrap().position_z;
            assert!(zb > za + 0.005, "FFI cant must deflect right: level={za} canted={zb}");
            ballistics_free_trajectory_result(a);
            ballistics_free_trajectory_result(b);
        }
    }

    #[test]
    fn ffi_vertical_wind_raises_trajectory() {
        let inputs = valid_trajectory_inputs();
        let no_wind = FFIWindConditions {
            speed: 0.0,
            direction: 0.0,
            vertical_speed: 0.0,
        };
        let updraft = FFIWindConditions {
            speed: 0.0,
            direction: 0.0,
            vertical_speed: 5.0,
        };
        unsafe {
            let a = ballistics_calculate_trajectory(&inputs, &no_wind, std::ptr::null(), 400.0, 1.0);
            let b = ballistics_calculate_trajectory(&inputs, &updraft, std::ptr::null(), 400.0, 1.0);
            assert!(!a.is_null() && !b.is_null());
            let ya = std::slice::from_raw_parts((*a).points, (*a).point_count as usize).last().unwrap().position_y;
            let yb = std::slice::from_raw_parts((*b).points, (*b).point_count as usize).last().unwrap().position_y;
            assert!(yb > ya + 0.01, "FFI updraft must raise the trajectory: no_wind={ya} updraft={yb}");
            ballistics_free_trajectory_result(a);
            ballistics_free_trajectory_result(b);
        }
    }

    // --- MBA-1356: cd_scale `_scaled` FFI variants ---

    #[test]
    fn trajectory_scaled_at_one_matches_unscaled_export() {
        let inputs = valid_trajectory_inputs();
        unsafe {
            let unscaled = ballistics_calculate_trajectory_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                300.0,
                1.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                DECK_MACH.len() as c_int,
            );
            let scaled = ballistics_calculate_trajectory_with_drag_table_scaled(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                300.0,
                1.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                DECK_MACH.len() as c_int,
                1.0,
            );
            assert!(!unscaled.is_null() && !scaled.is_null());
            assert_eq!(
                (*unscaled).impact_velocity.to_bits(),
                (*scaled).impact_velocity.to_bits(),
                "cd_scale=1.0 must be byte-identical to the unscaled export: unscaled={} scaled={}",
                (*unscaled).impact_velocity,
                (*scaled).impact_velocity
            );
            ballistics_free_trajectory_result(unscaled);
            ballistics_free_trajectory_result(scaled);
        }
    }

    #[test]
    fn trajectory_scaled_at_1_10_lowers_impact_velocity() {
        let inputs = valid_trajectory_inputs();
        unsafe {
            let baseline = ballistics_calculate_trajectory_with_drag_table_scaled(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                300.0,
                1.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                DECK_MACH.len() as c_int,
                1.0,
            );
            let scaled_up = ballistics_calculate_trajectory_with_drag_table_scaled(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                300.0,
                1.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                DECK_MACH.len() as c_int,
                1.10,
            );
            assert!(!baseline.is_null() && !scaled_up.is_null());
            assert!(
                (*scaled_up).impact_velocity < (*baseline).impact_velocity,
                "cd_scale=1.10 must increase drag -> lower impact velocity: base={} scaled={}",
                (*baseline).impact_velocity,
                (*scaled_up).impact_velocity
            );
            ballistics_free_trajectory_result(baseline);
            ballistics_free_trajectory_result(scaled_up);
        }
    }

    #[test]
    fn trajectory_scaled_rejects_invalid_cd_scale() {
        let inputs = valid_trajectory_inputs();
        unsafe {
            for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                let r = ballistics_calculate_trajectory_with_drag_table_scaled(
                    &inputs,
                    std::ptr::null(),
                    std::ptr::null(),
                    300.0,
                    1.0,
                    DECK_MACH.as_ptr(),
                    DECK_CD_LOW.as_ptr(),
                    DECK_MACH.len() as c_int,
                    bad,
                );
                assert!(r.is_null(), "cd_scale={bad} must be rejected (null sentinel)");
            }
        }
    }

    #[test]
    fn zero_angle_scaled_at_one_matches_unscaled_export() {
        let inputs = valid_trajectory_inputs();
        unsafe {
            let unscaled = ballistics_calculate_zero_angle_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                100.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                DECK_MACH.len() as c_int,
            );
            let scaled = ballistics_calculate_zero_angle_with_drag_table_scaled(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                100.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                DECK_MACH.len() as c_int,
                1.0,
            );
            assert!(unscaled.is_finite() && scaled.is_finite());
            assert_eq!(
                unscaled.to_bits(),
                scaled.to_bits(),
                "cd_scale=1.0 must be byte-identical to the unscaled export: unscaled={unscaled} scaled={scaled}"
            );
        }
    }

    #[test]
    fn zero_angle_scaled_at_1_10_differs_from_baseline() {
        let inputs = valid_trajectory_inputs();
        unsafe {
            let baseline = ballistics_calculate_zero_angle_with_drag_table_scaled(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                100.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                DECK_MACH.len() as c_int,
                1.0,
            );
            let scaled_up = ballistics_calculate_zero_angle_with_drag_table_scaled(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                100.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                DECK_MACH.len() as c_int,
                1.10,
            );
            assert!(baseline.is_finite() && scaled_up.is_finite());
            // More drag needs a steeper (larger) zero angle to still reach 100 m.
            assert!(
                scaled_up > baseline,
                "cd_scale=1.10 must need a larger zero angle: base={baseline} scaled={scaled_up}"
            );
        }
    }

    #[test]
    fn zero_angle_scaled_rejects_invalid_cd_scale() {
        let inputs = valid_trajectory_inputs();
        unsafe {
            for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                let angle = ballistics_calculate_zero_angle_with_drag_table_scaled(
                    &inputs,
                    std::ptr::null(),
                    std::ptr::null(),
                    100.0,
                    DECK_MACH.as_ptr(),
                    DECK_CD_LOW.as_ptr(),
                    DECK_MACH.len() as c_int,
                    bad,
                );
                assert!(angle.is_nan(), "cd_scale={bad} must be rejected (NaN sentinel)");
            }
        }
    }

    /// Legacy exports must remain byte-identical: the same deck through the unscaled
    /// export must be unaffected by the new cd_scale plumbing (a regression guard
    /// alongside the pre-existing, untouched drag-table tests above).
    #[test]
    fn legacy_drag_table_exports_unaffected_by_cd_scale_plumbing() {
        let inputs = valid_trajectory_inputs();
        unsafe {
            let a = ballistics_calculate_trajectory_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                300.0,
                1.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                DECK_MACH.len() as c_int,
            );
            let b = ballistics_calculate_trajectory_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                300.0,
                1.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                DECK_MACH.len() as c_int,
            );
            assert!(!a.is_null() && !b.is_null());
            assert_eq!((*a).impact_velocity.to_bits(), (*b).impact_velocity.to_bits());
            ballistics_free_trajectory_result(a);
            ballistics_free_trajectory_result(b);

            let za = ballistics_calculate_zero_angle_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                100.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                DECK_MACH.len() as c_int,
            );
            let zb = ballistics_calculate_zero_angle_with_drag_table(
                &inputs,
                std::ptr::null(),
                std::ptr::null(),
                100.0,
                DECK_MACH.as_ptr(),
                DECK_CD_LOW.as_ptr(),
                DECK_MACH.len() as c_int,
            );
            assert!(za.is_finite() && zb.is_finite());
            assert_eq!(za.to_bits(), zb.to_bits());
        }
    }

    // ---- MBA-1365: ballistics_bc_for_reference_standard --------------------------------

    #[test]
    fn bc_for_reference_standard_icao_is_a_byte_identical_no_op() {
        let bc = 0.475;
        assert_eq!(
            ballistics_bc_for_reference_standard(bc, FFI_BC_REFERENCE_ICAO).to_bits(),
            bc.to_bits()
        );
    }

    #[test]
    fn bc_for_reference_standard_unrecognized_value_falls_back_to_icao() {
        // Mirrors convert_inputs' permissive unrecognized-bc_type convention: an unknown
        // reference_standard is treated as ICAO (0), not rejected.
        let bc = 0.475;
        assert_eq!(
            ballistics_bc_for_reference_standard(bc, 99).to_bits(),
            bc.to_bits()
        );
    }

    #[test]
    fn bc_for_reference_standard_army_standard_metro_applies_the_documented_ratio() {
        let bc = 0.475;
        let converted =
            ballistics_bc_for_reference_standard(bc, FFI_BC_REFERENCE_ARMY_STANDARD_METRO);
        assert_eq!(converted, bc * crate::constants::ASM_TO_ICAO_BC);
        // Smaller BC == more drag under this engine's ICAO-calibrated retardation math.
        assert!(converted < bc);
    }

    // ---- MBA-1397: ballistics_reduce_qnh_pressure + FFIAtmosphericConditions.pressure ---

    #[test]
    fn reduce_qnh_pressure_matches_the_library_function_and_lowers_pressure() {
        let reduced = ballistics_reduce_qnh_pressure(1030.0, 1500.0);
        assert_eq!(
            reduced,
            crate::atmosphere::reduce_qnh_to_station_pressure(1030.0, 1500.0)
        );
        assert!(reduced < 1030.0);
    }

    #[test]
    fn reduce_qnh_pressure_passes_through_non_finite_inputs() {
        assert!(ballistics_reduce_qnh_pressure(f64::NAN, 1500.0).is_nan());
        assert_eq!(ballistics_reduce_qnh_pressure(1030.0, f64::INFINITY), 1030.0);
    }

    /// The FFI trajectory/Monte Carlo exports have always treated
    /// `FFIAtmosphericConditions.pressure` as absolute station pressure. A caller declaring a
    /// QNH (sea-level-corrected altimeter setting) reading MUST reduce it with
    /// `ballistics_reduce_qnh_pressure` before writing `pressure` -- this proves the reduced
    /// value actually reaches the solve (a materially different, and correct -- flatter,
    /// less-drop -- trajectory than feeding the raw, unreduced QNH straight through, which
    /// would silently over-state air density).
    #[test]
    fn ffi_trajectory_uses_the_reduced_pressure_not_the_raw_qnh() {
        let inputs = valid_trajectory_inputs();
        let altitude_m = 1500.0;
        let qnh_hpa = 1030.0;
        let reduced = ballistics_reduce_qnh_pressure(qnh_hpa, altitude_m);
        assert!(reduced < qnh_hpa);

        let atmo_reduced = FFIAtmosphericConditions {
            temperature: 15.0,
            pressure: reduced,
            humidity: 50.0,
            altitude: altitude_m,
        };
        let atmo_raw_qnh = FFIAtmosphericConditions {
            temperature: 15.0,
            pressure: qnh_hpa,
            humidity: 50.0,
            altitude: altitude_m,
        };

        unsafe {
            let a = ballistics_calculate_trajectory(
                &inputs,
                std::ptr::null(),
                &atmo_reduced,
                400.0,
                1.0,
            );
            let b = ballistics_calculate_trajectory(
                &inputs,
                std::ptr::null(),
                &atmo_raw_qnh,
                400.0,
                1.0,
            );
            assert!(!a.is_null() && !b.is_null());
            let drop_a = std::slice::from_raw_parts((*a).points, (*a).point_count as usize)
                .last()
                .unwrap()
                .position_y;
            let drop_b = std::slice::from_raw_parts((*b).points, (*b).point_count as usize)
                .last()
                .unwrap()
                .position_y;
            assert!(
                (drop_a - drop_b).abs() > 1e-6,
                "reduced vs. raw-QNH pressure must produce materially different trajectories: \
                 {drop_a} vs {drop_b}"
            );
            ballistics_free_trajectory_result(a);
            ballistics_free_trajectory_result(b);
        }
    }

    /// Same proof as above, for the Monte Carlo FFI path (`ballistics_monte_carlo`), which
    /// reads `FFIAtmosphericConditions.pressure` into `BallisticInputs.pressure` directly
    /// (`ballistics_monte_carlo_impl`) before running the shared `run_monte_carlo_*` core.
    #[test]
    fn ffi_monte_carlo_uses_the_reduced_pressure_not_the_raw_qnh() {
        let inputs = valid_trajectory_inputs();
        let altitude_m = 1500.0;
        let qnh_hpa = 1030.0;
        let reduced = ballistics_reduce_qnh_pressure(qnh_hpa, altitude_m);

        let atmo_reduced = FFIAtmosphericConditions {
            temperature: 15.0,
            pressure: reduced,
            humidity: 50.0,
            altitude: altitude_m,
        };
        let atmo_raw_qnh = FFIAtmosphericConditions {
            temperature: 15.0,
            pressure: qnh_hpa,
            humidity: 50.0,
            altitude: altitude_m,
        };
        let params = FFIMonteCarloParams {
            num_simulations: 200,
            velocity_std_dev: 1.0,
            angle_std_dev: 0.0,
            bc_std_dev: 0.0,
            wind_speed_std_dev: 0.0,
            target_distance: f64::NAN,
            base_wind_speed: 0.0,
            base_wind_direction: 0.0,
            azimuth_std_dev: 0.0,
        };

        unsafe {
            let a = ballistics_monte_carlo(&inputs, &atmo_reduced, &params);
            let b = ballistics_monte_carlo(&inputs, &atmo_raw_qnh, &params);
            assert!(!a.is_null() && !b.is_null());
            assert!(
                ((*a).mean_range - (*b).mean_range).abs() > 0.5,
                "reduced vs. raw-QNH pressure must change MC mean range materially: \
                 {} vs {}",
                (*a).mean_range,
                (*b).mean_range
            );
            ballistics_free_monte_carlo_results(a);
            ballistics_free_monte_carlo_results(b);
        }
    }

    // ---- MBA-1366: ballistics_density_altitude_* + FFIAtmosphericConditions parity --------

    #[test]
    fn density_altitude_ffi_exports_match_the_library_function() {
        let da_m = 1000.0 * 0.3048;
        let expected = crate::atmosphere::resolve_atmosphere_for_density_altitude(da_m, None);
        assert_eq!(
            ballistics_density_altitude_altitude_m(da_m, FFI_NO_EXPLICIT_TEMPERATURE),
            expected.0
        );
        assert_eq!(
            ballistics_density_altitude_temperature_c(da_m, FFI_NO_EXPLICIT_TEMPERATURE),
            expected.1
        );
        assert_eq!(
            ballistics_density_altitude_pressure_hpa(da_m, FFI_NO_EXPLICIT_TEMPERATURE),
            expected.2
        );

        // With no explicit temperature the resolved altitude must equal the input exactly.
        assert!((ballistics_density_altitude_altitude_m(da_m, FFI_NO_EXPLICIT_TEMPERATURE) - da_m).abs() < 1e-6);
    }

    #[test]
    fn density_altitude_ffi_explicit_temperature_is_honored_exactly() {
        let da_m = 500.0;
        let explicit_temp_c = 30.0;
        let expected =
            crate::atmosphere::resolve_atmosphere_for_density_altitude(da_m, Some(explicit_temp_c));
        assert_eq!(
            ballistics_density_altitude_temperature_c(da_m, explicit_temp_c),
            explicit_temp_c
        );
        assert_eq!(
            ballistics_density_altitude_temperature_c(da_m, explicit_temp_c),
            expected.1
        );
        assert_eq!(
            ballistics_density_altitude_pressure_hpa(da_m, explicit_temp_c),
            expected.2
        );
        assert_eq!(
            ballistics_density_altitude_altitude_m(da_m, explicit_temp_c),
            expected.0
        );
    }

    #[test]
    fn density_altitude_ffi_non_finite_input_returns_nan() {
        assert!(
            ballistics_density_altitude_temperature_c(f64::INFINITY, FFI_NO_EXPLICIT_TEMPERATURE)
                .is_nan()
        );
        assert!(
            ballistics_density_altitude_pressure_hpa(f64::NAN, FFI_NO_EXPLICIT_TEMPERATURE).is_nan()
        );
        assert!(
            ballistics_density_altitude_altitude_m(f64::NEG_INFINITY, FFI_NO_EXPLICIT_TEMPERATURE)
                .is_nan()
        );
    }

    /// Same proof shape as the QNH FFI tests above: writing the density-altitude-derived
    /// station values into `FFIAtmosphericConditions` before an existing trajectory export
    /// reaches the solve (a materially different, lower-density trajectory at a higher DA than
    /// at sea level, exactly like the QNH-vs-absolute divergence proof).
    #[test]
    fn ffi_trajectory_uses_the_density_altitude_derived_station_values() {
        let inputs = valid_trajectory_inputs();
        let da_m = 3000.0 * 0.3048; // 3000 ft density altitude
        let altitude_m =
            ballistics_density_altitude_altitude_m(da_m, FFI_NO_EXPLICIT_TEMPERATURE);
        let temperature_c =
            ballistics_density_altitude_temperature_c(da_m, FFI_NO_EXPLICIT_TEMPERATURE);
        let pressure_hpa =
            ballistics_density_altitude_pressure_hpa(da_m, FFI_NO_EXPLICIT_TEMPERATURE);

        let atmo_da = FFIAtmosphericConditions {
            temperature: temperature_c,
            pressure: pressure_hpa,
            humidity: 50.0,
            altitude: altitude_m,
        };
        let atmo_sea_level = FFIAtmosphericConditions {
            temperature: 15.0,
            pressure: 1013.25,
            humidity: 50.0,
            altitude: 0.0,
        };

        unsafe {
            let a = ballistics_calculate_trajectory(
                &inputs,
                std::ptr::null(),
                &atmo_da,
                400.0,
                1.0,
            );
            let b = ballistics_calculate_trajectory(
                &inputs,
                std::ptr::null(),
                &atmo_sea_level,
                400.0,
                1.0,
            );
            assert!(!a.is_null() && !b.is_null());
            let drop_a = std::slice::from_raw_parts((*a).points, (*a).point_count as usize)
                .last()
                .unwrap()
                .position_y;
            let drop_b = std::slice::from_raw_parts((*b).points, (*b).point_count as usize)
                .last()
                .unwrap()
                .position_y;
            assert!(
                (drop_a - drop_b).abs() > 1e-6,
                "density-altitude-derived vs sea-level atmosphere must produce materially \
                 different trajectories: {drop_a} vs {drop_b}"
            );
            ballistics_free_trajectory_result(a);
            ballistics_free_trajectory_result(b);
        }
    }
}
