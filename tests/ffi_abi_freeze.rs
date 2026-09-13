//! Freezes the C ABI: the exported symbol set, each export's signature, and the layout of
//! every `repr(C)` type those signatures point at.
//!
//! The `ballistics_*` exports are a shipped contract. iOS/Swift consume them through an
//! xcframework, and downstream callers link the staticlib/cdylib against headers they wrote by
//! hand. Nothing in the crate asserted any of that shape, so a rename, a dropped export, an
//! extra parameter, or one new field in the middle of `FFIBallisticInputs` would compile,
//! test, and ship green while silently breaking or misreading every existing caller.
//!
//! Three layers, because no single one of them sees everything:
//!
//! 1. `exported_symbol_set_is_frozen` scans the source text for `#[no_mangle]` items and
//!    diffs the names against the list below. This is the only layer that notices an export
//!    being ADDED, and the only one that notices `#[no_mangle]` being dropped from a function
//!    that still exists (which un-exports the symbol while leaving every Rust reference to it
//!    compiling fine).
//! 2. `FROZEN_*` below coerce each export to an explicitly written `extern "C"` function
//!    pointer. These are compile-time: a rename, a changed parameter type or count, a changed
//!    return type, or a safe/`unsafe` flip fails the BUILD, at the export, before any test
//!    runs.
//! 3. `frozen_repr_c!` pins each `repr(C)` struct against a mirror declaration -- size, align,
//!    and every field's byte offset by name. Layer 2 cannot see through a pointer: adding a
//!    field to `FFIBallisticInputs` leaves `*const FFIBallisticInputs` looking identical while
//!    every existing caller's struct is reinterpreted. Comparing against a mirror rather than
//!    against hardcoded byte counts keeps this honest on the 32-bit targets in the build
//!    fleet, where pointer-bearing structs are legitimately a different size.
//!
//! Adding an export is a normal, backward-compatible thing to do -- append the name to
//! `FROZEN_EXPORTS` and add its signature constant. Renaming, removing, or changing the shape
//! of an existing one is not: it is a breaking change to a published contract and the failure
//! here is the prompt to treat it as one.
//!
//! `tests/test_ffi.c` remains the complementary check from the other side (it re-declares the
//! structs field-by-field in C and is exercised by the gcc invocation in CLAUDE.md); this file
//! is what runs in `cargo test`.

#![cfg(feature = "ffi")]

use std::mem::{align_of, offset_of, size_of};
use std::os::raw::{c_char, c_double, c_int};
use std::path::Path;

use ballistics_engine::ffi::{
    FFIAtmosphericConditions, FFIBallisticInputs, FFIMonteCarloParams, FFIMonteCarloResults,
    FFIReticleHold, FFITrajectoryPoint, FFITrajectoryResult, FFITrajectorySample,
    FFIWindConditions,
};

// ---------------------------------------------------------------------------------------
// Layer 1: the exported symbol set
// ---------------------------------------------------------------------------------------

/// Every `#[no_mangle]` symbol this crate exports, sorted (the test asserts the sort, so the
/// diff on a change stays readable). `src/ffi.rs` holds the 18-export numeric C ABI;
/// `src/bridge/ffi.rs` holds the three `ballistics_bridge_*` symbols of the JSON bridge, the
/// mobile embedding surface. Both are the same published contract and both are frozen here.
const FROZEN_EXPORTS: &[&str] = &[
    "ballistics_bc_for_reference_standard",
    "ballistics_bridge_call",
    "ballistics_bridge_call_n",
    "ballistics_bridge_free",
    "ballistics_calculate_trajectory",
    "ballistics_calculate_trajectory_with_drag_table",
    "ballistics_calculate_trajectory_with_drag_table_scaled",
    "ballistics_calculate_zero_angle",
    "ballistics_calculate_zero_angle_with_drag_table",
    "ballistics_calculate_zero_angle_with_drag_table_scaled",
    "ballistics_density_altitude_altitude_m",
    "ballistics_density_altitude_pressure_hpa",
    "ballistics_density_altitude_temperature_c",
    "ballistics_free_monte_carlo_results",
    "ballistics_free_trajectory_result",
    "ballistics_get_version",
    "ballistics_hold_point_in_reticle",
    "ballistics_monte_carlo",
    "ballistics_monte_carlo_with_direction_std_dev",
    "ballistics_quick_trajectory",
    "ballistics_reduce_qnh_pressure",
];

/// The files that may contain `#[no_mangle]`. Listed explicitly rather than walked: a new file
/// exporting C symbols is itself a contract change that should land deliberately, and
/// `symbol_exports_live_only_in_the_known_files` fails if one appears elsewhere.
const EXPORTING_FILES: &[&str] = &["src/ffi.rs", "src/bridge/ffi.rs"];

/// Pull the exported symbol names out of one source file: for each `#[no_mangle]`, the name of
/// the next `extern "C" fn` declared after it.
///
/// Deliberately a text scan and not a parse. The thing being frozen is what the LINKER sees,
/// and the attribute-to-item association is what produces that; a scan keeps this test from
/// needing a syn dependency to guard a file that changes a few times a year.
fn exported_symbols_in(source: &str) -> Vec<String> {
    let mut names = Vec::new();
    let mut armed = false;
    for line in source.lines() {
        let line = line.trim();
        // `#[cfg(test)]` mirror structs and doc comments mentioning the attribute must not arm
        // the scan, so require the attribute to be the whole line.
        if line == "#[no_mangle]" || line == "#[unsafe(no_mangle)]" {
            armed = true;
            continue;
        }
        if !armed {
            continue;
        }
        if let Some(rest) = line.split_once("extern \"C\" fn ") {
            let name: String = rest
                .1
                .chars()
                .take_while(|c| c.is_alphanumeric() || *c == '_')
                .collect();
            assert!(
                !name.is_empty(),
                "could not read an export name out of: {line}"
            );
            names.push(name);
            armed = false;
        }
        // Anything else between the attribute and the function (other attributes, doc
        // comments, lint suppressions) is skipped without disarming -- `#[no_mangle]` is not
        // always the last attribute on these items.
    }
    assert!(
        !armed,
        "a #[no_mangle] in this file is not followed by an `extern \"C\" fn`"
    );
    names
}

fn read_source(relative: &str) -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join(relative);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
}

#[test]
fn exported_symbol_set_is_frozen() {
    let mut found: Vec<String> = EXPORTING_FILES
        .iter()
        .flat_map(|f| exported_symbols_in(&read_source(f)))
        .collect();
    found.sort();

    let expected: Vec<String> = FROZEN_EXPORTS.iter().map(|s| s.to_string()).collect();
    assert!(
        expected.windows(2).all(|w| w[0] < w[1]),
        "FROZEN_EXPORTS must stay sorted and free of duplicates"
    );

    let added: Vec<&String> = found.iter().filter(|n| !expected.contains(n)).collect();
    let removed: Vec<&String> = expected.iter().filter(|n| !found.contains(n)).collect();

    assert!(
        added.is_empty() && removed.is_empty(),
        "the exported C symbol set changed.\n  \
         added (new symbols not in FROZEN_EXPORTS): {added:?}\n  \
         removed (frozen symbols no longer exported, or no longer #[no_mangle]): {removed:?}\n\
         Adding an export is backward-compatible -- list it in FROZEN_EXPORTS and give it a \
         FROZEN_* signature constant. Removing or renaming one BREAKS every linked caller \
         (iOS xcframework, JNI, hand-written headers) and needs to be a deliberate, announced \
         change, not a green test."
    );
}

#[test]
fn symbol_exports_live_only_in_the_known_files() {
    let src = Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
    let mut offenders = Vec::new();
    let mut stack = vec![src.clone()];
    while let Some(dir) = stack.pop() {
        let Ok(entries) = std::fs::read_dir(&dir) else {
            continue;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                stack.push(path);
                continue;
            }
            if path.extension().and_then(|e| e.to_str()) != Some("rs") {
                continue;
            }
            let relative = path
                .strip_prefix(src.parent().expect("src has a parent"))
                .expect("path is under the crate root")
                .to_string_lossy()
                .replace('\\', "/");
            if EXPORTING_FILES.contains(&relative.as_str()) {
                continue;
            }
            let text = std::fs::read_to_string(&path).unwrap_or_default();
            if text.lines().any(|l| l.trim() == "#[no_mangle]") {
                offenders.push(relative);
            }
        }
    }
    assert!(
        offenders.is_empty(),
        "these files export C symbols but are not listed in EXPORTING_FILES, so \
         exported_symbol_set_is_frozen never looked at them: {offenders:?}"
    );
}

// ---------------------------------------------------------------------------------------
// Layer 2: signatures
//
// Each constant coerces the export to an explicitly spelled function-pointer type. These are
// checked by the COMPILER, so a signature drift is a build error at the offending export
// rather than a failing assertion.
// ---------------------------------------------------------------------------------------

const FROZEN_CALCULATE_TRAJECTORY: unsafe extern "C" fn(
    *const FFIBallisticInputs,
    *const FFIWindConditions,
    *const FFIAtmosphericConditions,
    c_double,
    c_double,
) -> *mut FFITrajectoryResult = ballistics_engine::ffi::ballistics_calculate_trajectory;

const FROZEN_CALCULATE_TRAJECTORY_WITH_DRAG_TABLE: unsafe extern "C" fn(
    *const FFIBallisticInputs,
    *const FFIWindConditions,
    *const FFIAtmosphericConditions,
    c_double,
    c_double,
    *const c_double,
    *const c_double,
    c_int,
) -> *mut FFITrajectoryResult =
    ballistics_engine::ffi::ballistics_calculate_trajectory_with_drag_table;

#[allow(clippy::type_complexity)]
const FROZEN_CALCULATE_TRAJECTORY_WITH_DRAG_TABLE_SCALED: unsafe extern "C" fn(
    *const FFIBallisticInputs,
    *const FFIWindConditions,
    *const FFIAtmosphericConditions,
    c_double,
    c_double,
    *const c_double,
    *const c_double,
    c_int,
    c_double,
) -> *mut FFITrajectoryResult =
    ballistics_engine::ffi::ballistics_calculate_trajectory_with_drag_table_scaled;

const FROZEN_FREE_TRAJECTORY_RESULT: unsafe extern "C" fn(*mut FFITrajectoryResult) =
    ballistics_engine::ffi::ballistics_free_trajectory_result;

const FROZEN_CALCULATE_ZERO_ANGLE: unsafe extern "C" fn(
    *const FFIBallisticInputs,
    *const FFIWindConditions,
    *const FFIAtmosphericConditions,
    c_double,
) -> c_double = ballistics_engine::ffi::ballistics_calculate_zero_angle;

const FROZEN_CALCULATE_ZERO_ANGLE_WITH_DRAG_TABLE: unsafe extern "C" fn(
    *const FFIBallisticInputs,
    *const FFIWindConditions,
    *const FFIAtmosphericConditions,
    c_double,
    *const c_double,
    *const c_double,
    c_int,
) -> c_double = ballistics_engine::ffi::ballistics_calculate_zero_angle_with_drag_table;

#[allow(clippy::type_complexity)]
const FROZEN_CALCULATE_ZERO_ANGLE_WITH_DRAG_TABLE_SCALED: unsafe extern "C" fn(
    *const FFIBallisticInputs,
    *const FFIWindConditions,
    *const FFIAtmosphericConditions,
    c_double,
    *const c_double,
    *const c_double,
    c_int,
    c_double,
) -> c_double = ballistics_engine::ffi::ballistics_calculate_zero_angle_with_drag_table_scaled;

/// Safe, not `unsafe` -- pinned as declared, so adding `unsafe` (a source-compatibility break
/// for Rust callers) is caught too. A safe `extern "C" fn` would coerce INTO an `unsafe`
/// pointer type silently, which is why the safety is written out here rather than blanket
/// `unsafe` everywhere.
const FROZEN_QUICK_TRAJECTORY: extern "C" fn(
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
) -> c_double = ballistics_engine::ffi::ballistics_quick_trajectory;

const FROZEN_MONTE_CARLO: unsafe extern "C" fn(
    *const FFIBallisticInputs,
    *const FFIAtmosphericConditions,
    *const FFIMonteCarloParams,
) -> *mut FFIMonteCarloResults = ballistics_engine::ffi::ballistics_monte_carlo;

const FROZEN_MONTE_CARLO_WITH_DIRECTION_STD_DEV: unsafe extern "C" fn(
    *const FFIBallisticInputs,
    *const FFIAtmosphericConditions,
    *const FFIMonteCarloParams,
    c_double,
) -> *mut FFIMonteCarloResults =
    ballistics_engine::ffi::ballistics_monte_carlo_with_direction_std_dev;

const FROZEN_FREE_MONTE_CARLO_RESULTS: unsafe extern "C" fn(*mut FFIMonteCarloResults) =
    ballistics_engine::ffi::ballistics_free_monte_carlo_results;

const FROZEN_BC_FOR_REFERENCE_STANDARD: extern "C" fn(c_double, c_int) -> c_double =
    ballistics_engine::ffi::ballistics_bc_for_reference_standard;

const FROZEN_REDUCE_QNH_PRESSURE: extern "C" fn(c_double, c_double) -> c_double =
    ballistics_engine::ffi::ballistics_reduce_qnh_pressure;

const FROZEN_DENSITY_ALTITUDE_TEMPERATURE_C: extern "C" fn(c_double, c_double) -> c_double =
    ballistics_engine::ffi::ballistics_density_altitude_temperature_c;

const FROZEN_DENSITY_ALTITUDE_PRESSURE_HPA: extern "C" fn(c_double, c_double) -> c_double =
    ballistics_engine::ffi::ballistics_density_altitude_pressure_hpa;

const FROZEN_DENSITY_ALTITUDE_ALTITUDE_M: extern "C" fn(c_double, c_double) -> c_double =
    ballistics_engine::ffi::ballistics_density_altitude_altitude_m;

const FROZEN_HOLD_POINT_IN_RETICLE: unsafe extern "C" fn(
    c_double,
    c_double,
    c_double,
    *const c_double,
    c_int,
    c_int,
    c_double,
    *mut FFIReticleHold,
) -> c_int = ballistics_engine::ffi::ballistics_hold_point_in_reticle;

const FROZEN_GET_VERSION: extern "C" fn() -> *const c_char =
    ballistics_engine::ffi::ballistics_get_version;

#[cfg(feature = "bridge")]
const FROZEN_BRIDGE_CALL: unsafe extern "C" fn(*const c_char) -> *mut c_char =
    ballistics_engine::bridge::ffi::ballistics_bridge_call;

#[cfg(feature = "bridge")]
const FROZEN_BRIDGE_CALL_N: unsafe extern "C" fn(*const u8, usize) -> *mut c_char =
    ballistics_engine::bridge::ffi::ballistics_bridge_call_n;

#[cfg(feature = "bridge")]
const FROZEN_BRIDGE_FREE: unsafe extern "C" fn(*mut c_char) =
    ballistics_engine::bridge::ffi::ballistics_bridge_free;

/// The constants above do the real work at compile time; this exists so they are *used*
/// (an unused `const` is not an error, but an unreferenced one invites deletion) and so the
/// count of frozen signatures is asserted against the frozen symbol list.
#[test]
fn every_frozen_export_has_a_frozen_signature() {
    // Addresses, purely to reference each constant. Cast through a thin pointer so the
    // differing function types can share one array.
    // `mut` goes unused without the `bridge` feature, which compiles out the extend below.
    #[cfg_attr(not(feature = "bridge"), allow(unused_mut))]
    let mut signatures: Vec<*const ()> = vec![
        FROZEN_CALCULATE_TRAJECTORY as *const (),
        FROZEN_CALCULATE_TRAJECTORY_WITH_DRAG_TABLE as *const (),
        FROZEN_CALCULATE_TRAJECTORY_WITH_DRAG_TABLE_SCALED as *const (),
        FROZEN_FREE_TRAJECTORY_RESULT as *const (),
        FROZEN_CALCULATE_ZERO_ANGLE as *const (),
        FROZEN_CALCULATE_ZERO_ANGLE_WITH_DRAG_TABLE as *const (),
        FROZEN_CALCULATE_ZERO_ANGLE_WITH_DRAG_TABLE_SCALED as *const (),
        FROZEN_QUICK_TRAJECTORY as *const (),
        FROZEN_MONTE_CARLO as *const (),
        FROZEN_MONTE_CARLO_WITH_DIRECTION_STD_DEV as *const (),
        FROZEN_FREE_MONTE_CARLO_RESULTS as *const (),
        FROZEN_BC_FOR_REFERENCE_STANDARD as *const (),
        FROZEN_REDUCE_QNH_PRESSURE as *const (),
        FROZEN_DENSITY_ALTITUDE_TEMPERATURE_C as *const (),
        FROZEN_DENSITY_ALTITUDE_PRESSURE_HPA as *const (),
        FROZEN_DENSITY_ALTITUDE_ALTITUDE_M as *const (),
        FROZEN_HOLD_POINT_IN_RETICLE as *const (),
        FROZEN_GET_VERSION as *const (),
    ];

    #[cfg(feature = "bridge")]
    signatures.extend([
        FROZEN_BRIDGE_CALL as *const (),
        FROZEN_BRIDGE_CALL_N as *const (),
        FROZEN_BRIDGE_FREE as *const (),
    ]);

    // Without `bridge`, its three symbols are not compiled and cannot be pinned.
    let expected = if cfg!(feature = "bridge") {
        FROZEN_EXPORTS.len()
    } else {
        FROZEN_EXPORTS.len() - 3
    };
    assert_eq!(
        signatures.len(),
        expected,
        "{} export(s) are in FROZEN_EXPORTS with no FROZEN_* signature constant, so a change \
         to their parameters or return type would not be caught",
        expected.abs_diff(signatures.len())
    );
    assert!(
        signatures.iter().all(|p| !p.is_null()),
        "an export resolved to a null address"
    );
}

// ---------------------------------------------------------------------------------------
// Layer 3: repr(C) layouts behind the pointers
// ---------------------------------------------------------------------------------------

/// Declare a mirror of a `repr(C)` struct and assert the live type matches it in size, align,
/// and every field offset.
///
/// Compared against a mirror rather than against hardcoded byte counts on purpose: the build
/// fleet includes 32-bit targets, where a struct holding pointers is legitimately smaller.
/// The mirror moves with the target; a hardcoded number would not.
macro_rules! frozen_repr_c {
    ($live:ty, $mirror:ident, { $($field:ident : $ty:ty),+ $(,)? }) => {
        #[repr(C)]
        #[allow(dead_code)]
        struct $mirror { $($field: $ty),+ }

        impl $mirror {
            fn assert_live_layout_unchanged() {
                assert_eq!(
                    size_of::<$live>(), size_of::<$mirror>(),
                    concat!(stringify!($live), " changed size: a field was added, removed, or \
                             retyped. Every caller compiled against the old layout now \
                             misreads it.")
                );
                assert_eq!(
                    align_of::<$live>(), align_of::<$mirror>(),
                    concat!(stringify!($live), " changed alignment")
                );
                $(
                    assert_eq!(
                        offset_of!($live, $field), offset_of!($mirror, $field),
                        concat!(stringify!($live), ".", stringify!($field),
                                " moved to a different byte offset -- fields were reordered or \
                                 one was inserted ahead of it. Appending to the END is the only \
                                 backward-compatible way to extend this struct.")
                    );
                )+
            }
        }
    };
}

frozen_repr_c!(FFIBallisticInputs, FrozenBallisticInputs, {
    muzzle_velocity: c_double,
    muzzle_angle: c_double,
    bc_value: c_double,
    bullet_mass: c_double,
    bullet_diameter: c_double,
    bc_type: c_int,
    sight_height: c_double,
    target_distance: c_double,
    temperature: c_double,
    twist_rate: c_double,
    is_twist_right: c_int,
    shooting_angle: c_double,
    altitude: c_double,
    latitude: c_double,
    azimuth_angle: c_double,
    use_rk4: c_int,
    use_adaptive_rk45: c_int,
    enable_wind_shear: c_int,
    enable_trajectory_sampling: c_int,
    sample_interval: c_double,
    enable_pitch_damping: c_int,
    enable_precession_nutation: c_int,
    enable_spin_drift: c_int,
    enable_magnus: c_int,
    enable_coriolis: c_int,
    shot_azimuth: c_double,
    cant_angle: c_double,
    zero_poi_vertical: c_double,
    zero_poi_horizontal: c_double,
    sight_offset_lateral: c_double,
});

frozen_repr_c!(FFIWindConditions, FrozenWindConditions, {
    speed: c_double,
    direction: c_double,
    vertical_speed: c_double,
});

frozen_repr_c!(FFIAtmosphericConditions, FrozenAtmosphericConditions, {
    temperature: c_double,
    pressure: c_double,
    humidity: c_double,
    altitude: c_double,
});

frozen_repr_c!(FFITrajectorySample, FrozenTrajectorySample, {
    distance: c_double,
    time: c_double,
    velocity_mps: c_double,
    energy_joules: c_double,
    drop_meters: c_double,
    windage_meters: c_double,
    mach: c_double,
    spin_rate_rps: c_double,
});

frozen_repr_c!(FFITrajectoryPoint, FrozenTrajectoryPoint, {
    time: c_double,
    position_x: c_double,
    position_y: c_double,
    position_z: c_double,
    velocity_magnitude: c_double,
    kinetic_energy: c_double,
});

frozen_repr_c!(FFITrajectoryResult, FrozenTrajectoryResult, {
    max_range: c_double,
    max_height: c_double,
    time_of_flight: c_double,
    impact_velocity: c_double,
    impact_energy: c_double,
    points: *mut FFITrajectoryPoint,
    point_count: c_int,
    sampled_points: *mut FFITrajectorySample,
    sampled_point_count: c_int,
    min_pitch_damping: c_double,
    transonic_mach: c_double,
    final_pitch_angle: c_double,
    final_yaw_angle: c_double,
    max_yaw_angle: c_double,
    max_precession_angle: c_double,
});

frozen_repr_c!(FFIMonteCarloParams, FrozenMonteCarloParams, {
    num_simulations: c_int,
    velocity_std_dev: c_double,
    angle_std_dev: c_double,
    bc_std_dev: c_double,
    wind_speed_std_dev: c_double,
    target_distance: c_double,
    base_wind_speed: c_double,
    base_wind_direction: c_double,
    azimuth_std_dev: c_double,
});

frozen_repr_c!(FFIMonteCarloResults, FrozenMonteCarloResults, {
    ranges: *mut c_double,
    impact_velocities: *mut c_double,
    impact_positions_x: *mut c_double,
    impact_positions_y: *mut c_double,
    impact_positions_z: *mut c_double,
    num_results: c_int,
    mean_range: c_double,
    std_dev_range: c_double,
    mean_impact_velocity: c_double,
    std_dev_impact_velocity: c_double,
    hit_probability: c_double,
});

frozen_repr_c!(FFIReticleHold, FrozenReticleHold, {
    down_mil: c_double,
    right_mil: c_double,
    nearest_mark: c_int,
    nearest_mark_distance_mil: c_double,
    off_reticle: c_int,
    mark_scale: c_double,
});

#[test]
fn repr_c_layouts_are_frozen() {
    FrozenBallisticInputs::assert_live_layout_unchanged();
    FrozenWindConditions::assert_live_layout_unchanged();
    FrozenAtmosphericConditions::assert_live_layout_unchanged();
    FrozenTrajectorySample::assert_live_layout_unchanged();
    FrozenTrajectoryPoint::assert_live_layout_unchanged();
    FrozenTrajectoryResult::assert_live_layout_unchanged();
    FrozenMonteCarloParams::assert_live_layout_unchanged();
    FrozenMonteCarloResults::assert_live_layout_unchanged();
    FrozenReticleHold::assert_live_layout_unchanged();
}
