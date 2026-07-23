//! # Ballistics Engine
//!
//! High-performance ballistics trajectory calculation engine with comprehensive physics modeling.
//!
//! ## Interactive Web Demo
//!
//! Try the ballistics engine directly in your browser at [https://ballistics.rs/](https://ballistics.rs/)
//!
//! ## Features
//!
//! - Professional-grade trajectory calculations with multiple drag models
//! - Advanced physics including spin drift, Coriolis effect, and Magnus force
//! - Monte Carlo simulations for uncertainty analysis
//! - WebAssembly support for browser-based applications
//! - FFI bindings for iOS and Android development

// Re-export the main types and functions
pub use cli_api::{
    calculate_zero_angle, calculate_zero_angle_with_conditions, estimate_bc_fit,
    estimate_bc_from_trajectory, interpolate_powder_temp_curve,
    resolve_powder_adjusted_velocity, run_monte_carlo, run_monte_carlo_with_direction_std_dev,
    run_monte_carlo_with_wind, run_monte_carlo_with_wind_and_direction_std_dev,
    run_monte_carlo_with_wind_and_direction_std_dev_seeded, AtmosphericConditions,
    BallisticInputs, BallisticsError, BcEstimate, BcFitMode, MonteCarloParams,
    MonteCarloResults, TrajectoryPoint, TrajectoryResult, TrajectorySolver, WindConditions,
    DEFAULT_HIT_RADIUS_M, MAX_TRAJECTORY_POINTS, TARGET_NOT_REACHED_SENTINEL_M,
};
pub use atmosphere::{AtmoSegment, AtmoSock};
pub use drag_model::DragModel;
pub use moving_target::{
    calculate_lead, lead_from_tof, mover_ring, LeadComponents, LeadError, LeadSolution,
};
pub use solve_json::{
    decode_solve_request_v1, ResolvedSolveRequestV1, SolveErrorCodeV1, SolveErrorEnvelopeV1,
    SolveRequestV1, SolveSuccessV1, MAX_SOLVE_JSON_SAMPLES_V1, SOLVE_JSON_SCHEMA_VERSION_V1,
};
pub use solve_v1::solve_v1;
pub use trajectory_observation::{
    TrajectoryObservation, TrajectoryObservationError, TrajectoryObservationFlag,
    TrajectoryTermination,
};
pub use trajectory_sampling::MAX_TRAJECTORY_SAMPLES;

// Module declarations
pub mod cli_api;
pub mod moving_target;
mod drag_model;
// The C ABI. Gated behind the default-on `ffi` feature so a binary that links two versions of
// this crate can disable it on one edge and avoid duplicate #[no_mangle] symbols.
#[cfg(feature = "ffi")]
pub mod ffi;
pub mod solve_json;
pub mod solve_v1;
pub mod terminal_plot;
// MBA-1343: multi-observation velocity/BC truing core, shared by the CLI and the WASM terminal.
pub mod truing;
// MBA-1346: observation-range experiment design for identifiable MV/BC truing.
pub mod truing_plan;
// MBA-1353: opt-in uncertainty-aware joint MV/BC truing.
pub mod truing_uncertainty;
// MBA-1343 Phase B: WEZ (`monte-carlo --wez`) sweep core, shared by the CLI and the WASM terminal.
pub mod wez;
// MBA-1355: turret adjustment-unit conversions (SMOA/IPHY/clicks) and click-value parsing,
// shared by the CLI and the WASM terminal. No feature gate: must compile for wasm32.
pub mod adjustment;
// MBA-1409: cleanroom `.drg` (Doppler drag-curve text file) parser, shared by the CLI and
// the WASM terminal. No feature gate: must compile for wasm32; parses TEXT only (no
// std::fs — file I/O stays in main.rs/wasm.rs).
pub mod drag_file;
pub mod trajectory_observation;
#[cfg(target_arch = "wasm32")]
pub mod wasm;
#[cfg(test)]
mod wasm_tests;
// MBA-154: Make constants public for ballistics_rust wrapping
pub mod atmosphere;
pub mod constants;
pub mod drag;
pub mod wind;
// MBA-153: Make wind_shear public for ballistics_rust wrapping
pub mod wind_shear;
// MBA-154: Make derivatives public for ballistics_rust wrapping
pub mod derivatives;
pub mod trajectory_sampling;
// MBA-154: Make fast_trajectory public for ballistics_rust wrapping
pub mod fast_trajectory;
// MBA-155: Add advanced integration methods (RK4, RK45)
pub mod trajectory_integration;
// MBA-149 Phase 5 Priority 2: Export enhanced spin_drift
pub mod pitch_damping;
pub mod spin_decay;
pub mod spin_drift;
pub mod spin_drift_advanced;
// MBA-149 Phase 5 Priority 2: Export enhanced precession_nutation
pub mod precession_nutation;
// MBA-153: Make aerodynamic_jump public for ballistics_rust wrapping
pub mod aerodynamic_jump;
// MBA-149 Phase 5 Priority 2: Export enhanced angle_calculations
pub mod angle_calculations;
pub mod transonic_drag;
// MBA-153: Make reynolds public for ballistics_rust wrapping
pub mod reynolds;
// MBA-149 Phase 5 Priority 2: Export enhanced form_factor
pub mod form_factor;
// MBA-153: Make monte_carlo public for ballistics_rust wrapping
pub mod bc_estimation;
pub mod cluster_bc;
pub mod monte_carlo;
pub mod stability;
pub mod stability_advanced;

// Online mode: HTTP client for Flask API (feature-gated)
#[cfg(feature = "online")]
pub mod api_client;

// BC5D table auto-download (feature-gated)
#[cfg(feature = "online")]
pub mod bc_table_download;

// BC correction table for offline mode
pub mod bc_table;

// 5D BC correction tables (caliber-specific, ML-derived)
pub mod bc_table_5d;

// Import of third-party ballistic profile files (.a7p), feature-gated
#[cfg(feature = "profile-import")]
pub mod profile_import;

// Internal type alias for compatibility
pub(crate) type InternalBallisticInputs = BallisticInputs;

// BC segment data for velocity-dependent BC
#[derive(Debug, Clone)]
pub struct BCSegmentData {
    pub velocity_min: f64,
    pub velocity_max: f64,
    pub bc_value: f64,
}
