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
pub use drag_model::DragModel;
pub use cli_api::{
    BallisticInputs, TrajectorySolver, WindConditions, AtmosphericConditions,
    TrajectoryResult, TrajectoryPoint, MonteCarloParams, MonteCarloResults,
    run_monte_carlo, run_monte_carlo_with_wind, calculate_zero_angle, calculate_zero_angle_with_conditions,
    estimate_bc_from_trajectory, BallisticsError,
};

// Module declarations
mod drag_model;
pub mod cli_api;
pub mod ffi;
#[cfg(target_arch = "wasm32")]
pub mod wasm;
#[cfg(test)]
mod wasm_tests;
// MBA-154: Make constants public for ballistics_rust wrapping
pub mod constants;
pub mod drag;
mod drag_tables;
pub mod atmosphere;
pub mod wind;
// MBA-153: Make wind_shear public for ballistics_rust wrapping
pub mod wind_shear;
// MBA-154: Make derivatives public for ballistics_rust wrapping
pub mod derivatives;
// MBA-154: Make trajectory_solver public for ballistics_rust wrapping
pub mod trajectory_solver;
pub mod trajectory_sampling;
// MBA-154: Make fast_trajectory public for ballistics_rust wrapping
pub mod fast_trajectory;
// MBA-149 Phase 5 Priority 2: Export enhanced spin_drift
pub mod spin_drift;
pub mod spin_drift_advanced;
pub mod spin_decay;
pub mod pitch_damping;
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
pub mod monte_carlo;
pub mod bc_estimation;
pub mod stability;
pub mod stability_advanced;
pub mod cluster_bc;

// Internal type alias for compatibility
pub(crate) type InternalBallisticInputs = BallisticInputs;

// BC segment data for velocity-dependent BC
#[derive(Debug, Clone)]
pub struct BCSegmentData {
    pub velocity_min: f64,
    pub velocity_max: f64,
    pub bc_value: f64,
}