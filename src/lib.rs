//! # Ballistics Engine
//! 
//! High-performance ballistics trajectory calculation engine with comprehensive physics modeling.

// Re-export the main types and functions
pub use drag_model::DragModel;
pub use cli_api::{
    BallisticInputs, TrajectorySolver, WindConditions, AtmosphericConditions,
    TrajectoryResult, TrajectoryPoint, MonteCarloParams, MonteCarloResults,
    run_monte_carlo, calculate_zero_angle, calculate_zero_angle_with_conditions,
    estimate_bc_from_trajectory, BallisticsError,
};

// Module declarations
mod drag_model;
pub mod cli_api;
pub mod ffi;
mod constants;
mod drag;
mod drag_tables;
mod atmosphere;
mod wind;
mod wind_shear;
mod derivatives;
mod trajectory_solver;
mod trajectory_sampling;
mod fast_trajectory;
mod spin_drift;
mod spin_decay;
mod pitch_damping;
mod precession_nutation;
mod aerodynamic_jump;
mod angle_calculations;
mod transonic_drag;
mod reynolds;
mod form_factor;
mod monte_carlo;
mod bc_estimation;
mod stability;

// Internal type alias for compatibility
pub(crate) type InternalBallisticInputs = BallisticInputs;

// BC segment data for velocity-dependent BC
#[derive(Debug, Clone)]
pub struct BCSegmentData {
    pub velocity_min: f64,
    pub velocity_max: f64,
    pub bc_value: f64,
}