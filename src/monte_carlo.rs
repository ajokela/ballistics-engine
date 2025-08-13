/// Monte Carlo simulation support
/// 
/// This module provides statistical analysis through Monte Carlo simulations.
/// The main implementation is in cli_api.rs for the CLI tool.

use crate::InternalBallisticInputs as BallisticInputs;

// Stub for compatibility
pub struct MonteCarloResults {
    pub ranges: Vec<f64>,
    pub impact_velocities: Vec<f64>,
}