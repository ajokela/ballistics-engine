//! Request-level perturbation kernel (Phase 1 of the 0.33.0 decision-support train).
pub mod taxonomy;
pub use taxonomy::{axes_in_group, axis_meta, AxisKind, AxisMeta, InputAxis, InputGroup};

// 0.33.0 decision-support Task 5: read/write access to a single taxonomy axis on a canonical
// request. No feature gate: must compile for wasm32 (depends only on solve_json/solve_v1,
// both unconditional).
pub mod access;
pub use access::{read_axis, with_axis, AxisValue, KernelError};
