/// Drag model enum that works in both Python and pure Rust contexts
#[cfg(feature = "python")]
use pyo3::prelude::*;

#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "python", pyclass(eq, eq_int))]
pub enum DragModel {
    G1,
    G2,
    G5,
    G6,
    G7,
    G8,
    GI,
    GS,
}

impl DragModel {
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_uppercase().as_str() {
            "G1" => Some(DragModel::G1),
            "G2" => Some(DragModel::G2),
            "G5" => Some(DragModel::G5),
            "G6" => Some(DragModel::G6),
            "G7" => Some(DragModel::G7),
            "G8" => Some(DragModel::G8),
            "GI" => Some(DragModel::GI),
            "GS" => Some(DragModel::GS),
            _ => None,
        }
    }
}

impl std::fmt::Display for DragModel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}