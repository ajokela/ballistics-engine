//! Shared display-ready row type behind the CLI's card-shaped table/CSV/JSON surfaces
//! (0.33.0 decision-support Plan B Task 9).
//!
//! `come-ups`, `range-table`, `wind-card` and `compare` each grew their own
//! function-local row struct (`ComeUpRow`/`RangeRow`/`WindRow`/`LoadRow`) that said the
//! same handful of things -- range, drop, wind, velocity, energy, time -- a different way,
//! which blocked any shared card machinery between them. `CardRow` is the display-ready
//! superset: every surface populates only the fields it has ever had and leaves the rest
//! `None` / empty, so each surface's existing table/CSV/JSON writer keeps reading the
//! identical numbers it always did (pinned byte-identical by `tests/card_golden_cli.rs`).
//!
//! No feature gate: this module must compile for `wasm32-unknown-unknown` with no default
//! features (pure data, no `fs`, no `clap`, no `pdf`). Task 10 rewrites the PDF dope card on
//! `&[CardRow]`; Task 11 grows an adaptive-card engine in this module.

/// One card row. Display-ready values in the surface's chosen units (exactly what the
/// legacy per-surface structs stored), so rendering is unchanged; range is f64 metres-
/// or-display per the surface's existing convention — DO NOT re-convert anything.
#[derive(Debug, Clone)]
pub struct CardRow {
    pub range: f64,
    pub drop_linear: Option<f64>,
    pub drop_adj: Option<f64>,
    pub come_up: Option<f64>,
    pub wind_linear: Option<f64>,
    pub wind_adj: Option<f64>,
    pub velocity: Option<f64>,
    pub energy: Option<f64>,
    pub time: Option<f64>,
    pub lead_adj: Option<f64>,
    /// wind-card's per-speed drift columns; empty elsewhere.
    pub wind_columns: Vec<f64>,
}
