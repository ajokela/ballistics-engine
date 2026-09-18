//! [`ProfileData`] -> ArcherBC2 `.a7p` encoder, and the honest account of what a
//! `.a7p` file has no place to put.
//!
//! The inverse of `profile_import::a7p`. File layout (32 ASCII hex characters of
//! MD5 over the remainder, then a proto3 `Payload { Profile profile = 1; }`
//! message), the field numbers, and the fixed-point scale factors are the same
//! interoperability facts the parser was written from; nothing is vendored from
//! the upstream a7p package (GPL-3.0 as distributed; this crate is MIT OR
//! Apache-2.0).
//!
//! # The not-carried list is the point
//!
//! `.a7p` is somebody else's format, built around one rifle, one load and one
//! device. [`ProfileData`] carries a good deal it has no slot for. Exporting
//! therefore returns [`A7pExport::not_carried`] alongside the bytes: EVERY
//! `ProfileData` field this encoder cannot put in the file, named exactly as the
//! struct and the saved-profile JSON spell it, each with whether this particular
//! profile actually held a value and what that value was.
//!
//! The list is unconditional — a field that cannot be carried appears whether or
//! not it was set — because a caller that only ever hears about populated fields
//! cannot distinguish "this profile had no scope tracking factor" from "the
//! tracking factor was silently thrown away". `populated` makes that distinction
//! explicit, and [`A7pExport::dropped_fields`] narrows the list to the ones that
//! really did lose data on this export.
//!
//! [`CARRIED_FIELDS`] names the other side of the same partition, and a test
//! asserts the two together account for every serialized `ProfileData` key. A
//! field added to `ProfileData` later cannot quietly join the drop set: the test
//! fails until somebody puts it in one list or the other.
//!
//! Solve-time effect flags (spin drift, Magnus, Coriolis, aerodynamic jump,
//! wind shear model) are deliberately absent from both lists: they are not
//! `ProfileData` fields at all, so they were never part of a profile to lose.
//!
//! One case sits deliberately in [`A7pExport::warnings`] rather than in the
//! not-carried list, and a caller showing a drop list should show the warnings
//! beside it. A `.a7p` coef row means (BC, m/s) for G1/G7 and (Cd, Mach) for
//! CUSTOM, never both at once, so a profile carrying a `drag_curve` under a G7
//! `drag_model` (or `bc_segments` under CUSTOM) loses the half its own
//! `drag_model` does not select. That is a contradictory profile rather than an
//! unsupported field — `drag_model` already decides which half is live, and
//! both names stay in [`CARRIED_FIELDS`] because either one CAN travel — so the
//! warning names the dropped half instead of the list claiming the format has no
//! room for it.
//!
//! # What "valid" means here
//!
//! The minimum switch count, the value ranges and the mandatory fields enforced
//! below are the upstream `a7p` package's VALIDATOR's rules, established
//! black-box against it; the ArcherBC2 device itself may well be stricter, so
//! nothing here should be read as the full set of what a device will accept.
//!
//! # `switches` is in neither list, on purpose
//!
//! The partition above is over `ProfileData` fields — things a shooter has that
//! may or may not reach the file. `.a7p`'s `switches` (device zoom/range presets)
//! runs the other way: a DESTINATION field with no source on our side. It is not
//! carried and it is not dropped, because there was never anything of the
//! shooter's in it, so putting it in either list would answer a question nobody
//! asked and crowd out the fields that really did lose data.
//!
//! It still cannot be skipped. The ecosystem refuses a file with fewer than
//! [`MIN_SWITCHES`] of them, so every export writes that many placeholders and
//! names them as placeholders in [`A7pExport::warnings`] — see
//! [`PLACEHOLDER_SWITCHES`] for why they are the upstream tool's own factory
//! values rather than anything computed from the profile.
//!
//! # Do not extend the format
//!
//! There is no escape hatch here for the dropped fields and there must not be
//! one. Smuggling our data into unused field numbers would produce files that
//! Archer's own tools misread, and would undo the reason the parser was written
//! cleanroom in the first place.

use serde::Serialize;

use super::wire::{write_bytes_field, write_i32_field, write_packed_i32_field, write_string_field};
use crate::cli_api::UnitSystem;
use crate::constants::{FPS_TO_MPS, GRAMS_PER_GRAIN};
use crate::profile::ProfileData;
use crate::profile_import::wrap_payload;

// Fixed-point scale factors (physical value * SCALE = the stored integer). The
// importer holds the identical table; they are duplicated rather than shared
// because the reader's copy is private to its own subtree, and the round-trip
// test at the bottom of this file is what keeps the two in lockstep — a scale
// edited on one side alone fails it immediately.
const SCALE_TWIST: f64 = 100.0; // inches/turn
const SCALE_VELOCITY: f64 = 10.0; // m/s
const SCALE_DIMENSION: f64 = 1000.0; // inches
const SCALE_WEIGHT: f64 = 10.0; // grains
const SCALE_PRESSURE: f64 = 10.0; // hPa
const SCALE_COEF: f64 = 10_000.0; // BC or Cd
const SCALE_DISTANCE: f64 = 100.0; // meters
const SCALE_MACH: f64 = 10_000.0; // Mach (CUSTOM coef rows only)

// Profile field numbers, matching the parser's `match field.number` arms.
const F_PROFILE_NAME: u32 = 1;
const F_BULLET_NAME: u32 = 3;
const F_SIGHT_HEIGHT: u32 = 9;
const F_TWIST: u32 = 10;
const F_MUZZLE_VELOCITY: u32 = 11;
const F_AIR_TEMPERATURE: u32 = 15;
const F_AIR_PRESSURE: u32 = 16;
const F_AIR_HUMIDITY: u32 = 17;
const F_BULLET_DIAMETER: u32 = 20;
const F_BULLET_WEIGHT: u32 = 21;
const F_BULLET_LENGTH: u32 = 22;
const F_TWIST_DIR: u32 = 23;
const F_BC_TYPE: u32 = 24;
const F_SWITCHES: u32 = 25;
const F_DISTANCES: u32 = 26;
const F_COEF_ROWS: u32 = 27;
// Sub-message field numbers inside one coef row.
const F_ROW_BC_CD: u32 = 1;
const F_ROW_MV: u32 = 2;
// Sub-message field numbers inside one switch position. Confirmed black-box:
// the upstream package's public factory API was asked for files with known
// switch values and the resulting BYTES were decoded here — no schema file was
// read and nothing from that package is vendored. Same footing as the field
// numbers and scale factors the importer documents.
const F_SW_C_IDX: u32 = 1;
const F_SW_ZOOM: u32 = 3;
const F_SW_DISTANCE: u32 = 4;
// Fields 2 (reticle_idx) and 5 (distance_from) are identified and deliberately
// unused: both placeholder values are the proto3 default, so they are omitted
// like every other default. Named here so a later reader does not have to
// re-derive them.
// Payload wrapper.
const F_PAYLOAD_PROFILE: u32 = 1;

/// `.a7p` REQUIRES at least four switch positions. This is not a stylistic
/// preference of the ecosystem's: the upstream validator refuses a file with
/// three or fewer outright (`'[] is too short'`), so canonical proto3
/// default-omission — correct for every other field here — produces a file
/// nothing in the Archer world will open.
///
/// Verified black-box by bisection against the upstream package: 0, 1, 2 and 3
/// switches are all refused; 4 and 5 are accepted.
///
/// DO NOT "clean up" the placeholders below back to an empty list. That silently
/// makes every file this crate writes invalid, and nothing in our own round trip
/// would notice, because our parser is happy either way.
const MIN_SWITCHES: usize = 4;

/// The four placeholder switch positions, as `(zoom, distance)` in the file's
/// own `SCALE_DISTANCE` fixed point — 100 m, 200 m, 300 m, 1000 m.
///
/// PLACEHOLDERS, NOT DATA. A switch position is a device UI preset: a zoom level
/// and the range the reticle is set up for. `ProfileData` has no such concept,
/// so there is nothing of the shooter's to put here and nothing is derived from
/// their profile. These are the values the upstream package's OWN factory
/// writes, chosen precisely so that what we emit is the ecosystem's neutral
/// default rather than a range card we invented and presented as theirs — a
/// fabricated preset at a range the shooter never shot would look exactly like
/// one they had. Every export says so in [`A7pExport::warnings`].
const PLACEHOLDER_SWITCHES: [(i32, i32); MIN_SWITCHES] =
    [(1, 10_000), (2, 20_000), (3, 30_000), (4, 100_000)];

/// The `c_idx` the upstream factory writes on every placeholder position: the
/// sentinel for "no distance-list index selected", which is what makes these
/// positions inert rather than pointing at a range-card entry.
const PLACEHOLDER_SWITCH_C_IDX: i32 = 255;

// Value ranges the ArcherBC2 ecosystem enforces, in the FILE's own integer
// units. Derived black-box by bisecting the upstream validator (its public API
// in, accept/refuse out); nothing from that package is vendored. These are
// REFUSALS here rather than warnings, because a file outside them is one the
// recipient cannot open at all — that is not a lossy export, it is no export,
// and reporting it as a success with a warning would be the worst of both.
const BOUND_SIGHT_HEIGHT: (i32, i32) = (-5_000, 5_000); // mm
const BOUND_TWIST: (i32, i32) = (0, 10_000); // in/turn x100
const BOUND_VELOCITY: (i32, i32) = (100, 30_000); // m/s x10
const BOUND_TEMPERATURE: (i32, i32) = (-100, 100); // C
const BOUND_PRESSURE: (i32, i32) = (3_000, 15_000); // hPa x10
const BOUND_HUMIDITY: (i32, i32) = (0, 100); // %
const BOUND_DIAMETER: (i32, i32) = (1, 50_000); // in x1000
const BOUND_WEIGHT: (i32, i32) = (10, 65_535); // gr x10
const BOUND_LENGTH: (i32, i32) = (10, 200_000); // in x1000
const BOUND_DISTANCE: (i32, i32) = (100, 300_000); // m x100
                                                   // Deliberately unbounded here: coef rows carry no per-value limits the validator
                                                   // enforces, and the 200-entry cap on the distances list cannot be reached by an
                                                   // exporter that writes exactly one distance.
const BOUND_NONE: (i32, i32) = (i32::MIN, i32::MAX);

// Unit factors. Deliberately a local copy of the handful of factors the CLI's
// `UnitConverter` (main.rs) uses for the same quantities: this module must
// compile for wasm32 and for lib-only builds, where main.rs does not exist. Any
// change to a factor here without the matching change there is a bug in this
// file, not a licence to diverge.
const MM_PER_INCH: f64 = 25.4;
const METERS_PER_YARD: f64 = 0.9144;
const HPA_PER_INHG: f64 = 33.8639;

/// One `ProfileData` field the `.a7p` format has no place for.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct NotCarried {
    /// The `ProfileData` field name, spelled exactly as the struct and the saved
    /// profile JSON spell it. Never a prose label: a caller has to be able to
    /// look it up in the document it was handed.
    pub field: &'static str,
    /// Whether THIS profile actually held something in that field. `false` means
    /// nothing was lost here; `true` means data went missing and the shooter's
    /// recipient will not receive it.
    pub populated: bool,
    /// What was in the field, when `populated` — short enough to show a user.
    pub value: Option<String>,
    /// Why `.a7p` cannot take it.
    pub reason: &'static str,
}

/// The `ProfileData` fields this encoder DOES put in the file. The complement of
/// the not-carried list; the two partition the struct, and
/// `not_carried_and_carried_cover_every_profile_field` proves it.
///
/// `units` is in here because it is consumed rather than dropped: it decides how
/// every other value is read on the way out, and the file itself is written in
/// the format's own fixed units. `auto_zero` is in here because it shares the
/// file's single zero distance with `zero_distance`; when the two disagree there
/// is only one slot, `zero_distance` wins, and the export warns.
pub const CARRIED_FIELDS: &[&str] = &[
    "name",
    "velocity",
    "bc",
    "mass",
    "diameter",
    "drag_model",
    "twist_rate",
    "sight_height",
    "zero_distance",
    "units",
    "temperature",
    "pressure",
    "humidity",
    "bullet_name",
    "auto_zero",
    "twist_right",
    "bullet_length",
    "bc_segments",
    "drag_curve",
];

/// A `.a7p` file plus everything the caller has to be told about it.
#[derive(Debug, Clone)]
pub struct A7pExport {
    /// The complete `.a7p` file: MD5 hex envelope followed by the proto3 payload.
    pub bytes: Vec<u8>,
    /// Every `ProfileData` field `.a7p` cannot carry — see the module doc.
    pub not_carried: Vec<NotCarried>,
    /// Losses and approximations INSIDE the fields that were carried: values
    /// rounded onto the format's fixed-point grid, destination fields left at
    /// their format default because `ProfileData` has no equivalent, and
    /// disagreements the encoder had to resolve.
    pub warnings: Vec<String>,
}

impl A7pExport {
    /// The not-carried fields that this profile actually had data in — the list
    /// to put in front of a shooter, as opposed to the full drop surface.
    pub fn dropped_fields(&self) -> Vec<&'static str> {
        self.not_carried
            .iter()
            .filter(|n| n.populated)
            .map(|n| n.field)
            .collect()
    }
}

/// Why an export could not happen at all. Distinct from the lossy-but-successful
/// outcomes above: these are profiles that cannot be expressed as a `.a7p` file
/// without inventing something.
#[derive(Debug, Clone, PartialEq)]
pub enum A7pExportError {
    /// The profile's `units` string is neither `imperial` nor `metric`.
    Units(String),
    /// A drag model `.a7p` has no enum value for. Refused rather than written as
    /// G1, which would silently hand the recipient different physics under a
    /// familiar-looking label.
    DragModel(String),
    /// A field whose value cannot be encoded: non-finite, non-positive where the
    /// format requires a real measurement, or past the range of the int32 the
    /// format stores it in.
    Field {
        field: &'static str,
        message: String,
    },
}

impl std::fmt::Display for A7pExportError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            A7pExportError::Units(u) => write!(f, "{u}"),
            A7pExportError::DragModel(m) => write!(
                f,
                "drag model '{m}' has no .a7p equivalent — the format stores only G1, G7 \
                 or a CUSTOM Mach/Cd curve, and writing '{m}' as one of those would give \
                 the recipient different physics under the same name"
            ),
            A7pExportError::Field { field, message } => {
                write!(f, "{field}: {message}")
            }
        }
    }
}

impl std::error::Error for A7pExportError {}

/// Accumulates the encoder's warnings while values are converted.
struct Lossy {
    warnings: Vec<String>,
}

impl Lossy {
    /// Convert one physical value to the stored integer, recording a warning
    /// when it does not land on the format's fixed-point grid. `unit` names the
    /// unit the value is in AT THIS POINT (the format's, not the profile's), so
    /// a warning reads in the numbers actually written to the file.
    fn quantize(
        &mut self,
        field: &'static str,
        value: f64,
        scale: f64,
        unit: &str,
        bound: (i32, i32),
    ) -> Result<i32, A7pExportError> {
        if !value.is_finite() {
            return Err(A7pExportError::Field {
                field,
                message: format!("{value} is not a finite number"),
            });
        }
        let raw = value * scale;
        if raw < f64::from(i32::MIN) || raw > f64::from(i32::MAX) {
            return Err(A7pExportError::Field {
                field,
                message: format!(
                    "{value} {unit} is outside the range .a7p can store (the file holds \
                     {unit} x {scale} in a 32-bit integer)"
                ),
            });
        }
        let rounded = raw.round();
        let stored = rounded as i32;
        // The ecosystem's own limits, checked on the STORED integer because that
        // is the number the recipient's validator sees.
        let (min, max) = bound;
        if stored < min || stored > max {
            return Err(A7pExportError::Field {
                field,
                message: format!(
                    "{value} {unit} is outside what .a7p accepts: the format stores \
                     {unit} x {scale} and the ecosystem's validator requires that \
                     integer to be between {min} and {max} (this would be {stored})"
                ),
            });
        }
        let recovered = f64::from(stored) / scale;
        // Relative tolerance: the conversions above (mm -> inches, grams ->
        // grains) do not produce exact binary values even when the underlying
        // measurement sits precisely on the grid, and warning about the last bit
        // of a double would make this list noise instead of information.
        if (recovered - value).abs() > 1e-9 * value.abs().max(1.0) {
            let grid = if scale == 1.0 {
                format!("whole {unit}")
            } else {
                format!("1/{scale} {unit}")
            };
            self.warnings.push(format!(
                "{field}: {value} {unit} does not fit .a7p's {grid} grid; \
                 written as {recovered} {unit}"
            ));
        }
        Ok(stored)
    }
}

fn require_positive(field: &'static str, value: f64) -> Result<f64, A7pExportError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(A7pExportError::Field {
            field,
            message: format!("{value} is not a usable measurement (must be finite and > 0)"),
        });
    }
    Ok(value)
}

/// `Some(value.to_string())` when the field holds something, `None` when it does
/// not — which is precisely [`NotCarried::populated`].
fn carried_opt<T: std::fmt::Display>(v: &Option<T>) -> Option<String> {
    v.as_ref().map(|v| v.to_string())
}

/// An optional list counts as populated only when it is non-empty: `Some(vec![])`
/// lost nothing, and reporting it as a drop would be the same false alarm the
/// unconditional list exists to avoid.
fn non_empty<T>(v: &Option<Vec<T>>) -> Option<&Vec<T>> {
    v.as_ref().filter(|v| !v.is_empty())
}

/// Encode a saved profile as an ArcherBC2 `.a7p` file.
///
/// Lossy by construction — see the module doc. The returned
/// [`A7pExport::not_carried`] is the contract; callers are expected to show it,
/// not to discard it.
pub fn export_a7p(profile: &ProfileData) -> Result<A7pExport, A7pExportError> {
    let units = profile.unit_system().map_err(A7pExportError::Units)?;
    let mut lossy = Lossy {
        warnings: Vec::new(),
    };

    // --- profile values, converted into the units .a7p stores ----------------
    // The format is fixed-unit (metric atmosphere, imperial bullet dimensions),
    // so the profile's own `units` is resolved here and never travels.
    let velocity_mps = match units {
        UnitSystem::Metric => profile.velocity,
        UnitSystem::Imperial => profile.velocity * FPS_TO_MPS,
    };
    let velocity_mps = require_positive("velocity", velocity_mps)?;

    let weight_grains = match units {
        UnitSystem::Metric => profile.mass / GRAMS_PER_GRAIN,
        UnitSystem::Imperial => profile.mass,
    };
    let weight_grains = require_positive("mass", weight_grains)?;

    let diameter_in = match units {
        UnitSystem::Metric => profile.diameter / MM_PER_INCH,
        UnitSystem::Imperial => profile.diameter,
    };
    let diameter_in = require_positive("diameter", diameter_in)?;

    // REQUIRED, not optional: the ecosystem's validator enforces a minimum
    // b_length, so a file written without one is refused outright. There is
    // nothing honest to substitute — a bullet length drives the recipient's own
    // stability and spin-drift model, and a plausible-looking invented one would
    // be indistinguishable from a measured one — so this refuses by name and
    // tells the caller to set the field.
    let length_in = match profile.bullet_length {
        Some(v) => match units {
            UnitSystem::Metric => v / MM_PER_INCH,
            UnitSystem::Imperial => v,
        },
        None => {
            return Err(A7pExportError::Field {
                field: "bullet_length",
                message: "the .a7p format requires a bullet length and this profile has \
                          none; set bullet_length (a length cannot be invented — it drives \
                          the recipient's stability model)"
                    .to_string(),
            })
        }
    };
    // `twist_rate` and `sight_height` are DELIBERATELY handled differently, and
    // the difference is decided rather than missed. Both are `Option` here, both
    // become 0 in the file when absent, and the validator accepts 0 for both — so
    // no rule derived from the validator alone separates them. What separates
    // them is what a 0 MEANS to the recipient:
    //
    //   * `sc_height` = 0 is a 0 mm mount, which no rifle has. It is an absence
    //     encoded as a number the recipient cannot tell from a measurement, and it
    //     moves every solution they compute. Refused below, like `bullet_length`.
    //   * `r_twist` = 0 may legitimately be this format's way of saying "twist
    //     unknown", in which case refusing would reject profiles the ecosystem
    //     handles perfectly well. We do not know that it is. We have the
    //     ecosystem's VALIDATOR, which accepts 0, and not the device, so this
    //     warns and writes nothing rather than refusing.
    //
    // What would settle it: observing how ArcherBC2 itself treats a profile with
    // r_twist = 0 — whether it presents the twist as unknown or silently solves
    // with a zero. If it is the latter, `twist_rate` should join `sight_height`
    // and `bullet_length` as a refusal. Until somebody can watch the device do it,
    // this asymmetry stands on a stated unknown, not on an oversight.
    let twist_in = profile.twist_rate.map(|v| match units {
        UnitSystem::Metric => v / MM_PER_INCH,
        UnitSystem::Imperial => v,
    });
    if twist_in.is_none() {
        lossy.warnings.push(
            "twist_rate: unset, so the file's r_twist stays 0. The format may read that as \
             \"twist unknown\" or as a zero twist, and which it is decides whether the \
             recipient gets spin drift at all — set twist_rate if the rifle's twist is \
             known"
                .to_string(),
        );
    }
    let sight_height_mm = match profile.sight_height {
        Some(v) => match units {
            UnitSystem::Metric => v,
            UnitSystem::Imperial => v * MM_PER_INCH,
        },
        None => {
            return Err(A7pExportError::Field {
                field: "sight_height",
                message: "the .a7p format has no way to say a sight height is unknown, and \
                          the 0 mm it would otherwise carry is a mount no rifle has — the \
                          recipient cannot tell it from a measurement and every solution \
                          they compute moves. Set sight_height"
                    .to_string(),
            })
        }
    };
    let zero_distance_m = profile.zero_distance.map(|v| match units {
        UnitSystem::Metric => v,
        UnitSystem::Imperial => v * METERS_PER_YARD,
    });
    let temperature_c = match units {
        UnitSystem::Metric => profile.temperature,
        UnitSystem::Imperial => (profile.temperature - 32.0) * 5.0 / 9.0,
    };
    let pressure_hpa = match units {
        UnitSystem::Metric => profile.pressure,
        UnitSystem::Imperial => profile.pressure * HPA_PER_INHG,
    };

    // --- drag: bc_type + coef rows ------------------------------------------
    // The one branch where the two sides of the format genuinely disagree about
    // what a coef row means: (BC, m/s) for G1/G7, (Cd, Mach) for CUSTOM.
    let model = profile.drag_model.trim().to_ascii_uppercase();
    let (bc_type, rows): (i32, Vec<(i32, i32)>) = match model.as_str() {
        "G1" | "G7" => {
            let bc_type = if model == "G1" { 0 } else { 1 };
            // `bc_field` follows the data: a warning about a rounded coefficient
            // has to name the field the number actually came from, or the reader
            // goes looking in the wrong place.
            let (segments, bc_field): (Vec<(f64, f64)>, &'static str) = match &profile.bc_segments {
                // `velocity_mps` is pinned to SI regardless of the profile's
                // `units` (see ProfileBcSegment), so it needs no conversion.
                Some(segments) if !segments.is_empty() => (
                    segments.iter().map(|s| (s.bc, s.velocity_mps)).collect(),
                    "bc_segments.bc",
                ),
                // No schedule: the scalar BC at the muzzle, which is exactly how
                // the importer reconstructs a scalar BC (the fastest row).
                _ => (
                    vec![(require_positive("bc", profile.bc)?, velocity_mps)],
                    "bc",
                ),
            };
            // The file has no separate scalar-BC slot: a reader rebuilds one from
            // the fastest coef row (that is exactly what the importer does). When
            // a schedule is present and the profile's own scalar `bc` disagrees
            // with its fastest row, the recipient will read the row — say so
            // rather than let the number change in transit unannounced.
            if profile.bc_segments.as_ref().is_some_and(|s| !s.is_empty()) {
                if let Some(fastest) = segments
                    .iter()
                    .max_by(|a, b| a.1.total_cmp(&b.1))
                    .map(|&(bc, _)| bc)
                {
                    if (profile.bc - fastest).abs() > 1e-9 {
                        lossy.warnings.push(format!(
                            "bc: the profile's scalar BC {} differs from the fastest \
                             bc_segments row {fastest}; .a7p has no separate scalar-BC \
                             field, so the recipient will read {fastest}",
                            profile.bc
                        ));
                    }
                }
            }
            let mut rows = Vec::with_capacity(segments.len());
            for (bc, mps) in segments {
                let velocity_field = if bc_field == "bc" {
                    "velocity"
                } else {
                    "bc_segments.velocity_mps"
                };
                rows.push((
                    lossy.quantize(bc_field, bc, SCALE_COEF, "BC", BOUND_NONE)?,
                    lossy.quantize(velocity_field, mps, SCALE_VELOCITY, "m/s", BOUND_NONE)?,
                ));
            }
            (bc_type, rows)
        }
        "CUSTOM" => {
            let curve = profile
                .drag_curve
                .as_ref()
                .filter(|c| !c.is_empty())
                .ok_or(A7pExportError::Field {
                    field: "drag_curve",
                    message: "drag_model is CUSTOM but the profile carries no drag curve"
                        .to_string(),
                })?;
            if profile.bc != 0.0 {
                lossy.warnings.push(format!(
                    "bc: {} is not carried for a CUSTOM drag model — the file's coef rows ARE \
                     the drag law, and a reader derives no scalar BC from them",
                    profile.bc
                ));
            }
            let mut rows = Vec::with_capacity(curve.len());
            for point in curve {
                rows.push((
                    lossy.quantize("drag_curve.cd", point.cd, SCALE_COEF, "Cd", BOUND_NONE)?,
                    lossy.quantize(
                        "drag_curve.mach",
                        point.mach,
                        SCALE_MACH,
                        "Mach",
                        BOUND_NONE,
                    )?,
                ));
            }
            (2, rows)
        }
        other => return Err(A7pExportError::DragModel(other.to_string())),
    };
    if profile.drag_curve.is_some() && bc_type != 2 {
        lossy.warnings.push(format!(
            "drag_curve: dropped — drag_model is {model}, and a .a7p profile stores EITHER \
             G1/G7 BC rows OR a CUSTOM Mach/Cd curve in the same field, never both"
        ));
    }
    if profile.bc_segments.is_some() && bc_type == 2 {
        lossy.warnings.push(
            "bc_segments: dropped — drag_model is CUSTOM, so the file's coef rows carry the \
             Mach/Cd curve instead"
                .to_string(),
        );
    }

    // --- zero distance -------------------------------------------------------
    // .a7p keeps a list of range-card distances and zeroes at one INDEX into it.
    // A ProfileData has one zero distance and no card, so the list is that one
    // distance at index 0 (the index is then the proto3 default and is omitted).
    // REQUIRED for the same reason as `bullet_length` above: the ecosystem's
    // validator refuses an empty distances list, and the one distance we have to
    // put in it is the zero distance. Guessing one would hand the recipient a
    // rifle zeroed somewhere its owner never zeroed it.
    let zero_distance_m = zero_distance_m.ok_or(A7pExportError::Field {
        field: "zero_distance",
        message: "the .a7p format requires at least one range-card distance and this \
                  profile has no zero distance to supply it; set zero_distance"
            .to_string(),
    })?;
    let distances: Vec<i32> = vec![lossy.quantize(
        "zero_distance",
        zero_distance_m,
        SCALE_DISTANCE,
        "m",
        BOUND_DISTANCE,
    )?];
    // `auto_zero` shares the file's single zero distance with `zero_distance`.
    // When they disagree there is only one slot, and `zero_distance` wins.
    if let (Some(auto), Some(zero)) = (profile.auto_zero, profile.zero_distance) {
        if (auto - zero).abs() > 1e-9 {
            lossy.warnings.push(format!(
                "auto_zero: {auto} differs from zero_distance {zero}; .a7p stores one zero \
                 distance, and zero_distance is what was written"
            ));
        }
    }

    // Destination fields with no ProfileData source. Named for the same reason
    // the not-carried list exists, in the other direction: the recipient's
    // device will show these, and they will be the format's defaults.
    lossy.warnings.push(
        "the file's cartridge_name, short_name_top, short_name_bot, user_note, caliber and \
         device_uuid are left empty, and its c_zero_temperature, c_zero_p_temperature, \
         c_t_coeff, c_zero_w_pitch and zero_x/zero_y stay at the format default: a saved \
         profile has no equivalent for any of them, and filling them in would be fabrication"
            .to_string(),
    );
    // `switches` is the one destination field that could NOT be left at its
    // default, because the ecosystem refuses a file with fewer than four. It is
    // therefore filled with placeholders, and the caller is told so in the same
    // breath — the alternative to saying it out loud is a shooter's friend seeing
    // four range presets that look like the shooter's own.
    lossy.warnings.push(format!(
        "the file's {MIN_SWITCHES} switch positions (device zoom/range presets at 100, 200, \
         300 and 1000 m) are PLACEHOLDERS, not this shooter's: a saved profile has no such \
         concept, and .a7p is refused by the ecosystem with fewer than {MIN_SWITCHES} of \
         them. They are the upstream tool's own factory values, so nothing here was derived \
         from the profile or invented as range data"
    ));

    // --- encode --------------------------------------------------------------
    // proto3 canonical form omits scalar fields equal to their default, which is
    // what other tools in this ecosystem are most likely to have been tested
    // against. The parser reads present-and-zero identically, so this is a
    // compatibility choice rather than a semantic one.
    let mut body = Vec::new();
    if !profile.name.is_empty() {
        write_string_field(F_PROFILE_NAME, &profile.name, &mut body);
    }
    if let Some(bullet_name) = profile.bullet_name.as_deref().filter(|s| !s.is_empty()) {
        write_string_field(F_BULLET_NAME, bullet_name, &mut body);
    }
    // Unconditional, unlike the optional fields above: a sight height is required
    // to get this far, so there is no "absent" case left to omit, and writing an
    // explicit 0 for a mount that rounds to nothing keeps this field's presence
    // in the file matching its presence in the profile. A present-and-zero
    // scalar reads identically to an omitted one, so this costs a byte and no
    // compatibility.
    write_i32_field(
        F_SIGHT_HEIGHT,
        lossy.quantize(
            "sight_height",
            sight_height_mm,
            1.0,
            "mm",
            BOUND_SIGHT_HEIGHT,
        )?,
        &mut body,
    );
    if let Some(inches) = twist_in {
        let raw = lossy.quantize("twist_rate", inches, SCALE_TWIST, "in/turn", BOUND_TWIST)?;
        if raw != 0 {
            write_i32_field(F_TWIST, raw, &mut body);
        }
    }
    write_i32_field(
        F_MUZZLE_VELOCITY,
        lossy.quantize(
            "velocity",
            velocity_mps,
            SCALE_VELOCITY,
            "m/s",
            BOUND_VELOCITY,
        )?,
        &mut body,
    );
    let temperature_raw =
        lossy.quantize("temperature", temperature_c, 1.0, "C", BOUND_TEMPERATURE)?;
    if temperature_raw != 0 {
        write_i32_field(F_AIR_TEMPERATURE, temperature_raw, &mut body);
    }
    write_i32_field(
        F_AIR_PRESSURE,
        lossy.quantize(
            "pressure",
            pressure_hpa,
            SCALE_PRESSURE,
            "hPa",
            BOUND_PRESSURE,
        )?,
        &mut body,
    );
    let humidity_raw = lossy.quantize("humidity", profile.humidity, 1.0, "%", BOUND_HUMIDITY)?;
    if humidity_raw != 0 {
        write_i32_field(F_AIR_HUMIDITY, humidity_raw, &mut body);
    }
    write_i32_field(
        F_BULLET_DIAMETER,
        lossy.quantize(
            "diameter",
            diameter_in,
            SCALE_DIMENSION,
            "in",
            BOUND_DIAMETER,
        )?,
        &mut body,
    );
    write_i32_field(
        F_BULLET_WEIGHT,
        lossy.quantize("mass", weight_grains, SCALE_WEIGHT, "gr", BOUND_WEIGHT)?,
        &mut body,
    );
    write_i32_field(
        F_BULLET_LENGTH,
        lossy.quantize(
            "bullet_length",
            length_in,
            SCALE_DIMENSION,
            "in",
            BOUND_LENGTH,
        )?,
        &mut body,
    );
    // TwistDir: RIGHT = 0 (the proto3 default, omitted), LEFT = 1. An unset
    // `twist_right` means the profile never recorded a direction; the format has
    // no way to say that, so it becomes the format's own default of RIGHT and
    // the not-carried list is not the right place to say so — the warning is.
    match profile.twist_right {
        Some(false) => write_i32_field(F_TWIST_DIR, 1, &mut body),
        Some(true) => {}
        None => lossy.warnings.push(
            "twist_right: unset; .a7p has no \"unknown\" twist direction, so the file says \
             RIGHT (its own default)"
                .to_string(),
        ),
    }
    if bc_type != 0 {
        write_i32_field(F_BC_TYPE, bc_type, &mut body);
    }
    // Switch positions (field 25) precede the distances list so the message stays
    // in ascending field order, which is what a canonical serializer emits.
    for (zoom, distance) in PLACEHOLDER_SWITCHES {
        let mut switch = Vec::new();
        write_i32_field(F_SW_C_IDX, PLACEHOLDER_SWITCH_C_IDX, &mut switch);
        // reticle_idx (field 2) and distance_from (field 5) are both 0 here — the
        // first reticle, and "distance is a value rather than an index into the
        // distances list" — so they are omitted like every other proto3 default.
        write_i32_field(F_SW_ZOOM, zoom, &mut switch);
        write_i32_field(F_SW_DISTANCE, distance, &mut switch);
        write_bytes_field(F_SWITCHES, &switch, &mut body);
    }
    if !distances.is_empty() {
        write_packed_i32_field(F_DISTANCES, &distances, &mut body);
        // c_zero_distance_idx (field 14) would be 0 here — the only entry — which
        // is the proto3 default, so it is omitted like every other default above.
        // A reader that sees no index resolves it to 0 as well.
    }
    for (bc_cd, mv) in &rows {
        let mut row = Vec::new();
        if *bc_cd != 0 {
            write_i32_field(F_ROW_BC_CD, *bc_cd, &mut row);
        }
        if *mv != 0 {
            write_i32_field(F_ROW_MV, *mv, &mut row);
        }
        write_bytes_field(F_COEF_ROWS, &row, &mut body);
    }

    let mut payload = Vec::new();
    write_bytes_field(F_PAYLOAD_PROFILE, &body, &mut payload);

    Ok(A7pExport {
        bytes: wrap_payload(&payload),
        not_carried: not_carried_for(profile),
        warnings: lossy.warnings,
    })
}

/// Build the unconditional not-carried list for one profile. Every entry names a
/// real `ProfileData` field; see the module doc for why absent fields are listed
/// too.
fn not_carried_for(p: &ProfileData) -> Vec<NotCarried> {
    let mut out = Vec::new();
    let mut push = |field: &'static str, value: Option<String>, reason: &'static str| {
        out.push(NotCarried {
            field,
            populated: value.is_some(),
            value,
            reason,
        });
    };

    push(
        // `altitude` is not an Option, so "populated" has to mean something: it
        // means non-zero. At sea level there is nothing to lose.
        "altitude",
        (p.altitude != 0.0).then(|| p.altitude.to_string()),
        "the file records the atmosphere at zeroing as temperature/pressure/humidity only; \
         it has no altitude field, and folding altitude into the pressure would change a \
         measured number into a derived one",
    );
    push(
        "density_altitude",
        carried_opt(&p.density_altitude),
        "no density-altitude concept in the format — see `altitude`",
    );
    push(
        "pressure_reference",
        carried_opt(&p.pressure_reference),
        "the file's pressure field declares no reference, so it is written as-is; a QNH \
         profile therefore travels as if it were station pressure",
    );
    push(
        "bc_reference",
        carried_opt(&p.bc_reference),
        "the file names no standard atmosphere for its BCs; a non-ICAO BC travels without \
         the note that says what it is referenced to",
    );
    push(
        "created",
        carried_opt(&p.created),
        "the format has no creation timestamp",
    );
    push(
        "wind_speed",
        carried_opt(&p.wind_speed),
        "a .a7p profile describes a rifle and a load, not a shot; it has no wind fields",
    );
    push(
        "wind_direction",
        carried_opt(&p.wind_direction),
        "no wind fields — see `wind_speed`",
    );
    push(
        "shooting_angle",
        carried_opt(&p.shooting_angle),
        "no look-angle field; a .a7p profile carries no shot conditions",
    );
    push(
        "use_bc_segments",
        carried_opt(&p.use_bc_segments),
        "no such switch in the format: a file's coef rows ARE its BC schedule, so there is \
         nothing to turn on or off",
    );
    push(
        "dsf_points",
        non_empty(&p.dsf_points).map(|v| format!("{} drop-scale-factor point(s)", v.len())),
        "truing results have no place in the format; the recipient gets the untrued load",
    );
    push(
        "elevation_cf",
        carried_opt(&p.elevation_cf),
        "no scope tracking-correction field; the recipient's solution will not be corrected \
         for how this scope actually tracks",
    );
    push(
        "windage_cf",
        carried_opt(&p.windage_cf),
        "no scope tracking-correction field — see `elevation_cf`",
    );
    push(
        "elevation_click",
        carried_opt(&p.elevation_click),
        "the file stores zeroing as raw device click COUNTS and never the click size itself; \
         a graduation written here would be read against the recipient device's own",
    );
    push(
        "windage_click",
        carried_opt(&p.windage_click),
        "no click-size field — see `elevation_click`",
    );
    push(
        "zero_poi_up_m",
        carried_opt(&p.zero_poi_up_m),
        "the file's only zero-offset slots are zero_x/zero_y, which are device-scoped click \
         counts with no click size recorded; converting a linear offset into clicks with THIS \
         profile's graduation would be misread by any device graduated differently",
    );
    push(
        "zero_poi_right_m",
        carried_opt(&p.zero_poi_right_m),
        "device-scoped click counts only — see `zero_poi_up_m`",
    );
    push(
        "sight_offset_lateral_m",
        carried_opt(&p.sight_offset_lateral_m),
        "the format models sight height only; there is no lateral sight-to-bore offset",
    );
    push(
        "zero_sets",
        non_empty(&p.zero_sets).map(|sets| {
            sets.iter()
                .map(|z| z.name.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        }),
        "one profile is one rifle with one zero in this format; alternate zeros and per-load \
         dial corrections have nowhere to go",
    );
    push(
        "reticle",
        p.reticle
            .as_ref()
            .map(|r| format!("\"{}\" ({} mark(s))", r.name, r.marks.len())),
        "the format describes the optic by sight height and click counts, not by its reticle",
    );
    push(
        "clicks_per_revolution",
        carried_opt(&p.clicks_per_revolution),
        "no turret-mechanics fields in the format",
    );
    push(
        "zero_stop",
        carried_opt(&p.zero_stop),
        "no turret-mechanics fields — see `clicks_per_revolution`",
    );
    push(
        "elevation_travel_up_mil",
        carried_opt(&p.elevation_travel_up_mil),
        "no turret travel limits in the format",
    );
    push(
        "elevation_travel_down_mil",
        carried_opt(&p.elevation_travel_down_mil),
        "no turret travel limits — see `elevation_travel_up_mil`",
    );
    push(
        "windage_travel_left_mil",
        carried_opt(&p.windage_travel_left_mil),
        "no turret travel limits — see `elevation_travel_up_mil`",
    );
    push(
        "windage_travel_right_mil",
        carried_opt(&p.windage_travel_right_mil),
        "no turret travel limits — see `elevation_travel_up_mil`",
    );
    push(
        "turret_elevation_dialed_mil",
        carried_opt(&p.turret_elevation_dialed_mil),
        "no field for what the turrets currently read",
    );
    push(
        "turret_windage_dialed_mil",
        carried_opt(&p.turret_windage_dialed_mil),
        "no field for what the turrets currently read",
    );
    push(
        "hold_bound_up_mil",
        carried_opt(&p.hold_bound_up_mil),
        "no reticle hold bounds in the format — see `reticle`",
    );
    push(
        "hold_bound_down_mil",
        carried_opt(&p.hold_bound_down_mil),
        "no reticle hold bounds — see `reticle`",
    );
    push(
        "hold_bound_left_mil",
        carried_opt(&p.hold_bound_left_mil),
        "no reticle hold bounds — see `reticle`",
    );
    push(
        "hold_bound_right_mil",
        carried_opt(&p.hold_bound_right_mil),
        "no reticle hold bounds — see `reticle`",
    );

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::profile_import::{map_a7p_to_profile, parse_a7p, EnvelopeStatus};

    /// Every `ProfileData` field populated, so the completeness test below sees
    /// the whole serialized key set (the `skip_serializing_if = Option::is_none`
    /// fields only appear when they hold something). Metric, and every carried
    /// measurement sits exactly on the `.a7p` fixed-point grid so the round-trip
    /// assertions can be about the mapping rather than about rounding.
    const FULL_FIXTURE: &str = r#"{
        "name": "export-fixture",
        "velocity": 792.0,
        "bc": 0.381,
        "mass": 19.439673,
        "diameter": 8.5852,
        "drag_model": "G7",
        "twist_rate": 254.0,
        "sight_height": 90.0,
        "zero_distance": 100.0,
        "units": "metric",
        "temperature": 15.0,
        "pressure": 1000.0,
        "humidity": 50.0,
        "altitude": 300.0,
        "bullet_name": "300GR OTM",
        "created": "1755400000",
        "wind_speed": 4.4704,
        "wind_direction": 90.0,
        "shooting_angle": -5.0,
        "auto_zero": 100.0,
        "twist_right": false,
        "use_bc_segments": true,
        "bullet_length": 45.72,
        "elevation_click": "0.1mil",
        "windage_click": "0.25moa",
        "bc_segments": [
            {"bc": 0.381, "velocity_mps": 792.0},
            {"bc": 0.360, "velocity_mps": 500.0}
        ],
        "drag_curve": [
            {"mach": 0.5, "cd": 0.23},
            {"mach": 3.0, "cd": 0.28}
        ],
        "dsf_points": [{"mach": 0.9, "dsf": 1.04}],
        "bc_reference": "army-standard-metro",
        "pressure_reference": "qnh",
        "density_altitude": 500.0,
        "zero_poi_up_m": 0.05,
        "zero_poi_right_m": -0.02,
        "sight_offset_lateral_m": 0.01,
        "elevation_cf": 0.97,
        "windage_cf": 1.02,
        "zero_sets": [{"name": "suppressed", "zero_distance": 200.0, "poi_up_mil": -0.3}],
        "reticle": {
            "name": "mil-grid 0.5/10",
            "focal_plane": "ffp",
            "reference_magnification": 10.0,
            "marks": [{"down_mil": 0.0, "right_mil": 0.0, "kind": "center"}]
        },
        "clicks_per_revolution": 100,
        "zero_stop": true,
        "elevation_travel_up_mil": 26.0,
        "elevation_travel_down_mil": 4.0,
        "windage_travel_left_mil": 12.0,
        "windage_travel_right_mil": 12.0,
        "turret_elevation_dialed_mil": 5.4,
        "turret_windage_dialed_mil": -0.2,
        "hold_bound_up_mil": 5.0,
        "hold_bound_down_mil": 10.0,
        "hold_bound_left_mil": 6.0,
        "hold_bound_right_mil": 6.0
    }"#;

    fn full_profile() -> ProfileData {
        serde_json::from_str(FULL_FIXTURE).expect("fixture must load")
    }

    /// The same rifle with nothing in it that `.a7p` cannot take — the profile
    /// the round-trip assertions are about.
    fn carryable_profile() -> ProfileData {
        let mut p = full_profile();
        p.altitude = 0.0;
        p.created = None;
        p.wind_speed = None;
        p.wind_direction = None;
        p.shooting_angle = None;
        p.use_bc_segments = None;
        p.elevation_click = None;
        p.windage_click = None;
        // G7 + a CUSTOM curve is a contradiction the file cannot hold both halves
        // of; the round trip is about the G7 half.
        p.drag_curve = None;
        p.dsf_points = None;
        p.bc_reference = None;
        p.pressure_reference = None;
        p.density_altitude = None;
        p.zero_poi_up_m = None;
        p.zero_poi_right_m = None;
        p.sight_offset_lateral_m = None;
        p.elevation_cf = None;
        p.windage_cf = None;
        p.zero_sets = None;
        p.reticle = None;
        p.clicks_per_revolution = None;
        p.zero_stop = None;
        p.elevation_travel_up_mil = None;
        p.elevation_travel_down_mil = None;
        p.windage_travel_left_mil = None;
        p.windage_travel_right_mil = None;
        p.turret_elevation_dialed_mil = None;
        p.turret_windage_dialed_mil = None;
        p.hold_bound_up_mil = None;
        p.hold_bound_down_mil = None;
        p.hold_bound_left_mil = None;
        p.hold_bound_right_mil = None;
        p
    }

    fn reimport(export: &A7pExport) -> ProfileData {
        let doc = parse_a7p(&export.bytes).expect("our own file must parse");
        assert!(
            matches!(doc.envelope, EnvelopeStatus::Verified),
            "the MD5 envelope we wrote must verify"
        );
        assert!(
            doc.unknown_fields.is_empty(),
            "the encoder must not emit field numbers the parser does not know: {:?}",
            doc.unknown_fields
        );
        map_a7p_to_profile(&doc, None, None)
            .expect("our own file must map")
            .profile
    }

    fn close(a: f64, b: f64, what: &str) {
        assert!(
            (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(1.0),
            "{what}: {a} != {b}"
        );
    }

    /// The load-bearing one: everything `CARRIED_FIELDS` claims survives really
    /// does survive a trip out through the encoder and back in through the
    /// EXISTING importer — which also pins the two modules' fixed-point scale
    /// tables against each other.
    #[test]
    fn round_trip_through_the_importer_preserves_every_carried_field() {
        let source = carryable_profile();
        let export = export_a7p(&source).expect("export");
        let back = reimport(&export);

        assert_eq!(back.name, source.name);
        assert_eq!(back.drag_model, source.drag_model);
        assert_eq!(back.units, source.units);
        assert_eq!(back.bullet_name, source.bullet_name);
        assert_eq!(back.twist_right, source.twist_right);
        close(back.velocity, source.velocity, "velocity");
        close(back.bc, source.bc, "bc");
        close(back.mass, source.mass, "mass");
        close(back.diameter, source.diameter, "diameter");
        close(back.temperature, source.temperature, "temperature");
        close(back.pressure, source.pressure, "pressure");
        close(back.humidity, source.humidity, "humidity");
        close(
            back.twist_rate.unwrap(),
            source.twist_rate.unwrap(),
            "twist",
        );
        close(
            back.sight_height.unwrap(),
            source.sight_height.unwrap(),
            "sight_height",
        );
        close(
            back.bullet_length.unwrap(),
            source.bullet_length.unwrap(),
            "bullet_length",
        );
        close(
            back.zero_distance.unwrap(),
            source.zero_distance.unwrap(),
            "zero_distance",
        );
        close(
            back.auto_zero.unwrap(),
            source.auto_zero.unwrap(),
            "auto_zero",
        );

        let (from, to) = (
            source.bc_segments.as_ref().unwrap(),
            back.bc_segments.as_ref().unwrap(),
        );
        assert_eq!(to.len(), from.len(), "bc_segments length");
        for (a, b) in from.iter().zip(to.iter()) {
            close(b.bc, a.bc, "bc_segments.bc");
            close(b.velocity_mps, a.velocity_mps, "bc_segments.velocity_mps");
        }
        assert!(back.drag_curve.is_none());

        // A clean profile must not generate rounding warnings. The only two are the
        // standing notes about destination fields with no source: the ones left at
        // the format default, and the placeholder switch positions.
        assert_eq!(
            export.warnings.len(),
            2,
            "unexpected warnings: {:?}",
            export.warnings
        );
    }

    /// `CARRIED_FIELDS` and the not-carried list partition `ProfileData`. A field
    /// added to the struct later belongs to one side or the other, and until
    /// somebody says which, this fails — which is the whole guarantee behind
    /// "names every field it cannot carry".
    #[test]
    fn not_carried_and_carried_cover_every_profile_field() {
        let profile = full_profile();
        let serialized = serde_json::to_value(&profile).expect("serialize");
        let keys: std::collections::BTreeSet<String> = serialized
            .as_object()
            .expect("a profile is a JSON object")
            .keys()
            .cloned()
            .collect();

        let export = export_a7p(&profile).expect("export");
        let mut accounted: std::collections::BTreeSet<String> =
            CARRIED_FIELDS.iter().map(|f| f.to_string()).collect();
        for entry in &export.not_carried {
            assert!(
                accounted.insert(entry.field.to_string()),
                "{} is in BOTH the carried and not-carried lists",
                entry.field
            );
        }

        let unaccounted: Vec<&String> = keys.difference(&accounted).collect();
        assert!(
            unaccounted.is_empty(),
            "ProfileData fields in neither list (add them to CARRIED_FIELDS or to \
             not_carried_for): {unaccounted:?}"
        );
        let invented: Vec<&String> = accounted.difference(&keys).collect();
        assert!(
            invented.is_empty(),
            "listed names that are not ProfileData fields: {invented:?}"
        );

        // `switches` is a DESTINATION field with no ProfileData source, so it
        // belongs to neither half of this partition — see the module doc. Pinned
        // because the obvious "fix" when someone meets the placeholder warning is
        // to add it to the not-carried list, where it would claim a shooter lost
        // something they never had.
        assert!(
            !accounted.contains("switches"),
            "switches is not a ProfileData field and must not appear in either list"
        );
    }

    /// The not-carried list is not decoration: for a profile that really does
    /// use the fields `.a7p` has no room for, they are reported as populated,
    /// carry the value that was lost, and are named exactly as the saved profile
    /// JSON names them.
    #[test]
    fn not_carried_names_real_populated_fields_with_their_values() {
        let profile = full_profile();
        let serialized = serde_json::to_value(&profile).expect("serialize");
        let object = serialized.as_object().expect("object");

        let export = export_a7p(&profile).expect("export");
        let dropped = export.dropped_fields();
        assert!(!dropped.is_empty(), "this profile drops plenty");

        for name in &dropped {
            assert!(
                object.contains_key(*name),
                "{name} is not a ProfileData field"
            );
        }
        // Spot-check across the different kinds of loss the ticket names:
        // truing results, effects-adjacent optic calibration, the DOPE-card-like
        // multi-zero list, and the per-device turret record.
        for expected in [
            "altitude",
            "dsf_points",
            "elevation_cf",
            "zero_sets",
            "reticle",
            "elevation_click",
            "zero_poi_up_m",
            "bc_reference",
            "clicks_per_revolution",
            "hold_bound_up_mil",
        ] {
            assert!(dropped.contains(&expected), "{expected} must be reported");
        }

        let dsf = export
            .not_carried
            .iter()
            .find(|n| n.field == "dsf_points")
            .expect("dsf_points entry");
        assert!(dsf.populated);
        assert_eq!(dsf.value.as_deref(), Some("1 drop-scale-factor point(s)"));
        assert!(!dsf.reason.is_empty());

        // Every entry is present whether or not it held anything, so a caller can
        // tell "had none" from "lost it".
        let mut empty = full_profile();
        empty.dsf_points = None;
        let empty_export = export_a7p(&empty).expect("export");
        let dsf = empty_export
            .not_carried
            .iter()
            .find(|n| n.field == "dsf_points")
            .expect("dsf_points is listed even when unset");
        assert!(!dsf.populated);
        assert_eq!(dsf.value, None);
    }

    /// The file is fixed-unit, so an imperial profile and its metric twin must
    /// produce the identical bytes. This is what proves the local unit factors
    /// are applied at all — a missing conversion would be invisible in a
    /// metric-only round trip.
    #[test]
    fn an_imperial_profile_and_its_metric_twin_encode_identically() {
        let metric = carryable_profile();
        let mut imperial = metric.clone();
        imperial.units = "imperial".to_string();
        imperial.velocity = 792.0 / FPS_TO_MPS; // fps
        imperial.mass = 300.0; // grains
        imperial.diameter = 0.338; // inches
        imperial.bullet_length = Some(1.8); // inches
        imperial.twist_rate = Some(10.0); // inches/turn
        imperial.sight_height = Some(90.0 / MM_PER_INCH); // inches
        imperial.zero_distance = Some(100.0 / METERS_PER_YARD); // yards
        imperial.auto_zero = imperial.zero_distance;
        imperial.temperature = 15.0 * 9.0 / 5.0 + 32.0; // F
        imperial.pressure = 1000.0 / HPA_PER_INHG; // inHg

        let a = export_a7p(&metric).expect("metric export");
        let b = export_a7p(&imperial).expect("imperial export");
        assert_eq!(
            a.bytes, b.bytes,
            "the same rifle must encode identically whatever units it was saved in"
        );
    }

    #[test]
    fn a_custom_drag_curve_round_trips_as_mach_cd_rows() {
        let mut profile = carryable_profile();
        profile.drag_model = "CUSTOM".to_string();
        profile.bc = 0.0; // the importer's inert sentinel
        profile.bc_segments = None;
        profile.drag_curve = Some(vec![
            crate::profile::ProfileDragPoint {
                mach: 0.5,
                cd: 0.23,
            },
            crate::profile::ProfileDragPoint {
                mach: 3.0,
                cd: 0.28,
            },
        ]);

        let export = export_a7p(&profile).expect("export");
        let back = reimport(&export);
        assert_eq!(back.drag_model, "CUSTOM");
        assert_eq!(back.bc, 0.0);
        let curve = back.drag_curve.as_ref().expect("curve survives");
        assert_eq!(curve.len(), 2);
        close(curve[0].mach, 0.5, "mach");
        close(curve[0].cd, 0.23, "cd");
        close(curve[1].mach, 3.0, "mach");
        close(curve[1].cd, 0.28, "cd");
    }

    /// The drag models `.a7p` has no enum value for are refused, not relabelled.
    /// Silently writing a G5 load as G1 would hand the recipient different
    /// physics under a name they would have no reason to doubt.
    #[test]
    fn a_drag_model_the_format_cannot_name_is_refused() {
        let mut profile = carryable_profile();
        profile.drag_model = "G5".to_string();
        match export_a7p(&profile) {
            Err(A7pExportError::DragModel(m)) => assert_eq!(m, "G5"),
            other => panic!("expected a DragModel refusal, got {other:?}"),
        }
    }

    /// Losses INSIDE a carried field are warnings, never silence. `.a7p` stores
    /// sight height as whole millimetres, so 50.8 mm (a 2-inch mount) cannot
    /// survive and the caller has to hear about it.
    #[test]
    fn a_value_that_does_not_fit_the_fixed_point_grid_is_warned_about() {
        let mut profile = carryable_profile();
        profile.sight_height = Some(50.8);
        let export = export_a7p(&profile).expect("export");
        assert!(
            export
                .warnings
                .iter()
                .any(|w| w.starts_with("sight_height:") && w.contains("51")),
            "expected a sight_height quantization warning, got {:?}",
            export.warnings
        );
        let back = reimport(&export);
        assert_eq!(back.sight_height, Some(51.0));
    }

    #[test]
    fn an_unusable_measurement_is_a_named_field_error() {
        let mut profile = carryable_profile();
        profile.velocity = 0.0;
        match export_a7p(&profile) {
            Err(A7pExportError::Field { field, .. }) => assert_eq!(field, "velocity"),
            other => panic!("expected a Field error, got {other:?}"),
        }

        let mut profile = carryable_profile();
        profile.units = "furlongs".to_string();
        assert!(matches!(
            export_a7p(&profile),
            Err(A7pExportError::Units(_))
        ));
    }

    /// A CUSTOM profile with no curve has nothing to encode, and inventing a G1
    /// BC for it would be exactly the silent fabrication the importer's own
    /// CUSTOM handling refuses.
    #[test]
    fn custom_without_a_curve_is_refused() {
        let mut profile = carryable_profile();
        profile.drag_model = "CUSTOM".to_string();
        profile.drag_curve = None;
        match export_a7p(&profile) {
            Err(A7pExportError::Field { field, .. }) => assert_eq!(field, "drag_curve"),
            other => panic!("expected a drag_curve error, got {other:?}"),
        }
    }

    /// Count the switch positions (field 25) in an exported file, using a reader
    /// written here rather than the importer's — the importer only COUNTS
    /// switches, and this test has to be able to see inside one.
    fn switch_entries(bytes: &[u8]) -> Vec<Vec<(u32, u64)>> {
        fn varint(b: &[u8], i: &mut usize) -> u64 {
            let (mut v, mut shift) = (0u64, 0u32);
            loop {
                let byte = b[*i];
                *i += 1;
                v |= u64::from(byte & 0x7f) << shift;
                if byte & 0x80 == 0 {
                    return v;
                }
                shift += 7;
            }
        }
        fn fields(b: &[u8]) -> Vec<(u32, u64, &[u8])> {
            let mut i = 0usize;
            let mut out = Vec::new();
            while i < b.len() {
                let key = varint(b, &mut i);
                let (number, wire) = ((key >> 3) as u32, key & 7);
                match wire {
                    0 => {
                        let v = varint(b, &mut i);
                        out.push((number, v, &b[0..0]));
                    }
                    2 => {
                        let n = varint(b, &mut i) as usize;
                        out.push((number, 0, &b[i..i + n]));
                        i += n;
                    }
                    other => panic!("unexpected wire type {other}"),
                }
            }
            out
        }
        let payload = &bytes[32..];
        let profile = fields(payload)
            .into_iter()
            .find(|(n, _, _)| *n == F_PAYLOAD_PROFILE)
            .expect("payload carries a profile")
            .2;
        fields(profile)
            .into_iter()
            .filter(|(n, _, _)| *n == F_SWITCHES)
            .map(|(_, _, body)| {
                fields(body)
                    .into_iter()
                    .map(|(n, v, _)| (n, v))
                    .collect::<Vec<_>>()
            })
            .collect()
    }

    /// THE REGRESSION GUARD. `.a7p` is refused outright by the ecosystem's own
    /// validator when it carries fewer than four switch positions — verified
    /// black-box against the upstream package by bisection (0/1/2/3 refused,
    /// 4/5 accepted). Our own parser is happy either way, so nothing in the round
    /// trip above would catch it: without this test, "cleaning up" the
    /// placeholders back to canonical proto3 default-omission would silently make
    /// every file this crate writes unopenable, and every test would still pass.
    #[test]
    fn every_export_carries_the_minimum_four_switch_positions() {
        for profile in [carryable_profile(), full_profile()] {
            let export = export_a7p(&profile).expect("export");
            let switches = switch_entries(&export.bytes);
            assert!(
                switches.len() >= MIN_SWITCHES,
                "exported {} switch positions; the ecosystem refuses fewer than \
                 {MIN_SWITCHES}",
                switches.len()
            );
            for (i, entry) in switches.iter().enumerate() {
                let (zoom, distance) = PLACEHOLDER_SWITCHES[i];
                assert_eq!(
                    entry,
                    &vec![
                        (F_SW_C_IDX, PLACEHOLDER_SWITCH_C_IDX as u64),
                        (F_SW_ZOOM, zoom as u64),
                        (F_SW_DISTANCE, distance as u64),
                    ],
                    "switch {i}"
                );
            }
            // And they are declared as placeholders rather than passed off as the
            // shooter's own presets.
            assert!(
                export
                    .warnings
                    .iter()
                    .any(|w| w.contains("PLACEHOLDERS") && w.contains("switch")),
                "the placeholder switches must be declared: {:?}",
                export.warnings
            );
        }
    }

    /// The three fields `ProfileData` leaves optional that this export refuses to
    /// write without. Refused by name rather than filled in: each would otherwise
    /// reach the recipient as a number they cannot tell from a measurement — a
    /// bullet length drives their stability model, a zero distance is where their
    /// rifle will shoot, and a 0 mm mount is one no rifle has.
    #[test]
    fn fields_the_format_requires_are_refused_when_absent_not_invented() {
        for (field, clear) in [
            (
                "bullet_length",
                Box::new(|p: &mut ProfileData| p.bullet_length = None)
                    as Box<dyn Fn(&mut ProfileData)>,
            ),
            (
                "zero_distance",
                Box::new(|p: &mut ProfileData| p.zero_distance = None),
            ),
            (
                "sight_height",
                Box::new(|p: &mut ProfileData| p.sight_height = None),
            ),
        ] {
            let mut profile = carryable_profile();
            clear(&mut profile);
            match export_a7p(&profile) {
                Err(A7pExportError::Field { field: got, .. }) => assert_eq!(got, field),
                other => panic!("{field}: expected a refusal, got {other:?}"),
            }
        }
    }

    /// `twist_rate` is the DELIBERATE exception to the rule above, and this test
    /// exists to stop the three refusals and this one warning being tidied into a
    /// single shared rule.
    ///
    /// An absent twist still exports, with a warning. The reason is not that a
    /// zero twist matters less than a zero mount — it is that `r_twist` = 0 may be
    /// how this format says "twist unknown", which would make a refusal reject
    /// profiles the ecosystem handles fine. We have the ecosystem's validator,
    /// which accepts 0, and not the device that interprets it. See the comment at
    /// the conversion site for what observation would settle it and let this
    /// become a refusal.
    #[test]
    fn an_absent_twist_rate_warns_rather_than_refusing_and_that_is_deliberate() {
        let mut profile = carryable_profile();
        profile.twist_rate = None;
        let export = export_a7p(&profile).expect("an unknown twist must still export");
        assert!(
            export
                .warnings
                .iter()
                .any(|w| w.starts_with("twist_rate:") && w.contains("r_twist")),
            "an absent twist must be reported: {:?}",
            export.warnings
        );
        // ...and the file really does carry no twist, rather than a fabricated one.
        let back = reimport(&export);
        assert_eq!(back.twist_rate, Some(0.0));
    }

    /// Values our own arithmetic accepts but the ecosystem's validator does not
    /// are refused here rather than written into a file the recipient cannot
    /// open. Ranges derived black-box from the upstream validator.
    #[test]
    fn values_outside_what_the_ecosystem_accepts_are_refused_by_name() {
        // Each of these is a number our own arithmetic encodes without complaint —
        // finite, positive, well inside an i32 — and each lands outside a limit the
        // recipient's validator enforces.
        for (field, mutate) in [
            (
                "velocity",
                Box::new(|p: &mut ProfileData| p.velocity = 5.0) as Box<dyn Fn(&mut ProfileData)>,
            ),
            (
                "pressure",
                Box::new(|p: &mut ProfileData| p.pressure = 100.0),
            ),
            (
                "temperature",
                Box::new(|p: &mut ProfileData| p.temperature = 250.0),
            ),
            (
                "zero_distance",
                Box::new(|p: &mut ProfileData| p.zero_distance = Some(5000.0)),
            ),
            (
                "sight_height",
                Box::new(|p: &mut ProfileData| p.sight_height = Some(9000.0)),
            ),
        ] {
            let mut profile = carryable_profile();
            mutate(&mut profile);
            match export_a7p(&profile) {
                Err(A7pExportError::Field { field: got, .. }) => assert_eq!(got, field),
                other => panic!("{field}: expected a refusal, got {other:?}"),
            }
        }
    }
}
