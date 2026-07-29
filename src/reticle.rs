//! MBA-1361: reticle schema, parametric generators, and the hold-point-in-reticle API.
//!
//! Shooters who HOLD rather than dial need the answer expressed where they actually read
//! it: a point in their own reticle. This module is the engine slice of that — one shared
//! model (`serde`-serializable) that the CLI, the browser terminal, the FFI consumers and
//! the front ends can all speak, plus the coordinate transform from an angular firing
//! solution to a reticle coordinate. It is a transform plus a schema, not new physics: the
//! raw material (angular drop and drift in milliradians) is already a first-class output
//! everywhere in this crate.
//!
//! # Intellectual-property exclusions (deliberate, do not "fill in")
//!
//! Horus grid reticles and Time-of-Flight Wind Dots are actively patented, and Horus
//! monetizes app integration through its own licensed app. Therefore this module has, and
//! must keep having:
//!
//! * **no TREMOR-family / Horus grid layouts** — [`ReticleDescription::mil_grid`] builds a
//!   plain mil-hash CROSS (marks along the two stadia), never a filled two-dimensional
//!   grid, and [`ReticleDescription::tree`] is a generic parametric widening tree with no
//!   vendor geometry in it;
//! * **no wind-dot calibration** — nothing here maps a time of flight, a wind speed or a
//!   "wind hold number" onto a dot. Wind enters only as an angular deflection the caller
//!   already solved, in milliradians, exactly like elevation;
//! * **no vendor reticle catalog.** Manufacturer subtension sheets are published facts and
//!   are a legally viable catalog source, but curating one is a separate, per-vendor
//!   IP-reviewed data project (a tracked follow-up), not this module.
//!
//! # Angular conventions (the whole set, in one place)
//!
//! Every angle here is a **milliradian (mil)**, and every reticle coordinate is measured
//! **from the optical center**:
//!
//! * `down_mil` — POSITIVE is BELOW center. A holdover point is at positive `down_mil`.
//! * `right_mil` — POSITIVE is to the shooter's RIGHT of center.
//!
//! The hold point follows straight from that. If the bullet falls `d` mil below the line
//! of sight at some range, the shooter must place a reticle point `d` mil BELOW center on
//! the target — so `down_mil == drop_mil`. If the wind pushes the bullet `w` mil to the
//! RIGHT, the shooter must aim left by placing a point `w` mil to the RIGHT of center on
//! the target — so `right_mil == wind_mil`. [`hold_point_in_reticle`] therefore carries
//! the firing solution through unchanged and does the real work in the mark search; the
//! value of stating it here is that every surface now agrees on which way is which.
//!
//! # Focal plane
//!
//! Published optics-manual math, no more:
//!
//! * **FFP** (first focal plane): the reticle is magnified with the image, so a mark
//!   subtends the same angle at every magnification. Marks are used as authored.
//! * **SFP** (second focal plane): the reticle is a fixed angular size at the eyepiece
//!   while the target image scales, so a mark's TRUE subtension is
//!   `nominal * reference_magnification / magnification`. A 2 mil mark on a reticle
//!   calibrated at 10x covers 4 mil of target at 5x, and 1 mil at 20x.
//!
//! The hold point is a property of the trajectory, not of the optic, so it is always
//! TRUE angular. The mark search therefore scales the MARKS into true angular space and
//! compares there — never the other way round.

use std::error::Error;
use std::fmt;

use serde::{Deserialize, Serialize};

/// Largest mark count any generator will produce, and the most
/// [`ReticleDescription::validate`] will accept.
///
/// A reticle is a human-readable aiming device; nothing legitimate approaches this. The
/// cap exists so a hand-authored `--reticle-json` (or an FFI caller) cannot turn an
/// `O(marks)` search into an unbounded one, in the same spirit as the FFI drag-table
/// length guard (MBA-1407).
pub const MAX_RETICLE_MARKS: usize = 4096;

/// Fraction of each axis's mark span added as slack before a hold counts as
/// [`ReticleHold::off_reticle`]. See that field for the exact rule.
pub const OFF_RETICLE_MARGIN_FRACTION: f64 = 0.20;

/// Which focal plane the reticle is etched in.
///
/// Serialized as `"ffp"` / `"sfp"` — the spellings every optics catalog and every shooter
/// uses, rather than the Rust variant names.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum FocalPlane {
    /// First focal plane: subtensions are constant at every magnification.
    #[serde(rename = "ffp")]
    #[default]
    First,
    /// Second focal plane: subtensions scale as `reference_magnification / magnification`.
    #[serde(rename = "sfp")]
    Second,
}

impl FocalPlane {
    /// `"FFP"` / `"SFP"`, for tables and help text.
    pub fn label(self) -> &'static str {
        match self {
            FocalPlane::First => "FFP",
            FocalPlane::Second => "SFP",
        }
    }

    /// True when mark subtensions depend on magnification.
    pub fn is_magnification_dependent(self) -> bool {
        matches!(self, FocalPlane::Second)
    }
}

/// What a mark looks like. Purely descriptive — the hold math treats every kind
/// identically, and renderers use it to draw.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum MarkKind {
    /// A round dot.
    #[default]
    Dot,
    /// A short line across a stadium.
    Hash,
    /// A thick post (typically the lower stadium of a duplex).
    Post,
    /// The optical center / primary aiming point.
    Center,
}

impl MarkKind {
    /// Lower-case wire spelling, identical to the serde representation.
    pub fn as_str(self) -> &'static str {
        match self {
            MarkKind::Dot => "dot",
            MarkKind::Hash => "hash",
            MarkKind::Post => "post",
            MarkKind::Center => "center",
        }
    }
}

/// One aiming mark, positioned in NOMINAL angular units from the optical center.
///
/// "Nominal" means: as authored, i.e. the true subtension for an FFP reticle at any
/// magnification, and for an SFP reticle at its
/// [`ReticleDescription::reference_magnification`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReticleMark {
    /// Milliradians BELOW the optical center (negative = above).
    pub down_mil: f64,
    /// Milliradians RIGHT of the optical center (negative = left).
    pub right_mil: f64,
    /// How the mark is drawn. Does not affect the hold math.
    #[serde(default)]
    pub kind: MarkKind,
    /// Optional human label, e.g. a BDC mark's range (`"400 yd"`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
}

impl ReticleMark {
    /// A mark at `(down_mil, right_mil)` with no label.
    pub fn new(down_mil: f64, right_mil: f64, kind: MarkKind) -> Self {
        Self {
            down_mil,
            right_mil,
            kind,
            label: None,
        }
    }

    /// A mark at `(down_mil, right_mil)` carrying `label`.
    pub fn labeled(down_mil: f64, right_mil: f64, kind: MarkKind, label: impl Into<String>) -> Self {
        Self {
            down_mil,
            right_mil,
            kind,
            label: Some(label.into()),
        }
    }
}

/// A complete reticle: its focal plane, its calibration magnification, and its marks.
///
/// This is THE shared schema. It is deliberately permissive on unknown JSON keys (no
/// `deny_unknown_fields`) so a richer front-end description round-trips through the engine
/// without being rejected, and every field a solve needs is required.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReticleDescription {
    /// Display name, e.g. `"mil-grid 0.5/10"` or a vendor model designation.
    pub name: String,
    pub focal_plane: FocalPlane,
    /// The magnification at which [`ReticleMark`] subtensions are true. Meaningful only
    /// for [`FocalPlane::Second`]; ignored (and unvalidated) for FFP, where subtensions
    /// are magnification-independent by construction.
    pub reference_magnification: f64,
    pub marks: Vec<ReticleMark>,
}

/// Why a reticle operation was rejected. Typed rather than stringly, so front ends can
/// render their own wording and the FFI can map to a code.
#[derive(Debug, Clone, PartialEq)]
pub enum ReticleError {
    /// `magnification` was not finite and strictly positive. Checked on EVERY focal plane
    /// — an FFP hold does not depend on it, but zero magnification is not a physical
    /// optic and silently accepting it would mask a caller bug.
    NonPositiveMagnification { magnification: f64 },
    /// An SFP reticle carried a non-finite or non-positive
    /// [`ReticleDescription::reference_magnification`], which its subtension scaling
    /// divides by conceptually and multiplies by literally.
    NonPositiveReferenceMagnification { reference_magnification: f64 },
    /// The description carried no marks. A hold point has nothing to be near.
    NoMarks,
    /// The description carried more than [`MAX_RETICLE_MARKS`] marks.
    TooManyMarks { count: usize, max: usize },
    /// A mark's coordinates were not finite.
    NonFiniteMark { index: usize },
    /// The supplied firing solution (`drop_mil` / `wind_mil`) was not finite.
    NonFiniteHold { drop_mil: f64, wind_mil: f64 },
    /// A generator parameter violated its rule.
    InvalidGeneratorParameter {
        parameter: &'static str,
        value: f64,
        rule: &'static str,
    },
}

impl fmt::Display for ReticleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReticleError::NonPositiveMagnification { magnification } => write!(
                f,
                "magnification must be finite and greater than zero (got {magnification})"
            ),
            ReticleError::NonPositiveReferenceMagnification {
                reference_magnification,
            } => write!(
                f,
                "an SFP reticle's reference magnification must be finite and greater than \
                 zero (got {reference_magnification})"
            ),
            ReticleError::NoMarks => {
                write!(f, "the reticle description carries no marks")
            }
            ReticleError::TooManyMarks { count, max } => write!(
                f,
                "the reticle description carries {count} marks, more than the supported maximum of {max}"
            ),
            ReticleError::NonFiniteMark { index } => write!(
                f,
                "reticle mark {index} has non-finite coordinates"
            ),
            ReticleError::NonFiniteHold { drop_mil, wind_mil } => write!(
                f,
                "the hold must be finite (got drop {drop_mil} mil, wind {wind_mil} mil)"
            ),
            ReticleError::InvalidGeneratorParameter {
                parameter,
                value,
                rule,
            } => write!(f, "{parameter} must be {rule} (got {value})"),
        }
    }
}

impl Error for ReticleError {}

/// Where a firing solution lands in a reticle (MBA-1361).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReticleHold {
    /// TRUE angular milliradians BELOW center — equal to the supplied angular drop. See
    /// the module's conventions section for why this is an identity and not a transform.
    pub down_mil: f64,
    /// TRUE angular milliradians RIGHT of center — equal to the supplied angular wind
    /// deflection.
    pub right_mil: f64,
    /// Index into [`ReticleDescription::marks`] of the mark nearest the hold, in TRUE
    /// angular space (i.e. after SFP scaling). `None` only when the reticle has no marks,
    /// which [`hold_point_in_reticle`] rejects — so in practice this is always `Some`.
    pub nearest_mark: Option<usize>,
    /// Euclidean distance (mil) from the hold to that mark, measured in TRUE angular
    /// space. `0.0` when the hold lands exactly on a mark.
    pub nearest_mark_distance_mil: f64,
    /// True when the hold falls outside the marks' bounding box grown by
    /// [`OFF_RETICLE_MARGIN_FRACTION`] of that box's span, PER AXIS.
    ///
    /// Precisely: let `[lo, hi]` be the min/max of the TRUE-angular mark coordinates on an
    /// axis and `m = OFF_RETICLE_MARGIN_FRACTION * (hi - lo)`; the hold is off-reticle
    /// when it lies outside `[lo - m, hi + m]` on either axis. A degenerate axis (all
    /// marks share a coordinate, e.g. a pure BDC ladder with no windage marks) has
    /// `m == 0`, so ANY deviation on that axis reads as off-reticle. That is deliberate:
    /// such a reticle genuinely offers nothing to hold on in that direction.
    pub off_reticle: bool,
    /// The SFP subtension scale actually applied to the marks
    /// (`reference_magnification / magnification`); exactly `1.0` for FFP.
    pub mark_scale: f64,
}

/// The true-angular position of a mark after focal-plane scaling.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScaledMark {
    pub down_mil: f64,
    pub right_mil: f64,
}

impl ReticleDescription {
    /// Validate the description on its own terms: mark count, finiteness, and (SFP only)
    /// the reference magnification.
    pub fn validate(&self) -> Result<(), ReticleError> {
        if self.marks.is_empty() {
            return Err(ReticleError::NoMarks);
        }
        if self.marks.len() > MAX_RETICLE_MARKS {
            return Err(ReticleError::TooManyMarks {
                count: self.marks.len(),
                max: MAX_RETICLE_MARKS,
            });
        }
        for (index, mark) in self.marks.iter().enumerate() {
            if !mark.down_mil.is_finite() || !mark.right_mil.is_finite() {
                return Err(ReticleError::NonFiniteMark { index });
            }
        }
        if self.focal_plane.is_magnification_dependent()
            && (!self.reference_magnification.is_finite() || self.reference_magnification <= 0.0)
        {
            return Err(ReticleError::NonPositiveReferenceMagnification {
                reference_magnification: self.reference_magnification,
            });
        }
        Ok(())
    }

    /// The factor that converts NOMINAL mark subtensions to TRUE subtensions at
    /// `magnification`: `reference_magnification / magnification` for SFP, exactly `1.0`
    /// for FFP.
    ///
    /// Assumes [`Self::validate`] has passed and `magnification` is finite and positive.
    pub fn mark_scale(&self, magnification: f64) -> f64 {
        match self.focal_plane {
            FocalPlane::First => 1.0,
            FocalPlane::Second => self.reference_magnification / magnification,
        }
    }

    /// Every mark's TRUE angular position at `magnification`.
    pub fn scaled_marks(&self, magnification: f64) -> Result<Vec<ScaledMark>, ReticleError> {
        self.validate()?;
        require_positive_magnification(magnification)?;
        let scale = self.mark_scale(magnification);
        Ok(self
            .marks
            .iter()
            .map(|mark| ScaledMark {
                down_mil: mark.down_mil * scale,
                right_mil: mark.right_mil * scale,
            })
            .collect())
    }

    /// A plain mil-hash CROSS: marks every `spacing_mil` along the vertical and horizontal
    /// stadia out to `extent_mil`, plus a [`MarkKind::Center`] at the origin.
    ///
    /// This is NOT a filled two-dimensional grid — see the module header's IP exclusions.
    /// The result is FFP with a `reference_magnification` of 1.0 (unused for FFP);
    /// callers wanting an SFP grid set those two fields afterwards.
    pub fn mil_grid(spacing_mil: f64, extent_mil: f64) -> Result<Self, ReticleError> {
        require_generator_positive("spacing", spacing_mil)?;
        require_generator_positive("extent", extent_mil)?;
        if extent_mil < spacing_mil {
            return Err(ReticleError::InvalidGeneratorParameter {
                parameter: "extent",
                value: extent_mil,
                rule: "greater than or equal to the spacing",
            });
        }
        let steps = (extent_mil / spacing_mil).floor() as usize;
        // 1 center + 4 hashes per step; checked so a huge extent cannot wrap the count
        // below the cap and unleash the loop below.
        require_generated_size_checked(4usize.checked_mul(steps).and_then(|v| v.checked_add(1)))?;

        let mut marks = Vec::with_capacity(1 + 4 * steps);
        marks.push(ReticleMark::new(0.0, 0.0, MarkKind::Center));
        for step in 1..=steps {
            let offset = spacing_mil * step as f64;
            // Vertical stadium: below then above center.
            marks.push(ReticleMark::new(offset, 0.0, MarkKind::Hash));
            marks.push(ReticleMark::new(-offset, 0.0, MarkKind::Hash));
            // Horizontal stadium: right then left of center.
            marks.push(ReticleMark::new(0.0, offset, MarkKind::Hash));
            marks.push(ReticleMark::new(0.0, -offset, MarkKind::Hash));
        }
        Ok(Self {
            name: format!("mil-cross {spacing_mil}/{extent_mil}"),
            focal_plane: FocalPlane::First,
            reference_magnification: 1.0,
            marks,
        })
    }

    /// A generic parametric holdover tree: `rows` rows below center at `row_spacing_mil`
    /// intervals, each row `n` carrying windage dots at `±k * spread_step_mil` for
    /// `k` in `1..=n`, so the tree widens with depth. Plus a [`MarkKind::Center`].
    ///
    /// Generic geometry only — no vendor tree layout is reproduced here (module header).
    pub fn tree(
        rows: usize,
        row_spacing_mil: f64,
        spread_step_mil: f64,
    ) -> Result<Self, ReticleError> {
        if rows == 0 {
            return Err(ReticleError::InvalidGeneratorParameter {
                parameter: "rows",
                value: 0.0,
                rule: "at least 1",
            });
        }
        require_generator_positive("row-spacing", row_spacing_mil)?;
        require_generator_positive("spread-step", spread_step_mil)?;
        // 1 center + per row: the on-axis mark plus 2 windage dots per step. Checked so a
        // huge `rows` cannot wrap the product below the cap and unleash the loops below.
        require_generated_size_checked(
            rows.checked_add(1)
                .and_then(|r1| rows.checked_mul(r1))
                .and_then(|p| p.checked_add(rows))
                .and_then(|p| p.checked_add(1)),
        )?;

        let mut marks = Vec::new();
        marks.push(ReticleMark::new(0.0, 0.0, MarkKind::Center));
        for row in 1..=rows {
            let down = row_spacing_mil * row as f64;
            marks.push(ReticleMark::new(down, 0.0, MarkKind::Hash));
            for step in 1..=row {
                let spread = spread_step_mil * step as f64;
                marks.push(ReticleMark::new(down, spread, MarkKind::Dot));
                marks.push(ReticleMark::new(down, -spread, MarkKind::Dot));
            }
        }
        Ok(Self {
            name: format!("tree {rows}x{row_spacing_mil}/{spread_step_mil}"),
            focal_plane: FocalPlane::First,
            reference_magnification: 1.0,
            marks,
        })
    }

    /// A BDC ladder built from ALREADY-SOLVED drops: one labeled hash per
    /// `(range_m, drop_mil)` pair on the vertical stadium, plus a [`MarkKind::Center`].
    ///
    /// This generator deliberately does NOT run a solve. It is pure data assembly, so the
    /// caller stays in control of which load, atmosphere and zero the ladder describes,
    /// and this module keeps its "no physics" property.
    pub fn bdc_from_drops(drops: &[(f64, f64)]) -> Result<Self, ReticleError> {
        if drops.is_empty() {
            return Err(ReticleError::InvalidGeneratorParameter {
                parameter: "drops",
                value: 0.0,
                rule: "a non-empty list of (range, drop) pairs",
            });
        }
        require_generated_size(1 + drops.len())?;
        for &(range_m, drop_mil) in drops {
            if !range_m.is_finite() || range_m <= 0.0 {
                return Err(ReticleError::InvalidGeneratorParameter {
                    parameter: "drop range",
                    value: range_m,
                    rule: "finite and greater than zero",
                });
            }
            if !drop_mil.is_finite() {
                return Err(ReticleError::InvalidGeneratorParameter {
                    parameter: "drop",
                    value: drop_mil,
                    rule: "finite",
                });
            }
        }

        let mut marks = Vec::with_capacity(1 + drops.len());
        marks.push(ReticleMark::new(0.0, 0.0, MarkKind::Center));
        for &(range_m, drop_mil) in drops {
            marks.push(ReticleMark::labeled(
                drop_mil,
                0.0,
                MarkKind::Hash,
                format!("{range_m} m"),
            ));
        }
        Ok(Self {
            name: format!("bdc {} marks", drops.len()),
            focal_plane: FocalPlane::First,
            reference_magnification: 1.0,
            marks,
        })
    }
}

fn require_positive_magnification(magnification: f64) -> Result<(), ReticleError> {
    if !magnification.is_finite() || magnification <= 0.0 {
        return Err(ReticleError::NonPositiveMagnification { magnification });
    }
    Ok(())
}

fn require_generator_positive(parameter: &'static str, value: f64) -> Result<(), ReticleError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(ReticleError::InvalidGeneratorParameter {
            parameter,
            value,
            rule: "finite and greater than zero",
        });
    }
    Ok(())
}

fn require_generated_size(count: usize) -> Result<(), ReticleError> {
    if count > MAX_RETICLE_MARKS {
        return Err(ReticleError::TooManyMarks {
            count,
            max: MAX_RETICLE_MARKS,
        });
    }
    Ok(())
}

/// Same cap as [`require_generated_size`], but for a mark count assembled by `usize`
/// arithmetic that could overflow (a huge `--extent`/`--rows`). A `None` — the caller's
/// `checked_*` chain overflowed — is itself proof the count is far past the cap, so it is
/// rejected. Without this, `1 + 4 * steps` (or the tree product) wraps to a small value
/// and silently bypasses `MAX_RETICLE_MARKS`, then the generator loop runs unbounded.
fn require_generated_size_checked(count: Option<usize>) -> Result<(), ReticleError> {
    match count {
        Some(c) => require_generated_size(c),
        None => Err(ReticleError::TooManyMarks {
            count: usize::MAX,
            max: MAX_RETICLE_MARKS,
        }),
    }
}

/// Place a firing solution in a reticle (MBA-1361).
///
/// `drop_mil` is the angular drop below the line of sight (positive = below, i.e. the
/// come-up the shooter would otherwise dial) and `wind_mil` the angular wind deflection
/// (positive = the bullet goes RIGHT). `magnification` is the optic's CURRENT setting.
///
/// The returned hold coordinates are TRUE angular and equal to the inputs (see the module
/// header). The work this function does is the mark search: it scales the reticle's marks
/// into true angular space for the given magnification and focal plane, finds the nearest
/// one, and reports whether the hold has run off the marked part of the reticle.
///
/// # Errors
///
/// [`ReticleError::NonPositiveMagnification`] for a non-physical magnification (on EVERY
/// focal plane), [`ReticleError::NonPositiveReferenceMagnification`] for an SFP reticle
/// with no usable calibration magnification, [`ReticleError::NonFiniteHold`] for a
/// non-finite firing solution, and the [`ReticleDescription::validate`] errors for a
/// malformed description.
pub fn hold_point_in_reticle(
    drop_mil: f64,
    wind_mil: f64,
    magnification: f64,
    reticle: &ReticleDescription,
) -> Result<ReticleHold, ReticleError> {
    reticle.validate()?;
    require_positive_magnification(magnification)?;
    if !drop_mil.is_finite() || !wind_mil.is_finite() {
        return Err(ReticleError::NonFiniteHold { drop_mil, wind_mil });
    }

    let scale = reticle.mark_scale(magnification);
    let scaled: Vec<ScaledMark> = reticle
        .marks
        .iter()
        .map(|mark| ScaledMark {
            down_mil: mark.down_mil * scale,
            right_mil: mark.right_mil * scale,
        })
        .collect();

    // Nearest mark in TRUE angular space. Ties resolve to the LOWEST index (strict `<`),
    // which makes the answer independent of how the search is ordered.
    let mut nearest_index = 0usize;
    let mut nearest_distance = f64::INFINITY;
    for (index, mark) in scaled.iter().enumerate() {
        let distance = ((drop_mil - mark.down_mil).powi(2) + (wind_mil - mark.right_mil).powi(2))
            .sqrt();
        if distance < nearest_distance {
            nearest_distance = distance;
            nearest_index = index;
        }
    }

    let (down_lo, down_hi) = span(scaled.iter().map(|m| m.down_mil));
    let (right_lo, right_hi) = span(scaled.iter().map(|m| m.right_mil));
    let down_margin = OFF_RETICLE_MARGIN_FRACTION * (down_hi - down_lo);
    let right_margin = OFF_RETICLE_MARGIN_FRACTION * (right_hi - right_lo);
    let off_reticle = drop_mil < down_lo - down_margin
        || drop_mil > down_hi + down_margin
        || wind_mil < right_lo - right_margin
        || wind_mil > right_hi + right_margin;

    Ok(ReticleHold {
        down_mil: drop_mil,
        right_mil: wind_mil,
        nearest_mark: Some(nearest_index),
        nearest_mark_distance_mil: nearest_distance,
        off_reticle,
        mark_scale: scale,
    })
}

/// Min/max of a non-empty finite iterator. Callers have already validated finiteness.
fn span(values: impl Iterator<Item = f64>) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for value in values {
        if value < lo {
            lo = value;
        }
        if value > hi {
            hi = value;
        }
    }
    (lo, hi)
}

/// Rendering shape for the `reticle` command family.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReticleFormat {
    Table,
    Json,
}

/// Render a hold point, identically on every surface (MBA-1361).
///
/// This is THE formatter — the native CLI's `reticle hold` and the browser terminal's both
/// call it, so their output cannot drift apart. Same lesson as
/// [`crate::drag::format_reference_drag_curve`]: the `recoil` CSV header diverged between
/// the two surfaces precisely because each carried its own copy of the format strings.
///
/// Returned strings are newline-terminated; callers print or splice them verbatim.
pub fn format_reticle_hold(
    hold: &ReticleHold,
    reticle: &ReticleDescription,
    magnification: f64,
    format: ReticleFormat,
) -> String {
    let nearest = hold.nearest_mark.and_then(|index| reticle.marks.get(index));
    match format {
        ReticleFormat::Json => {
            let mut value = serde_json::json!({
                "reticle": reticle.name,
                "focal_plane": reticle.focal_plane.label(),
                "reference_magnification": reticle.reference_magnification,
                "magnification": magnification,
                "mark_scale": hold.mark_scale,
                "hold": {
                    "down_mil": hold.down_mil,
                    "right_mil": hold.right_mil,
                },
                "off_reticle": hold.off_reticle,
                "nearest_mark": serde_json::Value::Null,
            });
            if let (Some(index), Some(mark)) = (hold.nearest_mark, nearest) {
                let scale = hold.mark_scale;
                value["nearest_mark"] = serde_json::json!({
                    "index": index,
                    "kind": mark.kind.as_str(),
                    "label": mark.label,
                    "nominal_down_mil": mark.down_mil,
                    "nominal_right_mil": mark.right_mil,
                    "true_down_mil": mark.down_mil * scale,
                    "true_right_mil": mark.right_mil * scale,
                    "distance_mil": hold.nearest_mark_distance_mil,
                });
            }
            format!(
                "{}\n",
                serde_json::to_string_pretty(&value)
                    .unwrap_or_else(|_| "{\"error\":\"serialization failed\"}".to_string())
            )
        }
        ReticleFormat::Table => {
            let mut out = String::new();
            out.push_str("Reticle Hold Point\n");
            out.push_str("==================\n\n");
            out.push_str(&format!("Reticle:          {}\n", reticle.name));
            out.push_str(&format!(
                "Focal plane:      {}\n",
                reticle.focal_plane.label()
            ));
            if reticle.focal_plane.is_magnification_dependent() {
                out.push_str(&format!(
                    "Reference mag:    {:.2}x\n",
                    reticle.reference_magnification
                ));
                out.push_str(&format!("Magnification:    {magnification:.2}x\n"));
                out.push_str(&format!(
                    "Subtension scale: {:.4}x (marks read {:.4}x their etched value)\n",
                    hold.mark_scale, hold.mark_scale
                ));
            } else {
                out.push_str(&format!(
                    "Magnification:    {magnification:.2}x (FFP: subtensions are magnification-independent)\n"
                ));
            }
            out.push('\n');
            out.push_str(&format!("Hold down:  {:>8.3} mil\n", hold.down_mil));
            out.push_str(&format!("Hold right: {:>8.3} mil\n", hold.right_mil));
            out.push('\n');
            match nearest {
                Some(mark) => {
                    let scale = hold.mark_scale;
                    let label = mark.label.as_deref().unwrap_or("-");
                    out.push_str(&format!(
                        "Nearest mark:     #{} {} ({})\n",
                        hold.nearest_mark.unwrap_or(0),
                        mark.kind.as_str(),
                        label
                    ));
                    out.push_str(&format!(
                        "  at (down {:.3}, right {:.3}) mil true\n",
                        mark.down_mil * scale,
                        mark.right_mil * scale
                    ));
                    out.push_str(&format!(
                        "  distance from hold: {:.3} mil\n",
                        hold.nearest_mark_distance_mil
                    ));
                }
                None => out.push_str("Nearest mark:     none\n"),
            }
            if hold.off_reticle {
                out.push_str(
                    "\nWARNING: the hold falls outside the marked area of this reticle \
                     (dial instead, or use a reticle with more holdover).\n",
                );
            }
            out
        }
    }
}

/// Render a reticle description, identically on every surface (MBA-1361).
///
/// `-o json` emits the schema verbatim, so the output of `reticle generate ... -o json` is
/// exactly what `reticle hold --reticle-json` consumes.
pub fn format_reticle_description(reticle: &ReticleDescription, format: ReticleFormat) -> String {
    match format {
        ReticleFormat::Json => format!(
            "{}\n",
            serde_json::to_string_pretty(reticle)
                .unwrap_or_else(|_| "{\"error\":\"serialization failed\"}".to_string())
        ),
        ReticleFormat::Table => {
            let mut out = String::new();
            out.push_str(&format!("Reticle: {}\n", reticle.name));
            out.push_str(&format!(
                "Focal plane: {}",
                reticle.focal_plane.label()
            ));
            if reticle.focal_plane.is_magnification_dependent() {
                out.push_str(&format!(
                    "  Reference magnification: {:.2}x",
                    reticle.reference_magnification
                ));
            }
            out.push_str(&format!("  Marks: {}\n\n", reticle.marks.len()));
            out.push_str("   # Kind    Down(mil)  Right(mil)  Label\n");
            out.push_str("---- ------- ---------- ----------- --------------------\n");
            for (index, mark) in reticle.marks.iter().enumerate() {
                out.push_str(&format!(
                    "{:>4} {:<7} {:>10.3} {:>11.3}  {}\n",
                    index,
                    mark.kind.as_str(),
                    mark.down_mil,
                    mark.right_mil,
                    mark.label.as_deref().unwrap_or("-")
                ));
            }
            out
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sfp_two_mil_at_ten() -> ReticleDescription {
        ReticleDescription {
            name: "test sfp".to_string(),
            focal_plane: FocalPlane::Second,
            reference_magnification: 10.0,
            marks: vec![
                ReticleMark::new(0.0, 0.0, MarkKind::Center),
                ReticleMark::new(2.0, 0.0, MarkKind::Hash),
                ReticleMark::new(4.0, 0.0, MarkKind::Hash),
            ],
        }
    }

    fn ffp_ladder() -> ReticleDescription {
        ReticleDescription {
            name: "test ffp".to_string(),
            focal_plane: FocalPlane::First,
            reference_magnification: 1.0,
            marks: vec![
                ReticleMark::new(0.0, 0.0, MarkKind::Center),
                ReticleMark::new(2.0, 0.0, MarkKind::Hash),
                ReticleMark::new(4.0, 0.0, MarkKind::Hash),
                ReticleMark::new(2.0, 1.0, MarkKind::Dot),
                ReticleMark::new(2.0, -1.0, MarkKind::Dot),
            ],
        }
    }

    #[test]
    fn ffp_marks_are_invariant_across_magnification() {
        let reticle = ffp_ladder();
        let a = hold_point_in_reticle(2.3, 0.4, 4.0, &reticle).unwrap();
        let b = hold_point_in_reticle(2.3, 0.4, 25.0, &reticle).unwrap();
        assert_eq!(a, b, "FFP hold must not depend on magnification");
        assert_eq!(a.mark_scale, 1.0);
    }

    #[test]
    fn sfp_marks_scale_by_reference_over_current_magnification() {
        let reticle = sfp_two_mil_at_ten();

        // At the reference magnification the marks read their etched value.
        let scaled = reticle.scaled_marks(10.0).unwrap();
        assert_eq!(scaled[1].down_mil, 2.0);

        // Halving the magnification doubles what a mark covers: the 2 mil mark reads 4 mil.
        let scaled = reticle.scaled_marks(5.0).unwrap();
        assert_eq!(scaled[1].down_mil, 4.0);
        assert_eq!(scaled[2].down_mil, 8.0);

        // Doubling it halves them.
        let scaled = reticle.scaled_marks(20.0).unwrap();
        assert_eq!(scaled[1].down_mil, 1.0);
    }

    #[test]
    fn sfp_nearest_mark_is_measured_in_true_angular_space() {
        let reticle = sfp_two_mil_at_ten();
        // At 5x the etched 2 mil mark sits at 4 mil TRUE, so a 4 mil drop lands on it
        // exactly — the hold point itself is never rescaled.
        let hold = hold_point_in_reticle(4.0, 0.0, 5.0, &reticle).unwrap();
        assert_eq!(hold.nearest_mark, Some(1));
        assert_eq!(hold.nearest_mark_distance_mil, 0.0);
        assert_eq!(hold.down_mil, 4.0, "the hold stays TRUE angular");
        assert_eq!(hold.mark_scale, 2.0);

        // The same 4 mil drop at the reference magnification lands on the 4 mil mark.
        let hold = hold_point_in_reticle(4.0, 0.0, 10.0, &reticle).unwrap();
        assert_eq!(hold.nearest_mark, Some(2));
        assert_eq!(hold.nearest_mark_distance_mil, 0.0);
    }

    #[test]
    fn hold_on_a_mark_has_zero_distance() {
        let reticle = ffp_ladder();
        let hold = hold_point_in_reticle(2.0, 1.0, 10.0, &reticle).unwrap();
        assert_eq!(hold.nearest_mark, Some(3));
        assert_eq!(hold.nearest_mark_distance_mil, 0.0);
        assert!(!hold.off_reticle);
    }

    #[test]
    fn off_reticle_boundary_follows_the_documented_margin() {
        let reticle = ffp_ladder();
        // down span 0..4 => margin 0.8; right span -1..1 => margin 0.4.
        let inside = hold_point_in_reticle(4.8, 0.0, 10.0, &reticle).unwrap();
        assert!(!inside.off_reticle, "exactly on the margin counts as on-reticle");
        let outside = hold_point_in_reticle(4.80001, 0.0, 10.0, &reticle).unwrap();
        assert!(outside.off_reticle);

        let inside = hold_point_in_reticle(2.0, 1.4, 10.0, &reticle).unwrap();
        assert!(!inside.off_reticle);
        let outside = hold_point_in_reticle(2.0, 1.40001, 10.0, &reticle).unwrap();
        assert!(outside.off_reticle);

        // Above center is off the ladder too (span starts at 0.0).
        assert!(hold_point_in_reticle(-0.9, 0.0, 10.0, &reticle).unwrap().off_reticle);
    }

    #[test]
    fn sfp_off_reticle_uses_the_scaled_bounding_box() {
        let reticle = sfp_two_mil_at_ten();
        // At 5x the ladder reaches 8 mil TRUE (margin 1.6), so a 9 mil hold is still on.
        assert!(!hold_point_in_reticle(9.0, 0.0, 5.0, &reticle).unwrap().off_reticle);
        // At 20x it reaches only 2 mil TRUE (margin 0.4), so the same hold is far off.
        assert!(hold_point_in_reticle(9.0, 0.0, 20.0, &reticle).unwrap().off_reticle);
    }

    #[test]
    fn rejects_non_physical_magnification_on_both_planes() {
        for reticle in [ffp_ladder(), sfp_two_mil_at_ten()] {
            assert_eq!(
                hold_point_in_reticle(1.0, 0.0, 0.0, &reticle),
                Err(ReticleError::NonPositiveMagnification { magnification: 0.0 })
            );
            assert!(matches!(
                hold_point_in_reticle(1.0, 0.0, -3.0, &reticle),
                Err(ReticleError::NonPositiveMagnification { .. })
            ));
            assert!(matches!(
                hold_point_in_reticle(1.0, 0.0, f64::NAN, &reticle),
                Err(ReticleError::NonPositiveMagnification { .. })
            ));
        }
    }

    #[test]
    fn sfp_rejects_a_non_positive_reference_magnification_but_ffp_ignores_it() {
        let mut sfp = sfp_two_mil_at_ten();
        sfp.reference_magnification = 0.0;
        assert!(matches!(
            hold_point_in_reticle(1.0, 0.0, 10.0, &sfp),
            Err(ReticleError::NonPositiveReferenceMagnification { .. })
        ));

        let mut ffp = ffp_ladder();
        ffp.reference_magnification = 0.0;
        assert!(
            hold_point_in_reticle(1.0, 0.0, 10.0, &ffp).is_ok(),
            "FFP never consults the reference magnification"
        );
    }

    #[test]
    fn rejects_empty_and_non_finite_descriptions() {
        let mut reticle = ffp_ladder();
        reticle.marks.clear();
        assert_eq!(
            hold_point_in_reticle(1.0, 0.0, 10.0, &reticle),
            Err(ReticleError::NoMarks)
        );

        let mut reticle = ffp_ladder();
        reticle.marks[2].down_mil = f64::NAN;
        assert_eq!(
            hold_point_in_reticle(1.0, 0.0, 10.0, &reticle),
            Err(ReticleError::NonFiniteMark { index: 2 })
        );

        let reticle = ffp_ladder();
        assert!(matches!(
            hold_point_in_reticle(f64::NAN, 0.0, 10.0, &reticle),
            Err(ReticleError::NonFiniteHold { .. })
        ));
    }

    #[test]
    fn mil_grid_generates_a_cross_not_a_filled_grid() {
        let reticle = ReticleDescription::mil_grid(0.5, 2.0).unwrap();
        // 4 steps per arm plus the center.
        assert_eq!(reticle.marks.len(), 1 + 4 * 4);
        assert_eq!(reticle.marks[0].kind, MarkKind::Center);
        // Every non-center mark lies on exactly one axis — no (down, right) pair is both
        // non-zero, which is what distinguishes a cross from a grid.
        for mark in &reticle.marks[1..] {
            assert!(
                mark.down_mil == 0.0 || mark.right_mil == 0.0,
                "mil_grid must not produce off-axis marks"
            );
        }
        assert_eq!(reticle.focal_plane, FocalPlane::First);
    }

    #[test]
    fn tree_widens_one_step_per_row() {
        let reticle = ReticleDescription::tree(3, 1.0, 0.5).unwrap();
        // center + per row (1 on-axis + 2*row dots) = 1 + (1+2) + (1+4) + (1+6)
        assert_eq!(reticle.marks.len(), 1 + 3 + 2 * (1 + 2 + 3));
        let widest_in_row = |down: f64| {
            reticle
                .marks
                .iter()
                .filter(|m| m.down_mil == down)
                .map(|m| m.right_mil.abs())
                .fold(0.0_f64, f64::max)
        };
        assert_eq!(widest_in_row(1.0), 0.5);
        assert_eq!(widest_in_row(2.0), 1.0);
        assert_eq!(widest_in_row(3.0), 1.5);
    }

    #[test]
    fn bdc_from_drops_is_pure_data_assembly() {
        let reticle =
            ReticleDescription::bdc_from_drops(&[(300.0, 1.2), (400.0, 2.4), (500.0, 4.1)]).unwrap();
        assert_eq!(reticle.marks.len(), 4);
        assert_eq!(reticle.marks[0].kind, MarkKind::Center);
        assert_eq!(reticle.marks[1].down_mil, 1.2);
        assert_eq!(reticle.marks[1].label.as_deref(), Some("300 m"));
        assert_eq!(reticle.marks[3].down_mil, 4.1);
        // No windage component is invented.
        assert!(reticle.marks.iter().all(|m| m.right_mil == 0.0));
    }

    #[test]
    fn generators_reject_bad_parameters() {
        assert!(ReticleDescription::mil_grid(0.0, 5.0).is_err());
        assert!(ReticleDescription::mil_grid(0.5, 0.0).is_err());
        assert!(ReticleDescription::mil_grid(5.0, 1.0).is_err());
        assert!(ReticleDescription::tree(0, 1.0, 0.5).is_err());
        assert!(ReticleDescription::tree(3, -1.0, 0.5).is_err());
        assert!(ReticleDescription::bdc_from_drops(&[]).is_err());
        assert!(ReticleDescription::bdc_from_drops(&[(0.0, 1.0)]).is_err());
        assert!(ReticleDescription::bdc_from_drops(&[(100.0, f64::NAN)]).is_err());
        // The mark cap is enforced before a runaway grid is materialized.
        assert!(matches!(
            ReticleDescription::mil_grid(0.001, 10.0),
            Err(ReticleError::TooManyMarks { .. })
        ));
    }

    #[test]
    fn generator_sizes_cannot_overflow_past_the_cap() {
        // A step/row count large enough that the pre-cap `1 + 4*steps` (or the tree
        // product) WRAPS usize back below MAX_RETICLE_MARKS. Before the checked-size
        // guard these returned Ok and the generator loop ran unbounded (OOM/hang in
        // release, multiply-overflow panic in debug). Now they are rejected.
        let huge = (usize::MAX / 2) as f64; // exactly representable, cast without saturation
        assert!(matches!(
            ReticleDescription::mil_grid(1.0, huge),
            Err(ReticleError::TooManyMarks { .. })
        ));
        assert!(matches!(
            ReticleDescription::tree(usize::MAX / 2, 1.0, 0.5),
            Err(ReticleError::TooManyMarks { .. })
        ));
        // And the ordinary over-cap case (no overflow) still reports cleanly.
        assert!(matches!(
            ReticleDescription::tree(1000, 1.0, 0.5),
            Err(ReticleError::TooManyMarks { .. })
        ));
    }

    #[test]
    fn description_round_trips_through_serde() {
        let reticle = ReticleDescription {
            name: "round trip".to_string(),
            focal_plane: FocalPlane::Second,
            reference_magnification: 12.0,
            marks: vec![
                ReticleMark::new(0.0, 0.0, MarkKind::Center),
                ReticleMark::labeled(3.5, -1.5, MarkKind::Dot, "600 yd"),
                ReticleMark::new(6.0, 0.0, MarkKind::Post),
            ],
        };
        let json = serde_json::to_string(&reticle).unwrap();
        assert!(json.contains("\"sfp\""), "focal plane serializes as sfp: {json}");
        assert!(json.contains("\"center\""));
        let back: ReticleDescription = serde_json::from_str(&json).unwrap();
        assert_eq!(back, reticle);

        // An unlabeled mark emits no `label` key at all.
        let json = serde_json::to_string(&ReticleMark::new(1.0, 0.0, MarkKind::Hash)).unwrap();
        assert!(!json.contains("label"), "{json}");

        // Unknown keys are tolerated (front ends may carry render metadata).
        let permissive: ReticleDescription = serde_json::from_str(
            r#"{"name":"x","focal_plane":"ffp","reference_magnification":1.0,
                "marks":[{"down_mil":1.0,"right_mil":0.0,"kind":"hash"}],"stroke":"thin"}"#,
        )
        .unwrap();
        assert_eq!(permissive.marks.len(), 1);
    }

    #[test]
    fn hold_round_trips_through_serde() {
        let reticle = ffp_ladder();
        let hold = hold_point_in_reticle(3.1, 0.6, 10.0, &reticle).unwrap();
        let back: ReticleHold = serde_json::from_str(&serde_json::to_string(&hold).unwrap()).unwrap();
        assert_eq!(back, hold);
    }

    #[test]
    fn formatters_are_stable_and_json_is_the_schema() {
        let reticle = ReticleDescription::bdc_from_drops(&[(300.0, 1.2)]).unwrap();
        let json = format_reticle_description(&reticle, ReticleFormat::Json);
        assert!(json.ends_with('\n'));
        let back: ReticleDescription = serde_json::from_str(&json).unwrap();
        assert_eq!(back, reticle, "generate -o json feeds hold --reticle-json");

        let hold = hold_point_in_reticle(1.2, 0.0, 10.0, &reticle).unwrap();
        let table = format_reticle_hold(&hold, &reticle, 10.0, ReticleFormat::Table);
        assert!(table.contains("Hold down:"));
        assert!(table.contains("300 m"));
        assert!(!table.contains("Subtension scale"), "FFP hides SFP-only rows");

        let sfp = sfp_two_mil_at_ten();
        let hold = hold_point_in_reticle(4.0, 0.0, 5.0, &sfp).unwrap();
        let table = format_reticle_hold(&hold, &sfp, 5.0, ReticleFormat::Table);
        assert!(table.contains("Subtension scale: 2.0000x"));
    }
}
