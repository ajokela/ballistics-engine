//! MBA-1440: import Bero's "Ventum" reticle spec into the engine's [`ReticleDescription`].
//!
//! A reticle drawn in Bero's Ventum tool is a single JSON object: reticle-level metadata
//! (name, focal plane, calibration magnification, angular unit) plus a `spec` list of
//! drawing elements. This module is the one-way transform from that authoring format into
//! the engine's own [`ReticleDescription`], so a Ventum reticle can be hold-solved by the
//! existing [`hold_point_in_reticle`](crate::reticle::hold_point_in_reticle) without any of the CLI/FFI/WASM surfaces having to
//! learn a second schema. It is a converter, not new physics — every coordinate is carried
//! straight through, only unit-normalized to milliradians and repeat-expanded.
//!
//! # What is kept, what is dropped
//!
//! The Ventum `spec` mixes *hold-bearing* marks with *decoration*:
//!
//! * `dot` and `tick` are the holdable marks. A `dot` becomes a [`MarkKind::Dot`]; a `tick`
//!   becomes a [`MarkKind::Hash`] whose hold point is its own `(x, y)` anchor.
//! * `text` is a standalone label positioned in reticle coordinates. Following the schema's
//!   recommended "option (b)", each text is bound to the nearest hold-bearing mark within
//!   [`DEFAULT_TEXT_BIND_MIL`] and becomes that mark's `label`; unattached text is dropped.
//! * `line`, `circle`, `rect`, `grid`, and any unknown future element type are pure
//!   decoration and are dropped.
//!
//! **Nothing is dropped silently.** [`import_ventum_reticle_with_report`] returns a
//! [`VentumImportReport`] tallying every element that produced no hold point, so a caller can
//! say "this document drew 6 things I cannot aim with" instead of handing back a reticle that
//! is quietly missing hold points and looks merely sparse. [`import_ventum_reticle`] is that
//! function with the report discarded; it exists for callers that genuinely do not care.
//!
//! # Coordinate and unit conventions
//!
//! Ventum uses `+x = right`, `+y = down`, origin at the reticle center — identical to the
//! engine's [`ReticleMark::right_mil`] / [`ReticleMark::down_mil`], so no sign flips are
//! needed. All coordinates are expressed in the reticle's own `unit`; MOA input is
//! converted to milliradians with `1 MOA = 0.2908882 mil` ([`MOA_TO_MIL`]). Auto-numbered
//! ladder labels are NOT unit-converted — a label is the reticle-unit value it names.
//!
//! Confirmed by the format's author (2026-08-01): `+y` is down because a reticle is numbered
//! downward (a 4-unit holdover is `y: 4`; a stadia mark 5 up is `y: -5`), which is plain SVG
//! screen convention and happens to match shooter intuition.
//!
//! ## Arc angles
//!
//! Ventum draws arcs as a `circle` element carrying `start`/`end` angles. Those angles are
//! measured from 3 o'clock (`0° = +x`, right) and sweep **CLOCKWISE**:
//!
//! ```text
//!   0° = right      90° = down      180° = left      270° = UP
//! ```
//!
//! The one counterintuitive consequence, worth stating because it is the opposite of a
//! compass: since `+y` is down, the TOP of the reticle is 270°, not 90°. A horseshoe opening
//! downward is therefore `start: 200, end: 340`. This is the same system as the coordinates
//! (`+x` right, `+y` down, 0° at 3 o'clock, clockwise), so the format never mixes conventions.
//!
//! A point on an arc of radius `r` about the reticle center is therefore
//!
//! ```text
//!   x = r * cos(θ)          y = r * sin(θ)          (θ in degrees, +y DOWN)
//! ```
//!
//! Note the PLUS on the sine: in a y-up math convention this term is negated, and that single
//! sign is the whole trap. Verified against the author's reference diagram — its 0/90/180/270
//! markers and the `start: 200, end: 340` horseshoe endpoints reproduce exactly under the
//! formula above, and the arc's 140° clockwise sweep passes through 270° (the top), leaving the
//! gap at the bottom.
//!
//! ## Why an arc still contributes no mark (MBA-1441)
//!
//! A horseshoe's apex and its two tips are real aiming references on Vortex-style reticles,
//! so "does an arc carry a hold?" has no single answer — and the Ventum format gives the
//! document no way to say. A `circle` with a sweep is spelled identically whether it is a
//! ranging horseshoe or a decorative ring segment, and turning every one of them into three
//! marks would do two bad things at once: fabricate aiming points a reticle may not have, and
//! silently change what [`hold_point_in_reticle`](crate::reticle::hold_point_in_reticle)
//! answers for every Ventum reticle already imported since 0.32.0, since the nearest mark it
//! reports would start snapping to invented geometry. Inventing holds out of decoration is the
//! same sin as dropping holds without saying so.
//!
//! So an arc is still not a mark — but it is no longer a silence either. Each one is resolved
//! (using exactly the convention above, which is why the convention was recorded) into the
//! three points a shooter could actually index on, and handed back on the report as a
//! [`VentumArc`]: its two tips and its apex, in the engine's own `right_mil` / `down_mil`. A
//! caller who knows their horseshoe is hold-bearing turns them into marks with three
//! [`ReticleMark::new`] calls and never re-derives the clockwise/`+y`-down trap; a caller who
//! does not, at least learns the document drew something they cannot aim with.
//!
//! # Safety
//!
//! The mark cap ([`crate::reticle::MAX_RETICLE_MARKS`]) is enforced *during* repeat
//! expansion, so a hostile `repeat.n` (or a huge mirrored ladder) can never allocate an
//! unbounded vector: expansion stops the instant the running instance count would exceed
//! the cap, having materialized at most one instance past it. This is the same defense as
//! the 0.31.0 reticle generator size guards.

use serde::de::{self, SeqAccess, Visitor};
use serde::{Deserialize, Deserializer};
use std::collections::BTreeMap;
use std::fmt;

use crate::reticle::{
    FocalPlane, MarkKind, ReticleDescription, ReticleError, ReticleMark, MAX_RETICLE_MARKS,
};

/// Milliradians per minute of angle (`1 MOA = 0.2908882 mil`). Ventum reticles authored in
/// MOA are converted to the engine's milliradian marks with this factor.
pub const MOA_TO_MIL: f64 = 0.2908882;

/// Largest `(x, y)` distance (in milliradians) at which a `text` element binds to a
/// hold-bearing mark and becomes its label. Text farther than this from every mark is
/// treated as free-floating decoration and dropped.
pub const DEFAULT_TEXT_BIND_MIL: f64 = 1.0;

/// The tag [`VentumImportReport::dropped_element_types`] uses for a `circle` that carries a
/// sweep. The format has no `arc` element — an arc IS a `circle` with `start`/`end` — but the
/// two are worth telling apart in a report, because only one of them might have been a hold.
pub const ARC_TAG: &str = "arc";

/// How near a mirrored arc's angles must land to the original's for the reflection to count
/// as having mapped the arc onto itself (degrees). The comparison is of sums of the
/// document's own angle literals, so the only error to absorb is the last bit or two of the
/// subtraction — this is a float-equality guard, not a modelling tolerance.
const ARC_SYMMETRY_EPSILON_DEGREES: f64 = 1e-9;

/// A point on an imported arc, in the engine's own mark coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct VentumArcPoint {
    /// Milliradians right of the optical center (negative = left).
    pub right_mil: f64,
    /// Milliradians below the optical center (negative = above).
    pub down_mil: f64,
}

/// One `circle` element that carried a sweep — a horseshoe or other arc — resolved into the
/// points a shooter could index on, but NOT imported as marks. See the module documentation
/// for why the decision is left to the caller.
///
/// Every field is in milliradians from the optical center (angles excepted, which are the
/// document's own degrees), so
/// `ReticleMark::new(arc.apex.down_mil, arc.apex.right_mil, MarkKind::Dot)` is all it takes to
/// adopt one.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct VentumArc {
    /// The arc's center — the point its radius is measured from, NOT an aiming point.
    pub center: VentumArcPoint,
    /// Arc radius in milliradians.
    pub radius_mil: f64,
    /// The document's `start` angle, degrees, measured from 3 o'clock and sweeping clockwise
    /// (so 270° is the TOP of the reticle — see the module documentation).
    pub start_degrees: f64,
    /// The document's `end` angle, same convention.
    pub end_degrees: f64,
    /// How far the arc sweeps clockwise from `start` to `end`, in `(0, 360)` degrees. A
    /// `circle` whose angles describe a whole revolution is a ring, not an arc, and never
    /// appears here.
    pub sweep_degrees: f64,
    /// Where the arc begins: the visible tip at `start_degrees`.
    pub start_tip: VentumArcPoint,
    /// The sweep midpoint — a horseshoe's apex, the closed end opposite its gap. The one
    /// point on an arc a shooter can index without measuring.
    pub apex: VentumArcPoint,
    /// Where the arc ends: the visible tip at `end_degrees`.
    pub end_tip: VentumArcPoint,
}

/// What an import could not turn into a hold point, so a sparse reticle is explained rather
/// than merely suspicious. Returned by [`import_ventum_reticle_with_report`].
///
/// # How the counts are taken
///
/// `arc` and `circle` entries count *expanded* instances, because this report resolves each
/// instance's own geometry and a mirrored pair of horseshoes is genuinely two of them.
/// `line`, `rect`, `grid`, `text` and unknown types count elements *as the document writes
/// them* — a `repeat` on one of those counts once, since the report carries no geometry to
/// distinguish a dropped shape's copies from the shape itself.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct VentumImportReport {
    /// How many elements produced no hold point — the sum of [`Self::dropped_element_types`].
    /// Zero means the whole document is in the returned [`ReticleDescription`].
    pub dropped_elements: usize,
    /// Per-type counts of that drop, sorted by tag so the report is deterministic. Tags are
    /// the format's own element types, plus [`ARC_TAG`] for a `circle` carrying a sweep and
    /// `"unknown"` for an element type this module does not know. A `"text"` entry counts
    /// only UNBOUND text: text that landed on a mark became its label and was not dropped.
    pub dropped_element_types: Vec<(String, usize)>,
    /// Every dropped arc, resolved. See [`VentumArc`] and the module documentation.
    pub arcs: Vec<VentumArc>,
    /// Arcs that could not be resolved into points because the element omitted `r` or gave a
    /// non-finite coordinate. They are counted in [`Self::dropped_element_types`] under
    /// [`ARC_TAG`] like any other, so `arcs.len()` plus this equals that tally — the
    /// discrepancy is stated rather than left to be noticed.
    pub arcs_unresolved: usize,
}

impl VentumImportReport {
    /// Whether anything at all was dropped. A caller that only wants to know "should I warn
    /// the user" asks this.
    pub fn is_empty(&self) -> bool {
        self.dropped_elements == 0
    }

    /// The per-type tally as `"arc x1, line x3"`, for a one-line notice. Empty string when
    /// nothing was dropped.
    pub fn tally(&self) -> String {
        self.dropped_element_types
            .iter()
            .map(|(tag, count)| format!("{tag} x{count}"))
            .collect::<Vec<_>>()
            .join(", ")
    }
}

/// Import a Ventum reticle spec into a [`ReticleDescription`].
///
/// `json` is a single Ventum reticle object (see the module documentation for the fields).
/// The returned description is ready for [`hold_point_in_reticle`](crate::reticle::hold_point_in_reticle): its marks are in
/// milliradians from the optical center, `focal_plane` mirrors the spec's `plane`, and
/// `reference_magnification` carries the spec's `ref_magnification` for an SFP reticle
/// (`1.0` for FFP, where subtensions are magnification-independent).
///
/// This is a pure transform: it does not itself reject an empty or all-decoration reticle
/// (the resulting description simply carries no marks, which [`hold_point_in_reticle`](crate::reticle::hold_point_in_reticle) then
/// reports as [`ReticleError::NoMarks`]). It DOES enforce
/// [`crate::reticle::MAX_RETICLE_MARKS`] during expansion.
///
/// Use [`import_ventum_reticle_with_report`] when you want to know what the document drew
/// that could not become a hold point — which is most callers, since a reticle missing its
/// arcs looks sparse rather than wrong.
///
/// # Errors
///
/// Returns [`ReticleError::InvalidSpec`] if `json` cannot be parsed as a Ventum reticle,
/// and [`ReticleError::TooManyMarks`] if repeat expansion would exceed the mark cap.
pub fn import_ventum_reticle(json: &str) -> Result<ReticleDescription, ReticleError> {
    import_ventum_reticle_with_report(json).map(|(description, _)| description)
}

/// Import a Ventum reticle spec, also returning a [`VentumImportReport`] of everything the
/// hold/decoration split left behind.
///
/// Identical to [`import_ventum_reticle`] in every other respect; that function is this one
/// with the report discarded. The description is byte-for-byte what it has returned since
/// 0.32.0 — resolving arcs added a report, not a mark.
///
/// # Errors
///
/// As [`import_ventum_reticle`].
pub fn import_ventum_reticle_with_report(
    json: &str,
) -> Result<(ReticleDescription, VentumImportReport), ReticleError> {
    let reticle: VentumReticle =
        serde_json::from_str(json).map_err(|e| ReticleError::InvalidSpec(e.to_string()))?;

    let scale = match reticle.unit {
        Unit::Mil => 1.0,
        Unit::Moa => MOA_TO_MIL,
    };

    // Expand every element's `repeat` into concrete, unit-scaled instances first, with the
    // mark cap enforced during expansion.
    let mut dropped: BTreeMap<&'static str, usize> = BTreeMap::new();
    let instances = expand_elements(&reticle.spec.0, scale, &mut dropped)?;

    // Classify: dot/tick become marks; text is collected for nearest-mark binding; a circle
    // is decoration that gets counted, and an arc is decoration that gets counted AND
    // resolved.
    let mut marks: Vec<ReticleMark> = Vec::new();
    let mut texts: Vec<TextLabel> = Vec::new();
    let mut arcs: Vec<VentumArc> = Vec::new();
    let mut arcs_unresolved = 0usize;
    for instance in instances {
        match instance.role {
            Role::Dot => marks.push(ReticleMark::new(instance.down, instance.right, MarkKind::Dot)),
            Role::Tick => {
                marks.push(ReticleMark::new(instance.down, instance.right, MarkKind::Hash))
            }
            Role::Text(label) => texts.push(TextLabel {
                down: instance.down,
                right: instance.right,
                label,
            }),
            Role::Circle(circle) => {
                let tag = if circle.sweep_degrees().is_some() {
                    ARC_TAG
                } else {
                    "circle"
                };
                *dropped.entry(tag).or_insert(0) += 1;
                if tag == ARC_TAG {
                    match circle.resolve(instance.down, instance.right) {
                        Some(arc) => arcs.push(arc),
                        None => arcs_unresolved += 1,
                    }
                }
            }
        }
    }

    // Text that found no mark inside the bind radius is dropped, and that is a lost label,
    // not a lost shape — count it here, where binding is what decides it.
    let unbound = bind_text_labels(&mut marks, &texts);
    if unbound > 0 {
        *dropped.entry("text").or_insert(0) += unbound;
    }

    // `ref_magnification` is meaningful only for SFP; an FFP reticle's subtensions do not
    // depend on magnification, so we normalize its reference to 1.0 (as the engine does).
    let reference_magnification = match reticle.plane {
        FocalPlane::Second => reticle.ref_magnification,
        FocalPlane::First => 1.0,
    };

    let report = VentumImportReport {
        dropped_elements: dropped.values().sum(),
        dropped_element_types: dropped
            .into_iter()
            .map(|(tag, count)| (tag.to_string(), count))
            .collect(),
        arcs,
        arcs_unresolved,
    };

    Ok((
        ReticleDescription {
            name: reticle.name,
            focal_plane: reticle.plane,
            reference_magnification,
            marks,
        },
        report,
    ))
}

/// A `text` element after expansion: its position (already unit-scaled, in milliradians)
/// and the label string it would apply to the nearest mark.
struct TextLabel {
    down: f64,
    right: f64,
    label: String,
}

/// Bind each collected text to the nearest mark within [`DEFAULT_TEXT_BIND_MIL`], setting
/// that mark's label. Text with no mark inside the cap is dropped. When two texts bind to
/// the same mark the later one wins.
///
/// Returns how many texts bound to nothing, so MBA-1441's report can say a label was lost
/// rather than leave a mark looking as though the document never named it. A reticle with no
/// marks at all drops every text, which is the same answer by a shorter route.
fn bind_text_labels(marks: &mut [ReticleMark], texts: &[TextLabel]) -> usize {
    let mut unbound = 0usize;
    for text in texts {
        let mut nearest: Option<usize> = None;
        let mut nearest_distance = f64::INFINITY;
        for (index, mark) in marks.iter().enumerate() {
            let distance = ((text.right - mark.right_mil).powi(2)
                + (text.down - mark.down_mil).powi(2))
            .sqrt();
            if distance < nearest_distance {
                nearest_distance = distance;
                nearest = Some(index);
            }
        }
        match nearest {
            Some(index) if nearest_distance <= DEFAULT_TEXT_BIND_MIL => {
                marks[index].label = Some(text.label.clone());
            }
            _ => unbound += 1,
        }
    }
    unbound
}

/// The role an expanded point instance plays once classified.
#[derive(Debug, Clone, PartialEq)]
enum Role {
    Dot,
    Tick,
    /// A text label carrying its (possibly auto-numbered) string.
    Text(String),
    /// A `circle` element, with the sweep that makes it an arc when it has one. Not a hold —
    /// carried through expansion so MBA-1441's report can count and resolve each instance.
    Circle(CircleShape),
}

impl Role {
    /// This role as it appears in a [`Repeat`]'s mirrored twin.
    ///
    /// Only an arc changes: mirroring reflects its angles AND reverses its direction of
    /// travel, so the clockwise sweep `start -> end` becomes the clockwise sweep
    /// `f(end) -> f(start)`. Reflecting about the vertical axis (a stepped `x`, so
    /// `right -> -right`) maps `theta -> 180 - theta`; about the horizontal axis (a stepped
    /// `y`, so `down -> -down`) it maps `theta -> -theta`. Both fall straight out of
    /// `(cos t, sin t)` under the module's `+y`-down convention; neither is a guess.
    ///
    /// Dots, ticks and text are unchanged, which is what keeps mirrored ladder labels
    /// unsigned (see [`expand_point`]).
    fn mirrored(&self, axis: Axis) -> Role {
        match self {
            Role::Circle(circle) => Role::Circle(circle.mirrored(axis)),
            other => other.clone(),
        }
    }
}

/// A `circle` element's geometry, already unit-scaled. The center is carried by the
/// [`ExpandedInstance`] around this, so only the radius and the optional sweep live here.
#[derive(Debug, Clone, Copy, PartialEq)]
struct CircleShape {
    /// Radius in milliradians, or `None` when the element omitted `r` (or gave a value that
    /// is not a usable length). Such an element is still counted; it just cannot be resolved.
    radius_mil: Option<f64>,
    /// The document's `start`/`end` angles in degrees, when it declared both. `None` is a
    /// plain ring, which is decoration under any reading and gets no [`VentumArc`].
    angles: Option<(f64, f64)>,
}

impl CircleShape {
    /// How far this circle sweeps clockwise, in `(0, 360)` degrees, or `None` when it is a
    /// full ring — either because it declared no angles, or because the angles it declared
    /// describe a whole revolution (`start == end`, or a multiple of 360 apart).
    fn sweep_degrees(&self) -> Option<f64> {
        let (start, end) = self.angles?;
        if !start.is_finite() || !end.is_finite() {
            return None;
        }
        let sweep = (end - start).rem_euclid(360.0);
        // A zero remainder is a closed circle drawn as an arc, not a degenerate arc.
        (sweep > 0.0).then_some(sweep)
    }

    /// Whether reflecting this shape about `axis` maps it onto ITSELF, so a mirrored twin
    /// sharing its center would be a duplicate rather than a second shape.
    ///
    /// This is the arc's half of [`expand_point`]'s center dedupe. For a mark, "the mirror
    /// landed on me" is the whole test; for an arc it is only half of it, because reflection
    /// also reworks the angles ([`Role::mirrored`]). A ring maps onto itself under any
    /// reflection through its own center. A swept arc does so exactly when its two endpoints
    /// swap into each other: reflection about the vertical axis sends `theta -> 180 - theta`,
    /// so it needs `start + end == 180`; about the horizontal axis it sends `theta -> -theta`,
    /// so it needs `start + end == 0` — both modulo a full revolution. A horseshoe centered on
    /// the vertical axis (`start: 200, end: 340`, sum 540 ≡ 180) is symmetric about it and is
    /// deduped; the same horseshoe mirrored vertically opens the other way and is not.
    fn is_reflection_of_itself(&self, axis: Axis) -> bool {
        // A ring — no angles, or angles spanning a whole revolution — is its own reflection.
        let Some((start, end)) = self.angles.filter(|_| self.sweep_degrees().is_some()) else {
            return true;
        };
        let target = match axis {
            Axis::X => 180.0,
            Axis::Y => 0.0,
        };
        let offset = (start + end - target).rem_euclid(360.0);
        // Distance to the nearest multiple of 360, so 0 and 360 are both "no rotation".
        offset.min(360.0 - offset) <= ARC_SYMMETRY_EPSILON_DEGREES
    }

    /// This shape reflected for a mirrored [`Repeat`] twin. See [`Role::mirrored`] for the
    /// derivation; note the swap, which is the direction reversal.
    fn mirrored(&self, axis: Axis) -> CircleShape {
        let angles = self.angles.map(|(start, end)| match axis {
            Axis::X => (180.0 - end, 180.0 - start),
            Axis::Y => (-end, -start),
        });
        CircleShape {
            radius_mil: self.radius_mil,
            angles,
        }
    }

    /// Resolve this circle, centered at `(down, right)` milliradians, into the arc points a
    /// caller could adopt as holds. `None` when it is a ring, or when the element gave no
    /// usable radius or center — the caller counts those separately rather than reporting a
    /// zero-radius arc whose apex and tips all sit on the center.
    fn resolve(&self, down: f64, right: f64) -> Option<VentumArc> {
        let sweep = self.sweep_degrees()?;
        let radius = self.radius_mil?;
        if !radius.is_finite() || radius <= 0.0 || !down.is_finite() || !right.is_finite() {
            return None;
        }
        let (start, end) = self.angles?;

        // x = r*cos(theta), y = r*sin(theta) about the center, degrees, +y DOWN — the module
        // documentation's formula, verified against the format author's reference diagram.
        // The PLUS on the sine is the whole trap: a y-up convention negates it.
        let point = |degrees: f64| {
            let radians = degrees.to_radians();
            VentumArcPoint {
                right_mil: right + radius * radians.cos(),
                down_mil: down + radius * radians.sin(),
            }
        };

        Some(VentumArc {
            center: VentumArcPoint {
                right_mil: right,
                down_mil: down,
            },
            radius_mil: radius,
            start_degrees: start,
            end_degrees: end,
            sweep_degrees: sweep,
            start_tip: point(start),
            // Halfway along the sweep, travelling clockwise from `start` — for a horseshoe,
            // the closed end opposite its gap.
            apex: point(start + sweep / 2.0),
            end_tip: point(end),
        })
    }
}

/// One expanded point instance: a hold-bearing mark or a text label, positioned in
/// milliradians from center (unit scale already applied).
#[derive(Debug, Clone, PartialEq)]
struct ExpandedInstance {
    /// Milliradians below the optical center.
    down: f64,
    /// Milliradians right of the optical center.
    right: f64,
    role: Role,
}

/// Which point-element kind is being expanded, borrowing the base text for `text`.
enum PointRole<'a> {
    Dot,
    Tick,
    Text(&'a str),
    Circle(CircleShape),
}

impl PointRole<'_> {
    /// The [`Role`] for one instance. `ladder_label` is `Some` only for an auto-numbered
    /// text ladder copy; otherwise a text keeps its element's own string.
    fn role_with_label(&self, ladder_label: Option<String>) -> Role {
        match self {
            PointRole::Dot => Role::Dot,
            PointRole::Tick => Role::Tick,
            PointRole::Text(base) => Role::Text(ladder_label.unwrap_or_else(|| (*base).to_string())),
            PointRole::Circle(shape) => Role::Circle(*shape),
        }
    }

    /// The auto-numbered label for copy `i` of a text ladder, or `None` when this is not a
    /// labeled text ladder (dot/tick, or text without `label:true`).
    fn ladder_label(&self, repeat: &Repeat, i: u32) -> Option<String> {
        match self {
            PointRole::Text(_) if repeat.label => {
                Some(format_label_number(repeat.label_start + f64::from(i) * repeat.label_step))
            }
            _ => None,
        }
    }
}

/// Expand every element's `repeat` into concrete instances, applying `scale` to positions
/// and enforcing the mark cap along the way.
///
/// MBA-1441: nothing is dropped without a tally any more. `line`/`rect`/`grid`/unknown carry
/// no geometry this module reads, so each is counted once into `dropped` and skipped;
/// `circle` goes through the expander like a point element, because an arc's resolved apex
/// and tips are per-instance and a mirrored pair of horseshoes is two distinct arcs.
fn expand_elements(
    elements: &[Element],
    scale: f64,
    dropped: &mut BTreeMap<&'static str, usize>,
) -> Result<Vec<ExpandedInstance>, ReticleError> {
    let mut out: Vec<ExpandedInstance> = Vec::new();
    for element in elements {
        let (x, y, repeat, base) = match element {
            Element::Dot { x, y, repeat } => (*x, *y, repeat.as_ref(), PointRole::Dot),
            Element::Tick { x, y, repeat } => (*x, *y, repeat.as_ref(), PointRole::Tick),
            Element::Text {
                x, y, text, repeat,
            } => (*x, *y, repeat.as_ref(), PointRole::Text(text)),
            Element::Circle {
                x,
                y,
                r,
                start,
                end,
                repeat,
            } => (
                x.unwrap_or(0.0),
                y.unwrap_or(0.0),
                repeat.as_ref(),
                PointRole::Circle(CircleShape {
                    // The radius is a length in the reticle's unit, so it scales with the
                    // coordinates; the angles are angles and do not.
                    radius_mil: r.map(|r| r * scale),
                    angles: start.zip(*end),
                }),
            ),
            // Decoration this module reads no geometry from, and element types it has never
            // heard of. Counted, then skipped.
            Element::Line {} => {
                *dropped.entry("line").or_insert(0) += 1;
                continue;
            }
            Element::Rect {} => {
                *dropped.entry("rect").or_insert(0) += 1;
                continue;
            }
            Element::Grid {} => {
                *dropped.entry("grid").or_insert(0) += 1;
                continue;
            }
            Element::Unknown => {
                *dropped.entry("unknown").or_insert(0) += 1;
                continue;
            }
        };
        expand_point(x, y, repeat, scale, base, &mut out)?;
    }
    Ok(out)
}

/// Expand a single point element (`dot`/`tick`/`text`) — with or without a `repeat` — into
/// `out`, applying `scale` and enforcing the cap.
fn expand_point(
    x: f64,
    y: f64,
    repeat: Option<&Repeat>,
    scale: f64,
    base: PointRole<'_>,
    out: &mut Vec<ExpandedInstance>,
) -> Result<(), ReticleError> {
    // Ventum uses +x = right, +y = down — the same axes the engine's marks use — so a point
    // at (x, y) maps to (down = y, right = x), unit-scaled.
    let Some(repeat) = repeat else {
        // No repeat: a single instance at (x, y).
        let role = base.role_with_label(None);
        return push_capped(y * scale, x * scale, role, out);
    };

    for i in 0..repeat.n {
        let offset = f64::from(i) * repeat.step;
        // Step the chosen axis; the other axis is unchanged.
        let (cx, cy) = match repeat.axis {
            Axis::X => (x + offset, y),
            Axis::Y => (x, y + offset),
        };
        // The coordinate on the stepped axis decides whether a mirrored copy would land on
        // center (and so must be skipped, to avoid a duplicate).
        let stepped = match repeat.axis {
            Axis::X => cx,
            Axis::Y => cy,
        };

        let down = cy * scale;
        let right = cx * scale;
        let role = base.role_with_label(base.ladder_label(repeat, i));

        // MBA-1441: a mirrored twin is skipped only when it would be a DUPLICATE, and being
        // on the mirror line is not by itself enough to make it one. It is for a mark, whose
        // whole identity is its position — but an arc is also a set of angles, and reflection
        // reworks those, so a centered ASYMMETRIC arc's twin is a genuinely different shape
        // (`start: 290, end: 70` reflects to an apex on the other side of the same center).
        // Counting one where the document drew two is the same undercount this ticket is
        // named for. The old test is exactly the position half of this one, so nothing about
        // dot/tick/text expansion changes.
        let mirror_is_distinct = stepped != 0.0
            || match &role {
                Role::Circle(circle) => !circle.is_reflection_of_itself(repeat.axis),
                Role::Dot | Role::Tick | Role::Text(_) => false,
            };

        if repeat.mirror && mirror_is_distinct {
            // Emit the axis-negated twin as well, but never a duplicate at the center. Only
            // the stepped axis is negated. A mirrored labeled-text copy keeps the SAME label
            // as its positive twin (magnitude convention: the "5" on the right and the "5"
            // on the left both read 5). The written schema does not specify this, so it began
            // as a deliberate assumption; the format's author confirmed it on 2026-08-01
            // ("mirrored numbers are without signs"). Do not "fix" this into signed labels.
            let (mirror_down, mirror_right) = match repeat.axis {
                Axis::X => (down, -right),
                Axis::Y => (-down, right),
            };
            // MBA-1441: an arc's ANGLES have to be reflected along with its center, or a
            // mirrored horseshoe reports an apex on the wrong side of itself. Every other
            // role mirrors to itself, so this is a no-op for the marks.
            let mirror_role = role.mirrored(repeat.axis);
            push_capped(down, right, role, out)?;
            push_capped(mirror_down, mirror_right, mirror_role, out)?;
        } else {
            push_capped(down, right, role, out)?;
        }
    }
    Ok(())
}

/// Push one instance, rejecting the reticle the instant the running count would exceed
/// [`MAX_RETICLE_MARKS`]. Because this checks before every push, expansion allocates at
/// most `MAX_RETICLE_MARKS + 1` instances regardless of how large a `repeat.n` is, so a
/// runaway ladder returns promptly instead of hanging or exhausting memory. Text instances
/// count toward the same cap as marks, so no element kind can grow unbounded.
fn push_capped(
    down: f64,
    right: f64,
    role: Role,
    out: &mut Vec<ExpandedInstance>,
) -> Result<(), ReticleError> {
    if out.len() >= MAX_RETICLE_MARKS {
        return Err(ReticleError::TooManyMarks {
            count: out.len() + 1,
            max: MAX_RETICLE_MARKS,
        });
    }
    out.push(ExpandedInstance { down, right, role });
    Ok(())
}

/// Format a ladder label value so integers print without a fractional part (`5.0 -> "5"`)
/// while fractional values keep their shortest decimal form (`5.5 -> "5.5"`).
fn format_label_number(value: f64) -> String {
    if value.is_finite() && value == value.round() {
        (value as i64).to_string()
    } else {
        value.to_string()
    }
}

// ---------------------------------------------------------------------------------------
// Ventum input model (serde). Deliberately permissive: no `deny_unknown_fields`, so
// cosmetic keys (color, width, cx/cy, size, len, orient, r, max_extent, tube_diameter,
// notes, manufacturer, ...) are ignored rather than rejected. "Ignored" includes the
// cosmetic keys this module DOES read: a `circle`'s x/y/r/start/end go through
// [`lenient_f64`], so a value of the wrong JSON type degrades that one field to `None`
// instead of failing the whole document. Only a hold-bearing coordinate — a `dot`/`tick`
// `x`/`y`, or a `text`'s string — is strict, because there is no sane fallback for a mark
// whose position cannot be read.
// ---------------------------------------------------------------------------------------

/// A whole Ventum reticle: metadata plus its drawing `spec`.
#[derive(Debug, Deserialize)]
struct VentumReticle {
    #[serde(default)]
    name: String,
    /// `"ffp"` / `"sfp"` — deserialized directly into the engine's [`FocalPlane`].
    #[serde(default)]
    plane: FocalPlane,
    /// The magnification the subtensions are true at (SFP only).
    #[serde(default = "default_ref_magnification")]
    ref_magnification: f64,
    /// The angular unit for every coordinate in `spec`.
    #[serde(default)]
    unit: Unit,
    #[serde(default)]
    spec: Spec,
}

fn default_ref_magnification() -> f64 {
    1.0
}

/// The angular unit of a Ventum reticle's coordinates.
#[derive(Debug, Clone, Copy, Default, Deserialize)]
#[serde(rename_all = "snake_case")]
enum Unit {
    #[default]
    Mil,
    Moa,
}

/// The `spec` field: a JSON array of [`Element`]s, given EITHER as a JSON array OR as a
/// string containing that array (the schema shows both forms). The custom deserializer
/// accepts either and parses the string form with `serde_json`.
#[derive(Debug, Default)]
struct Spec(Vec<Element>);

impl<'de> Deserialize<'de> for Spec {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct SpecVisitor;

        impl<'de> Visitor<'de> for SpecVisitor {
            type Value = Spec;

            fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str("a JSON array of reticle elements, or a string containing that array")
            }

            fn visit_str<E>(self, v: &str) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                let elements: Vec<Element> =
                    serde_json::from_str(v).map_err(de::Error::custom)?;
                Ok(Spec(elements))
            }

            fn visit_seq<A>(self, seq: A) -> Result<Self::Value, A::Error>
            where
                A: SeqAccess<'de>,
            {
                let elements =
                    Vec::<Element>::deserialize(de::value::SeqAccessDeserializer::new(seq))?;
                Ok(Spec(elements))
            }
        }

        deserializer.deserialize_any(SpecVisitor)
    }
}

/// Read a number the module treats as cosmetic geometry, never failing on one it cannot use.
///
/// `Option<f64>` alone absorbs an ABSENT or `null` field but still *rejects* a present field
/// of the wrong type — and one rejected field fails the entire document, not the element that
/// carried it. That is the wrong trade for a key this module only ever reads as a nicety: a
/// Ventum tool is free to write `"r": "2mil"`, `"r": {"v": 2, "unit": "mil"}` or anything else
/// its UI finds convenient, and none of that is a reason to refuse a reticle whose dots are
/// perfectly good. So the value is read as arbitrary JSON and kept only if it is a finite
/// number; everything else becomes `None`, i.e. a counted-but-unresolved decoration, which is
/// the same answer the field being absent has always given.
fn lenient_f64<'de, D>(deserializer: D) -> Result<Option<f64>, D::Error>
where
    D: Deserializer<'de>,
{
    let value = serde_json::Value::deserialize(deserializer)?;
    Ok(value.as_f64().filter(|v| v.is_finite()))
}

/// One drawing element, internally tagged on `"type"`. Unknown cosmetic fields on each
/// variant are ignored; unknown element *types* map to [`Element::Unknown`].
#[derive(Debug, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
enum Element {
    /// A filled dot (hold-bearing). `r` and any other cosmetic fields are ignored.
    Dot {
        x: f64,
        y: f64,
        #[serde(default)]
        repeat: Option<Repeat>,
    },
    /// A short hash mark (hold-bearing); its hold point is its `(x, y)` anchor. `len` and
    /// `orient` are ignored.
    Tick {
        x: f64,
        y: f64,
        #[serde(default)]
        repeat: Option<Repeat>,
    },
    /// A standalone text label positioned in reticle coordinates. `size` is ignored.
    Text {
        x: f64,
        y: f64,
        text: String,
        #[serde(default)]
        repeat: Option<Repeat>,
    },
    /// Decoration — dropped. Its geometry fields are ignored (empty struct variant so
    /// serde skips every field).
    Line {},
    /// A ring, or — when it carries `start`/`end` — an arc such as a horseshoe. Decoration
    /// either way: still dropped, but no longer in silence (MBA-1441). Every field is
    /// optional AND leniently typed ([`lenient_f64`]), so a cosmetic ring that names none of
    /// them — or spells one as `"2mil"`, an object, or `null` — still parses, exactly as it
    /// did when this variant read nothing at all.
    Circle {
        /// Center, `+x` right. `cx` is the schema's other spelling of the same field.
        #[serde(default, alias = "cx", deserialize_with = "lenient_f64")]
        x: Option<f64>,
        /// Center, `+y` down.
        #[serde(default, alias = "cy", deserialize_with = "lenient_f64")]
        y: Option<f64>,
        /// Radius, in the reticle's own unit.
        #[serde(default, deserialize_with = "lenient_f64")]
        r: Option<f64>,
        /// Sweep start in degrees from 3 o'clock, clockwise (see the module documentation —
        /// 270 is the TOP of the reticle). Present only on an arc.
        #[serde(default, deserialize_with = "lenient_f64")]
        start: Option<f64>,
        /// Sweep end, same convention. An arc needs BOTH: one angle alone describes no sweep.
        #[serde(default, deserialize_with = "lenient_f64")]
        end: Option<f64>,
        #[serde(default)]
        repeat: Option<Repeat>,
    },
    /// Decoration — dropped.
    Rect {},
    /// Decoration — dropped.
    Grid {},
    /// Any unknown future element type — dropped.
    #[serde(other)]
    Unknown,
}

/// The `repeat` operator: stamp `n` copies of an element along one axis, optionally
/// mirrored, optionally auto-numbering a text ladder's labels.
#[derive(Debug, Deserialize)]
struct Repeat {
    /// Which axis to step along.
    axis: Axis,
    /// Spacing between copies, in the reticle's unit.
    step: f64,
    /// Number of copies.
    n: u32,
    /// When true, also emit the axis-negated copy (a `±` symmetric ladder) — except where
    /// that copy would be a duplicate of the original, which for a mark means landing on the
    /// mirror line and for an arc additionally means the reflection leaving its angles
    /// unchanged (see [`CircleShape::is_reflection_of_itself`]).
    #[serde(default)]
    mirror: bool,
    /// (text only) When true, rewrite each copy's label to `label_start + i*label_step`.
    #[serde(default)]
    label: bool,
    /// (text ladder) The value of the first copy's auto-number.
    #[serde(rename = "labelStart", default)]
    label_start: f64,
    /// (text ladder) The auto-number increment per copy.
    #[serde(rename = "labelStep", default = "default_label_step")]
    label_step: f64,
}

fn default_label_step() -> f64 {
    1.0
}

/// The axis a [`Repeat`] steps along.
#[derive(Debug, Clone, Copy, Deserialize)]
#[serde(rename_all = "snake_case")]
enum Axis {
    X,
    Y,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reticle::hold_point_in_reticle;

    /// Bero's exact MBR dot-tree spec (the array form), as three `dot` rows each stamped by
    /// a horizontal mirrored `repeat`.
    const MBR_SPEC: &str = r#"[{"type":"dot","y":4,"x":1,"r":0.12,"repeat":{"axis":"x","step":1,"n":9,"mirror":true}},{"type":"dot","y":8,"x":1,"r":0.12,"repeat":{"axis":"x","step":1,"n":13,"mirror":true}},{"type":"dot","y":12,"x":1,"r":0.12,"repeat":{"axis":"x","step":1,"n":17,"mirror":true}}]"#;

    /// Collect the `right_mil` values of every mark on a given `down_mil` row, sorted.
    fn row_rights(desc: &ReticleDescription, down: f64) -> Vec<f64> {
        let mut rights: Vec<f64> = desc
            .marks
            .iter()
            .filter(|m| m.down_mil == down)
            .map(|m| m.right_mil)
            .collect();
        rights.sort_by(|a, b| a.partial_cmp(b).unwrap());
        rights
    }

    /// The set `{±1, ±2, ..., ±max}`, sorted ascending — the expected spread of one MBR row.
    fn symmetric_spread(max: i32) -> Vec<f64> {
        let mut v: Vec<f64> = Vec::new();
        for k in 1..=max {
            v.push(-f64::from(k));
        }
        for k in 1..=max {
            v.push(f64::from(k));
        }
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v
    }

    fn assert_mbr(desc: &ReticleDescription) {
        // 78 marks: row y=4 spans ±1..±9 (18), y=8 spans ±1..±13 (26), y=12 spans ±1..±17 (34).
        assert_eq!(desc.marks.len(), 78, "MBR expands to 78 holdable marks");
        assert_eq!(row_rights(desc, 4.0), symmetric_spread(9), "row down=4 -> ±1..±9");
        assert_eq!(row_rights(desc, 8.0), symmetric_spread(13), "row down=8 -> ±1..±13");
        assert_eq!(row_rights(desc, 12.0), symmetric_spread(17), "row down=12 -> ±1..±17");
        assert_eq!(row_rights(desc, 4.0).len(), 18);
        assert_eq!(row_rights(desc, 8.0).len(), 26);
        assert_eq!(row_rights(desc, 12.0).len(), 34);
        // Every mark is a dot, none sits on the vertical axis, and none exceeds its row span.
        assert!(desc.marks.iter().all(|m| m.kind == MarkKind::Dot), "all MBR marks are dots");
        assert!(
            desc.marks.iter().all(|m| m.right_mil != 0.0),
            "no MBR mark lands on the center line"
        );
        assert!(
            row_rights(desc, 4.0).iter().all(|&r| r.abs() <= 9.0),
            "row down=4 has nothing beyond ±9"
        );
    }

    #[test]
    fn mbr_array_form_expands_to_the_full_dot_tree() {
        let json = format!(r#"{{"name":"MBR","plane":"ffp","unit":"mil","spec":{MBR_SPEC}}}"#);
        let desc = import_ventum_reticle(&json).unwrap();
        assert_eq!(desc.name, "MBR");
        assert_eq!(desc.focal_plane, FocalPlane::First);
        assert_mbr(&desc);
    }

    #[test]
    fn mbr_string_form_parses_identically_to_the_array_form() {
        // `spec` given as a STRING containing the JSON array.
        let string_form = serde_json::json!({
            "name": "MBR", "plane": "ffp", "unit": "mil", "spec": MBR_SPEC
        })
        .to_string();
        let array_form = format!(r#"{{"name":"MBR","plane":"ffp","unit":"mil","spec":{MBR_SPEC}}}"#);

        let from_string = import_ventum_reticle(&string_form).unwrap();
        let from_array = import_ventum_reticle(&array_form).unwrap();

        assert_mbr(&from_string);
        assert_eq!(
            from_string, from_array,
            "the string-wrapped spec must import identically to the array spec"
        );
    }

    #[test]
    fn sfp_carries_reference_magnification_but_ffp_normalizes_it() {
        let sfp = import_ventum_reticle(
            r#"{"name":"SE","plane":"sfp","ref_magnification":6,"unit":"mil",
                "spec":[{"type":"dot","x":0,"y":1}]}"#,
        )
        .unwrap();
        assert_eq!(sfp.focal_plane, FocalPlane::Second);
        assert_eq!(sfp.reference_magnification, 6.0);

        // An FFP reticle ignores any ref_magnification and normalizes to 1.0.
        let ffp = import_ventum_reticle(
            r#"{"name":"F","plane":"ffp","ref_magnification":6,"unit":"mil",
                "spec":[{"type":"dot","x":0,"y":1}]}"#,
        )
        .unwrap();
        assert_eq!(ffp.focal_plane, FocalPlane::First);
        assert_eq!(ffp.reference_magnification, 1.0);
    }

    #[test]
    fn moa_coordinates_convert_to_milliradians() {
        let desc = import_ventum_reticle(
            r#"{"name":"M","unit":"moa","spec":[{"type":"dot","x":0,"y":2}]}"#,
        )
        .unwrap();
        assert_eq!(desc.marks.len(), 1);
        assert!(
            (desc.marks[0].down_mil - 2.0 * MOA_TO_MIL).abs() < 1e-6,
            "2 MOA -> {} mil (~0.581776)",
            desc.marks[0].down_mil
        );
    }

    #[test]
    fn text_binds_to_the_nearest_mark_and_far_text_is_dropped() {
        let desc = import_ventum_reticle(
            r#"{"name":"T","unit":"mil","spec":[
                {"type":"tick","x":0,"y":5},
                {"type":"text","x":0.5,"y":5,"text":"5"},
                {"type":"text","x":10,"y":10,"text":"far"}
            ]}"#,
        )
        .unwrap();
        // One tick -> one mark; the nearby "5" binds, the far "far" is dropped.
        assert_eq!(desc.marks.len(), 1);
        assert_eq!(desc.marks[0].down_mil, 5.0);
        assert_eq!(desc.marks[0].kind, MarkKind::Hash);
        assert_eq!(desc.marks[0].label.as_deref(), Some("5"));
        assert!(
            desc.marks.iter().all(|m| m.label.as_deref() != Some("far")),
            "distant text must not bind"
        );
    }

    #[test]
    fn mirror_never_duplicates_a_center_mark() {
        // A mirrored ladder whose stepped coordinate is 0 emits a single mark (no twin).
        let center = import_ventum_reticle(
            r#"{"name":"C","unit":"mil","spec":[
                {"type":"tick","x":0,"y":1,"repeat":{"axis":"x","step":1,"n":1,"mirror":true}}
            ]}"#,
        )
        .unwrap();
        assert_eq!(center.marks.len(), 1, "no duplicate at the center");
        assert_eq!(center.marks[0].right_mil, 0.0);

        // Off center, the mirror emits both ± copies.
        let off = import_ventum_reticle(
            r#"{"name":"O","unit":"mil","spec":[
                {"type":"tick","x":1,"y":1,"repeat":{"axis":"x","step":1,"n":1,"mirror":true}}
            ]}"#,
        )
        .unwrap();
        assert_eq!(off.marks.len(), 2, "mirror of an off-center mark yields ±1");
        let mut rights: Vec<f64> = off.marks.iter().map(|m| m.right_mil).collect();
        rights.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(rights, vec![-1.0, 1.0]);
    }

    #[test]
    fn text_ladder_auto_numbers_each_copy() {
        // Factored expander: a labeled text ladder renumbers each copy from labelStart by
        // labelStep (defaults 0 and 1).
        let elements = vec![Element::Text {
            x: 0.0,
            y: 1.0,
            text: "base".to_string(),
            repeat: Some(Repeat {
                axis: Axis::Y,
                step: 1.0,
                n: 3,
                mirror: false,
                label: true,
                label_start: 0.0,
                label_step: 1.0,
            }),
        }];
        let mut dropped = BTreeMap::new();
        let instances = expand_elements(&elements, 1.0, &mut dropped).unwrap();
        assert!(dropped.is_empty(), "a text ladder drops nothing");
        assert_eq!(
            instances,
            vec![
                ExpandedInstance { down: 1.0, right: 0.0, role: Role::Text("0".to_string()) },
                ExpandedInstance { down: 2.0, right: 0.0, role: Role::Text("1".to_string()) },
                ExpandedInstance { down: 3.0, right: 0.0, role: Role::Text("2".to_string()) },
            ]
        );
    }

    #[test]
    fn label_numbers_print_integers_without_a_decimal() {
        assert_eq!(format_label_number(5.0), "5");
        assert_eq!(format_label_number(0.0), "0");
        assert_eq!(format_label_number(-3.0), "-3");
        assert_eq!(format_label_number(5.5), "5.5");
    }

    #[test]
    fn a_runaway_repeat_is_capped_promptly() {
        let json = r#"{"name":"X","unit":"mil","spec":[
            {"type":"dot","x":0,"y":1,"repeat":{"axis":"x","step":1,"n":100000}}
        ]}"#;
        // Must reject with TooManyMarks without materializing 100000 marks (the cap is
        // enforced during expansion, so this returns immediately).
        assert!(matches!(
            import_ventum_reticle(json),
            Err(ReticleError::TooManyMarks { .. })
        ));
    }

    #[test]
    fn malformed_json_reports_invalid_spec() {
        assert!(matches!(
            import_ventum_reticle("{not json"),
            Err(ReticleError::InvalidSpec(_))
        ));
        // A `spec` string that is not itself a JSON array is also an invalid spec.
        assert!(matches!(
            import_ventum_reticle(r#"{"name":"X","spec":"not-an-array"}"#),
            Err(ReticleError::InvalidSpec(_))
        ));
    }

    #[test]
    fn decoration_and_unknown_elements_are_dropped() {
        let desc = import_ventum_reticle(
            r#"{"name":"D","unit":"mil","spec":[
                {"type":"line","x1":-5,"y1":0,"x2":5,"y2":0,"width":0.05},
                {"type":"circle","x":0,"y":0,"r":1.0},
                {"type":"rect","x":-1,"y":-1,"w":2,"h":2},
                {"type":"grid","x0":-5,"y0":-5,"x1":5,"y1":5,"step":1},
                {"type":"future_shape","x":0,"y":0},
                {"type":"dot","x":0,"y":2,"color":"red","width":0.1}
            ]}"#,
        )
        .unwrap();
        // Only the dot survives; every decoration and the unknown type are dropped.
        assert_eq!(desc.marks.len(), 1);
        assert_eq!(desc.marks[0].down_mil, 2.0);
        assert_eq!(desc.marks[0].kind, MarkKind::Dot);
    }

    // ----------------------------------------------------------------------------------
    // MBA-1441: dropped elements are reported, and an arc is resolved rather than guessed.
    // ----------------------------------------------------------------------------------

    /// The author's reference horseshoe: `start: 200, end: 340` opens DOWNWARD, so its apex
    /// is at 270 degrees, which is the TOP of the reticle — the one consequence of the
    /// clockwise/`+y`-down convention that is opposite to a compass. This test is the
    /// convention's regression guard: flip the sine's sign and the apex lands at the bottom.
    #[test]
    fn a_horseshoe_apex_resolves_to_the_top_of_the_reticle() {
        let (_, report) = import_ventum_reticle_with_report(
            r#"{"name":"H","unit":"mil","spec":[
                {"type":"dot","x":0,"y":1},
                {"type":"circle","x":0,"y":0,"r":2,"start":200,"end":340}
            ]}"#,
        )
        .unwrap();

        assert_eq!(report.arcs.len(), 1, "the horseshoe is reported");
        let arc = report.arcs[0];
        assert!((arc.sweep_degrees - 140.0).abs() < 1e-9, "200 -> 340 clockwise is 140 deg");

        // Apex at 270 deg, radius 2: (cos 270, sin 270) = (0, -1) with +y DOWN, so 2 mil UP.
        assert!(arc.apex.right_mil.abs() < 1e-9, "apex is on the vertical axis");
        assert!(
            (arc.apex.down_mil + 2.0).abs() < 1e-9,
            "apex must be 2 mil ABOVE center (down_mil -2), not below: {:?}",
            arc.apex
        );
        // Both tips sit below the apex and straddle it, leaving the gap at the bottom.
        assert!(arc.start_tip.right_mil < 0.0, "the 200 deg tip is on the left");
        assert!(arc.end_tip.right_mil > 0.0, "the 340 deg tip is on the right");
        for tip in [arc.start_tip, arc.end_tip] {
            assert!(
                tip.down_mil > arc.apex.down_mil,
                "a tip of a downward-opening horseshoe is below its apex: {tip:?}"
            );
        }
    }

    /// An arc is reported, not imported: the description is exactly what it was before this
    /// existed. Inventing three marks per horseshoe would silently move the nearest mark for
    /// every Ventum reticle imported since 0.32.0.
    #[test]
    fn an_arc_is_reported_but_never_becomes_a_mark() {
        let spec = |arc: &str| {
            format!(
                r#"{{"name":"A","unit":"mil","spec":[{{"type":"tick","x":0,"y":5}}{arc}]}}"#
            )
        };
        let without = import_ventum_reticle(&spec("")).unwrap();
        let with = import_ventum_reticle(&spec(
            r#",{"type":"circle","x":0,"y":0,"r":2,"start":200,"end":340}"#,
        ))
        .unwrap();
        assert_eq!(
            without, with,
            "adding an arc must not change a single imported mark"
        );
    }

    /// A full ring is decoration under any reading and gets no arc entry — but it is still
    /// counted, under its own tag, so "circle x1" and "arc x1" mean different things.
    #[test]
    fn a_ring_is_counted_as_a_circle_and_an_arc_as_an_arc() {
        let (_, report) = import_ventum_reticle_with_report(
            r#"{"name":"R","unit":"mil","spec":[
                {"type":"dot","x":0,"y":1},
                {"type":"circle","x":0,"y":0,"r":2},
                {"type":"circle","x":0,"y":0,"r":3,"start":0,"end":360},
                {"type":"circle","x":0,"y":0,"r":4,"start":200,"end":340}
            ]}"#,
        )
        .unwrap();
        assert_eq!(
            report.dropped_element_types,
            vec![("arc".to_string(), 1), ("circle".to_string(), 2)],
            "a whole-revolution sweep is a ring, not an arc"
        );
        assert_eq!(report.arcs.len(), 1);
        assert_eq!(report.dropped_elements, 3);
    }

    /// The tally names every type the document drew that produced no hold, and unbound text
    /// is in it: a label that bound to nothing is a lost label, not a lost shape.
    #[test]
    fn the_report_tallies_every_dropped_type_including_unbound_text() {
        let (desc, report) = import_ventum_reticle_with_report(
            r#"{"name":"D","unit":"mil","spec":[
                {"type":"dot","x":0,"y":2},
                {"type":"text","x":0.2,"y":2,"text":"near"},
                {"type":"text","x":9,"y":9,"text":"far"},
                {"type":"line","x1":-5,"y1":0,"x2":5,"y2":0},
                {"type":"line","x1":0,"y1":-5,"x2":0,"y2":5},
                {"type":"rect","x":-1,"y":-1,"w":2,"h":2},
                {"type":"grid","x0":-5,"y0":-5,"x1":5,"y1":5,"step":1},
                {"type":"future_shape","x":0,"y":0}
            ]}"#,
        )
        .unwrap();
        assert_eq!(desc.marks.len(), 1);
        assert_eq!(desc.marks[0].label.as_deref(), Some("near"));
        assert_eq!(
            report.dropped_element_types,
            vec![
                ("grid".to_string(), 1),
                ("line".to_string(), 2),
                ("rect".to_string(), 1),
                ("text".to_string(), 1),
                ("unknown".to_string(), 1),
            ]
        );
        assert_eq!(report.dropped_elements, 6);
        assert_eq!(report.tally(), "grid x1, line x2, rect x1, text x1, unknown x1");
        assert!(!report.is_empty());
    }

    /// A document that imports whole says so: an empty report is the signal that nothing was
    /// lost, and it is what the MBR dot tree produces.
    #[test]
    fn a_fully_representable_reticle_reports_nothing_dropped() {
        let json = format!(r#"{{"name":"MBR","plane":"ffp","unit":"mil","spec":{MBR_SPEC}}}"#);
        let (_, report) = import_ventum_reticle_with_report(&json).unwrap();
        assert!(report.is_empty(), "MBR is all dots: {report:?}");
        assert_eq!(report, VentumImportReport::default());
        assert_eq!(report.tally(), "");
    }

    /// A mirrored pair of horseshoes is two arcs, and the twin's ANGLES are reflected along
    /// with its center. Without that, a mirrored horseshoe would report an apex on the wrong
    /// side of its own arc.
    #[test]
    fn a_mirrored_arc_reflects_its_angles_not_just_its_center() {
        // A horseshoe centered 3 right, opening toward the center line (apex at 0 deg, i.e.
        // to the RIGHT of its own center), mirrored across the vertical axis.
        let (_, report) = import_ventum_reticle_with_report(
            r#"{"name":"M","unit":"mil","spec":[
                {"type":"circle","x":3,"y":0,"r":1,"start":290,"end":70,
                 "repeat":{"axis":"x","step":1,"n":1,"mirror":true}}
            ]}"#,
        )
        .unwrap();
        assert_eq!(report.arcs.len(), 2, "mirror emits both horseshoes");
        let right = report.arcs.iter().find(|a| a.center.right_mil > 0.0).unwrap();
        let left = report.arcs.iter().find(|a| a.center.right_mil < 0.0).unwrap();

        // The original's apex is at 0 deg: 1 mil right of its own center, so 4 mil right.
        assert!((right.apex.right_mil - 4.0).abs() < 1e-9, "{:?}", right.apex);
        // The twin must be the mirror image: 1 mil LEFT of its center at -3, so -4 — not -2,
        // which is what an unreflected sweep would give.
        assert!(
            (left.apex.right_mil + 4.0).abs() < 1e-9,
            "mirrored apex must reflect too (got {:?}, an unreflected sweep gives -2)",
            left.apex
        );
        // Reflection preserves the sweep's extent and keeps both apexes on the same row.
        assert!((left.sweep_degrees - right.sweep_degrees).abs() < 1e-9);
        assert!((left.apex.down_mil - right.apex.down_mil).abs() < 1e-9);
    }

    /// An arc the module cannot resolve (no radius) is still counted, and the discrepancy
    /// between the tally and `arcs.len()` is stated rather than left to be noticed.
    #[test]
    fn an_arc_without_a_radius_is_counted_as_unresolved() {
        let (_, report) = import_ventum_reticle_with_report(
            r#"{"name":"U","unit":"mil","spec":[
                {"type":"circle","x":0,"y":0,"start":200,"end":340},
                {"type":"circle","x":0,"y":0,"r":2,"start":200,"end":340}
            ]}"#,
        )
        .unwrap();
        assert_eq!(report.dropped_element_types, vec![("arc".to_string(), 2)]);
        assert_eq!(report.arcs.len(), 1);
        assert_eq!(report.arcs_unresolved, 1);
        assert_eq!(
            report.arcs.len() + report.arcs_unresolved,
            report.dropped_elements,
            "every counted arc is either resolved or explicitly unresolved"
        );
    }

    /// An MOA arc's center and radius scale to milliradians; its angles do not. Scaling the
    /// angles would rotate the horseshoe by a factor of 3.44.
    #[test]
    fn moa_scales_an_arcs_lengths_but_never_its_angles() {
        let (_, report) = import_ventum_reticle_with_report(
            r#"{"name":"M","unit":"moa","spec":[
                {"type":"circle","x":0,"y":0,"r":2,"start":200,"end":340}
            ]}"#,
        )
        .unwrap();
        let arc = report.arcs[0];
        assert!((arc.radius_mil - 2.0 * MOA_TO_MIL).abs() < 1e-9);
        assert_eq!(arc.start_degrees, 200.0);
        assert_eq!(arc.end_degrees, 340.0);
        assert!((arc.sweep_degrees - 140.0).abs() < 1e-9);
        // Apex still at 270 deg (straight up), now 2 MOA above center.
        assert!((arc.apex.down_mil + 2.0 * MOA_TO_MIL).abs() < 1e-9);
    }

    /// `cx`/`cy` is the schema's other spelling of a circle's center, and a circle naming
    /// neither is centered on the reticle — both must still parse, because they always did.
    #[test]
    fn a_circle_accepts_cx_cy_and_defaults_to_the_center() {
        let (_, report) = import_ventum_reticle_with_report(
            r#"{"name":"C","unit":"mil","spec":[
                {"type":"circle","cx":1,"cy":2,"r":1,"start":0,"end":180},
                {"type":"circle","r":1,"start":0,"end":180}
            ]}"#,
        )
        .unwrap();
        assert_eq!(report.arcs.len(), 2);
        assert_eq!(report.arcs[0].center.right_mil, 1.0);
        assert_eq!(report.arcs[0].center.down_mil, 2.0);
        assert_eq!(report.arcs[1].center, VentumArcPoint::default());
    }

    /// A cosmetic key must never fail the import, whatever type the authoring tool wrote it
    /// as. `r` carrying a unit suffix (`"2mil"`) or a richer object is the realistic case, and
    /// giving `circle` typed fields made it reject the WHOLE document — the dots that have
    /// imported since 0.32.0 would have gone with the ring.
    #[test]
    fn a_non_numeric_radius_is_ignored_rather_than_rejected() {
        for radius in [r#""2mil""#, r#"{"v":2,"unit":"mil"}"#, "[2]", "true", "null"] {
            let json = format!(
                r#"{{"name":"L","unit":"mil","spec":[
                    {{"type":"dot","x":0,"y":1}},
                    {{"type":"circle","x":0,"y":0,"r":{radius}}}
                ]}}"#
            );
            let (desc, report) = import_ventum_reticle_with_report(&json)
                .unwrap_or_else(|e| panic!("r: {radius} must not fail the import, got {e:?}"));

            // The document imported: the sibling dot is still a hold point...
            assert_eq!(desc.marks.len(), 1, "r: {radius} — the dot must survive");
            assert_eq!(desc.marks[0].down_mil, 1.0);
            // ...and the ring is COUNTED as a decoration, not quietly absent.
            assert_eq!(
                report.dropped_element_types,
                vec![("circle".to_string(), 1)],
                "r: {radius} — the ring is counted under its own tag"
            );
            assert_eq!(report.dropped_elements, 1);
        }
    }

    /// The same leniency on the angles. An arc whose `start`/`end` are not numbers describes
    /// no sweep, so it degrades to a plain ring — counted as `circle`, never an error.
    #[test]
    fn non_numeric_angles_degrade_the_arc_to_a_ring_rather_than_failing() {
        let (desc, report) = import_ventum_reticle_with_report(
            r#"{"name":"A","unit":"mil","spec":[
                {"type":"dot","x":0,"y":1},
                {"type":"circle","x":0,"y":0,"r":2,"start":"200deg","end":{"deg":340}}
            ]}"#,
        )
        .unwrap();
        assert_eq!(desc.marks.len(), 1, "the dot must survive");
        assert_eq!(
            report.dropped_element_types,
            vec![("circle".to_string(), 1)],
            "no usable sweep -> a ring, not an arc"
        );
        assert!(report.arcs.is_empty());
        assert_eq!(report.arcs_unresolved, 0);
    }

    /// An arc with real angles but an unreadable radius keeps its `arc` tag and lands in the
    /// stated-discrepancy bucket, exactly like one that omitted `r` altogether.
    #[test]
    fn an_arc_with_an_unreadable_radius_is_counted_as_unresolved() {
        let (_, report) = import_ventum_reticle_with_report(
            r#"{"name":"U","unit":"mil","spec":[
                {"type":"circle","x":0,"y":0,"r":"2mil","start":200,"end":340}
            ]}"#,
        )
        .unwrap();
        assert_eq!(report.dropped_element_types, vec![("arc".to_string(), 1)]);
        assert_eq!(report.arcs_unresolved, 1);
        assert!(report.arcs.is_empty());
        assert_eq!(report.arcs.len() + report.arcs_unresolved, report.dropped_elements);
    }

    /// A non-numeric center falls back to the reticle center, the same answer an absent
    /// `x`/`y` has always given, and the arc still resolves around it.
    #[test]
    fn a_non_numeric_center_falls_back_to_the_reticle_center() {
        let (_, report) = import_ventum_reticle_with_report(
            r#"{"name":"C","unit":"mil","spec":[
                {"type":"circle","x":"left","y":null,"r":2,"start":200,"end":340}
            ]}"#,
        )
        .unwrap();
        assert_eq!(report.arcs.len(), 1);
        assert_eq!(report.arcs[0].center, VentumArcPoint::default());
    }

    /// Leniency must not cost a real arc its geometry: a genuine numeric circle in the SAME
    /// document as a cosmetic one still resolves to a drawn arc, apex and tips included.
    #[test]
    fn a_numeric_arc_still_resolves_beside_a_cosmetic_one() {
        let (desc, report) = import_ventum_reticle_with_report(
            r#"{"name":"H","unit":"mil","spec":[
                {"type":"dot","x":0,"y":1},
                {"type":"circle","x":0,"y":0,"r":"2mil"},
                {"type":"circle","x":0,"y":0,"r":2,"start":200,"end":340}
            ]}"#,
        )
        .unwrap();
        assert_eq!(desc.marks.len(), 1, "the dot is unaffected by either circle");
        assert_eq!(
            report.dropped_element_types,
            vec![("arc".to_string(), 1), ("circle".to_string(), 1)],
            "the cosmetic ring and the real horseshoe are counted apart"
        );
        assert_eq!(report.dropped_elements, 2);

        // The real horseshoe is fully resolved: radius, sweep, and an apex straight up.
        assert_eq!(report.arcs.len(), 1);
        assert_eq!(report.arcs_unresolved, 0);
        let arc = report.arcs[0];
        assert_eq!(arc.radius_mil, 2.0);
        assert!((arc.sweep_degrees - 140.0).abs() < 1e-9);
        assert!(
            (arc.apex.down_mil + 2.0).abs() < 1e-9 && arc.apex.right_mil.abs() < 1e-9,
            "apex at 270 deg is 2 mil straight above center, got {:?}",
            arc.apex
        );
    }

    /// A CENTERED asymmetric arc with `mirror:true` draws two different shapes, because
    /// reflection reworks the angles as well as the center. Both must be counted, and their
    /// apexes must land on opposite sides — reporting one would be the undercount MBA-1441
    /// exists to prevent.
    #[test]
    fn a_centered_asymmetric_arc_still_gets_its_mirrored_twin() {
        let (_, report) = import_ventum_reticle_with_report(
            r#"{"name":"A","unit":"mil","spec":[
                {"type":"circle","x":0,"y":0,"r":2,"start":290,"end":70,
                 "repeat":{"axis":"x","step":1,"n":1,"mirror":true}}
            ]}"#,
        )
        .unwrap();
        assert_eq!(report.dropped_element_types, vec![("arc".to_string(), 2)]);
        assert_eq!(report.arcs.len(), 2, "the centered twin is a second shape");

        // Original apex is at 290 + 140/2 = 0 deg, i.e. 2 mil RIGHT of center; the twin's is
        // at 180 deg, 2 mil LEFT. Same center, same sweep, opposite sides.
        let (first, second) = (report.arcs[0], report.arcs[1]);
        assert!((first.apex.right_mil - 2.0).abs() < 1e-9, "{:?}", first.apex);
        assert!((second.apex.right_mil + 2.0).abs() < 1e-9, "{:?}", second.apex);
        assert!((first.sweep_degrees - second.sweep_degrees).abs() < 1e-9);
        assert_eq!(first.center, second.center);
    }

    /// The dedupe is not abandoned, only narrowed to what it was always for. An arc that the
    /// reflection maps onto ITSELF — a horseshoe centered on the axis it is mirrored across,
    /// or any ring — is still one shape and is still counted once.
    #[test]
    fn a_centered_symmetric_arc_is_still_deduped() {
        // start 200 + end 340 = 540 = 180 (mod 360): symmetric about the vertical axis.
        let (_, horseshoe) = import_ventum_reticle_with_report(
            r#"{"name":"S","unit":"mil","spec":[
                {"type":"circle","x":0,"y":0,"r":2,"start":200,"end":340,
                 "repeat":{"axis":"x","step":1,"n":1,"mirror":true}}
            ]}"#,
        )
        .unwrap();
        assert_eq!(horseshoe.dropped_element_types, vec![("arc".to_string(), 1)]);
        assert_eq!(horseshoe.arcs.len(), 1);

        // ...and a plain ring maps onto itself under any reflection through its own center.
        let (_, ring) = import_ventum_reticle_with_report(
            r#"{"name":"R","unit":"mil","spec":[
                {"type":"circle","x":0,"y":0,"r":2,
                 "repeat":{"axis":"y","step":1,"n":1,"mirror":true}}
            ]}"#,
        )
        .unwrap();
        assert_eq!(ring.dropped_element_types, vec![("circle".to_string(), 1)]);

        // The SAME horseshoe mirrored across the HORIZONTAL axis is not symmetric about it —
        // it opens the other way — so that one is two shapes.
        let (_, flipped) = import_ventum_reticle_with_report(
            r#"{"name":"F","unit":"mil","spec":[
                {"type":"circle","x":0,"y":0,"r":2,"start":200,"end":340,
                 "repeat":{"axis":"y","step":1,"n":1,"mirror":true}}
            ]}"#,
        )
        .unwrap();
        assert_eq!(flipped.dropped_element_types, vec![("arc".to_string(), 2)]);
        assert!((flipped.arcs[0].apex.down_mil + 2.0).abs() < 1e-9, "apex up");
        assert!((flipped.arcs[1].apex.down_mil - 2.0).abs() < 1e-9, "twin apex down");
    }

    #[test]
    fn imported_mbr_reticle_hold_solves() {
        let json = format!(r#"{{"name":"MBR","plane":"ffp","unit":"mil","spec":{MBR_SPEC}}}"#);
        let desc = import_ventum_reticle(&json).unwrap();
        // (down 4, right 2) is a real MBR mark, so a matching hold lands on it.
        let hold = hold_point_in_reticle(4.0, 2.0, 1.0, &desc).unwrap();
        assert!(hold.nearest_mark.is_some());
        assert!(
            hold.nearest_mark_distance_mil < 0.5,
            "hold should sit on/near a real mark (distance {})",
            hold.nearest_mark_distance_mil
        );
        assert!(!hold.off_reticle, "a mid-tree hold is on the reticle");
    }
}
