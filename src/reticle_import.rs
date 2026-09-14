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
//! So an arc is still not a mark — but it is no longer a silence either. Each one whose
//! geometry this module can read is resolved (using exactly the convention above, which is why
//! the convention was recorded) into the three points a shooter could actually index on, and
//! handed back on the report as a [`VentumArc`]: its two tips and its apex, in the engine's own
//! `right_mil` / `down_mil`. One whose radius, sweep or center the document wrote in a form
//! this module cannot read is counted in [`VentumImportReport::arcs_unresolved`] instead —
//! reporting points computed from a fallback value nobody wrote would be inventing the very
//! holds the paragraph above refuses to invent. A caller who knows their horseshoe is
//! hold-bearing turns the resolved ones into marks with three [`ReticleMark::new`] calls and
//! never re-derives the clockwise/`+y`-down trap; a caller who does not, at least learns the
//! document drew something they cannot aim with.
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
/// `arc`, `circle` and `text` entries count *expanded* instances. This report resolves each
/// arc instance's own geometry, and a mirrored pair of horseshoes is genuinely two of them; an
/// unbound `text` is a lost LABEL rather than a lost shape, so a three-copy ladder that binds
/// to nothing has lost three of them. The one exception is a `circle` whose `repeat` this
/// module could not read: it is expanded as though it had none and so counts once, and
/// [`Self::circle_repeats_unreadable`] says how many elements that happened to, so the
/// shortfall is declared rather than hidden.
///
/// `line`, `rect`, `grid` and unknown types count elements *as the document writes them* — a
/// `repeat` on one of those counts once, since expansion reads no geometry from them and the
/// report carries none to tell a dropped shape's copies from the shape itself.
///
/// Pinned by `tests::a_repeat_counts_per_instance_for_text_and_once_for_a_dropped_shape`.
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
    /// Arcs the document declared but this module could not resolve into points.
    ///
    /// An arc's apex and tips are built from three things — a center, a radius that is a
    /// positive length, and a sweep — and the importer computes nothing without
    /// all three. Anything that leaves one of them unavailable lands the arc here: a key the
    /// document omitted, a value spelled in a type this module cannot read, or a radius that
    /// reads perfectly well as a number and is not a length (`0`, or negative). Each of those
    /// cases is driven through the importer by
    /// `tests::an_arc_resolves_only_from_a_center_a_positive_radius_and_a_sweep`.
    ///
    /// They are the SAME loss and land here for the same reason: an arc is reported so a
    /// shooter can copy its apex and tips onto their reticle, and a coordinate derived from a
    /// value nobody wrote is worse than no coordinate at all. They are counted in
    /// [`Self::dropped_element_types`] under [`ARC_TAG`] like any other arc, so `arcs.len()`
    /// plus this equals that tally — the discrepancy is stated rather than left to be
    /// noticed.
    pub arcs_unresolved: usize,
    /// How many `circle` elements declared a `repeat` this module could not read, so each was
    /// expanded as though it had none and counted ONCE where the document may draw many.
    ///
    /// A circle's repeat is read leniently (see the serde section below) because a circle is
    /// never a hold point and must not be able to fail an import. That leniency costs the
    /// tally its accuracy for those elements, and this is the accounting for it: the count is
    /// low by an unknown amount, and the caller is told so instead of being handed a total
    /// that quietly is not one.
    pub circle_repeats_unreadable: usize,
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
    let mut circle_repeats_unreadable = 0usize;
    let instances = expand_elements(
        &reticle.spec.0,
        scale,
        &mut dropped,
        &mut circle_repeats_unreadable,
    )?;

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
                let tag = if circle.is_arc() { ARC_TAG } else { "circle" };
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
        circle_repeats_unreadable,
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
    /// Radius in milliradians, or `None` when the element omitted `r` or wrote one this
    /// module could not read as a number. A number it COULD read is kept exactly as written,
    /// `0` and negatives included; [`Self::resolve`] is the one place a radius is judged as a
    /// length. Either way the element is still counted; it just cannot be resolved.
    radius_mil: Option<f64>,
    /// The document's `start`/`end` angles in degrees, when it declared both AS NUMBERS.
    /// `None` is a plain ring, which is decoration under any reading and gets no [`VentumArc`]
    /// — unless [`Self::sweep_unreadable`] says the document declared a sweep this module
    /// could not read, in which case it is an arc whose angles are lost, not a ring.
    /// Both members are finite by construction ([`Cosmetic::Value`]).
    angles: Option<(f64, f64)>,
    /// The document gave BOTH `start` and `end` a value, and at least one of those values is
    /// not a number.
    ///
    /// This is the difference between "no sweep was declared" and "a sweep was declared and I
    /// cannot read it", and collapsing the two is how a horseshoe turns into a ring in
    /// silence. One angle alone still describes no sweep — the format needs both — so a
    /// circle that declares only `start` is a ring whatever type that `start` has.
    ///
    /// "Gave a value" excludes `null`, which this model reads as absence everywhere
    /// ([`Cosmetic::Absent`]). So `{"start": null, "end": "340deg"}` is a ring too: the
    /// document supplied one angle, not two, and one angle is no sweep however it is spelled.
    /// Pinned by `tests::null_is_read_as_absence_not_as_a_bad_value`.
    sweep_unreadable: bool,
    /// The document wrote `x`/`cx` or `y`/`cy` as something that is not a number, so the
    /// center this shape is resolved about is the fallback origin rather than the one the
    /// document meant. On an arc it would put the reported apex and tips somewhere nobody
    /// wrote, so it makes the arc unresolved. A ring is never resolved and so loses no points
    /// to it — but it does lose the mirror dedupe, which cannot prove a twin is a duplicate
    /// about a center it does not have, and is then counted as the two shapes it drew.
    center_unreadable: bool,
}

impl CircleShape {
    /// Whether this element is an ARC — a `circle` the document gave a sweep — rather than a
    /// plain ring. True both for a sweep this module read and for one it could not, because
    /// the tag a document earns must not depend on whether its angles happened to be spelled
    /// in a type this module understands.
    fn is_arc(&self) -> bool {
        self.sweep_degrees().is_some() || self.sweep_unreadable
    }

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
    /// also reworks the angles ([`Role::mirrored`]). Given geometry this module can read — the
    /// body below is where unreadable geometry is turned down, since it proves nothing either
    /// way — a ring maps onto itself under any reflection through its own center, and a swept
    /// arc does so exactly when its two endpoints
    /// swap into each other: reflection about the vertical axis sends `theta -> 180 - theta`,
    /// so it needs `start + end == 180`; about the horizontal axis it sends `theta -> -theta`,
    /// so it needs `start + end == 0` — both modulo a full revolution. A horseshoe centered on
    /// the vertical axis (`start: 200, end: 340`, sum 540 ≡ 180) is symmetric about it and is
    /// deduped; the same horseshoe mirrored vertically opens the other way and is not.
    fn is_reflection_of_itself(&self, axis: Axis) -> bool {
        // Geometry the document wrote and this module could not read answers NOTHING here.
        // The dedupe drops a twin only when it can prove the twin is a duplicate, and with
        // the angles or the center unreadable there is no proof either way — so the twin is
        // kept and the document is credited with the two shapes it drew. (On an arc both
        // copies then land in `arcs_unresolved`, where the caller is told the geometry was
        // lost; a ring is only ever counted, so for one this is the whole of the difference.)
        // Guessing "duplicate" here would undercount in silence, which is the failure this
        // whole report exists to end.
        if self.sweep_unreadable || self.center_unreadable {
            return false;
        }
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
            sweep_unreadable: self.sweep_unreadable,
            center_unreadable: self.center_unreadable,
        }
    }

    /// Resolve this circle, centered at `(down, right)` milliradians, into the arc points a
    /// caller could adopt as holds. `None` when it is a ring, or when the element declared
    /// geometry this module could not read — the caller counts those separately rather than
    /// reporting points derived from a fallback the document never wrote.
    fn resolve(&self, down: f64, right: f64) -> Option<VentumArc> {
        // A center the document declared and this module could not read is exactly as fatal
        // to the answer as a missing radius: both feed the same `point()` below, and a
        // coordinate computed from a value nobody wrote is a hold this report would be
        // inventing. That half of this guard is the one that does the work.
        //
        // The `sweep_unreadable` half decides nothing on its own — an unreadable angle leaves
        // `angles` empty, so `sweep_degrees()?` on the next line would bail regardless
        // (removing it fails no test; that was checked). It is stated anyway because the two
        // flags mean the same thing to a reader, and because `is_arc` above genuinely depends
        // on the distinction.
        if self.sweep_unreadable || self.center_unreadable {
            return None;
        }
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
///
/// `circle_repeats_unreadable` counts the circles whose `repeat` this module could not read
/// and therefore expanded as though absent, so the caller learns the tally is low rather than
/// being handed a total that quietly is not one.
fn expand_elements(
    elements: &[Element],
    scale: f64,
    dropped: &mut BTreeMap<&'static str, usize>,
    circle_repeats_unreadable: &mut usize,
) -> Result<Vec<ExpandedInstance>, ReticleError> {
    let mut out: Vec<ExpandedInstance> = Vec::new();
    for element in elements {
        let (x, y, repeat, base) = match element {
            Element::Dot { x, y, repeat } => (*x, *y, repeat.as_ref(), PointRole::Dot),
            Element::Tick { x, y, repeat } => (*x, *y, repeat.as_ref(), PointRole::Tick),
            Element::Text {
                x, y, text, repeat,
            } => (*x, *y, repeat.as_ref(), PointRole::Text(text)),
            Element::Circle(CircleFields {
                x,
                y,
                r,
                start,
                end,
                repeat,
            }) => {
                if repeat.is_unreadable() {
                    *circle_repeats_unreadable += 1;
                }
                (
                    x.value().unwrap_or(0.0),
                    y.value().unwrap_or(0.0),
                    repeat.as_repeat(),
                    PointRole::Circle(CircleShape {
                        // The radius is a length in the reticle's unit, so it scales with the
                        // coordinates; the angles are angles and do not.
                        radius_mil: r.value().map(|r| r * scale),
                        angles: start.value().zip(end.value()),
                        // The format needs BOTH angles for a sweep, so "the document declared
                        // an arc" means both keys are there; it is unreadable when either of
                        // them is something this module cannot read as a number.
                        sweep_unreadable: start.is_declared()
                            && end.is_declared()
                            && (start.is_unreadable() || end.is_unreadable()),
                        center_unreadable: x.is_unreadable() || y.is_unreadable(),
                    }),
                )
            }
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
// cosmetic keys (color, width, size, len, orient, max_extent, tube_diameter, notes,
// manufacturer, ...) are ignored rather than rejected. A `circle` goes further, because it
// is the one decoration this module reads geometry from and reading it must not cost a
// document its dots: no key of a `circle` refuses the document over how it is written — not
// for its type (`"r": "2mil"`), not for which of the schema's two spellings of a center it
// uses (`x` beside `cx`), and not for being written twice. See `CircleFields`.
//
// Among the drawing ELEMENTS, strict is what a hold is built from: a `dot`'s or `tick`'s
// `x`/`y`, a `text`'s `x`/`y` and its string, and the `repeat` that stamps copies of any of
// those. There is no sane fallback for a mark whose position cannot be read, and a mark's
// `repeat` quietly degrading to a single copy would drop hold points without saying so. That
// boundary is stated as code and swept in both directions by
// `tests::strictness_is_exactly_what_a_hold_is_built_from`, which drives every element type
// this model can produce against every key this importer reads (plus two it ignores). A
// field that starts refusing a document, or stops, turns it red rather than leaving this
// paragraph wrong.
//
// Reticle-level metadata (`name`, `plane`, `unit`, `ref_magnification`) is strict as well,
// `name` included even though no hold depends on it. That is unchanged since 0.32.0 and is
// not what this leniency is about: nothing there is per-element, so a bad value is the
// document saying something wrong about itself rather than one decoration spoiling the rest.
//
// What a `circle` can still cost a document, neither of them a key it carries:
//
//   * The mark cap. `circle` instances now go through repeat expansion, which enforces
//     MAX_RETICLE_MARKS over everything it materializes, so a document at the cap is refused
//     where 0.32.0 dropped its rings before expansion began. This is the one way this branch
//     found to refuse a document 0.32.0 imported; it is a cap decision rather than a parsing
//     one, is pinned by `tests::a_circle_counts_against_the_mark_cap_unlike_0_32_0`, and is
//     still owed a fix.
//   * A refusal that lands before this model is consulted at all: serde_json's PARSER on a
//     numeric literal outside `f64`'s range (`"r": 1e400`) or JSON nested past the recursion
//     limit, and serde's tag reader on an element that writes `type` twice. Nothing here can
//     see any of them, and neither could 0.32.0. They are not about the circle —
//     `tests::a_refusal_that_is_not_about_a_key_falls_on_every_element_type_alike` shows them
//     falling identically on element types that read no fields whatsoever.
//
// That list is what a corpus of hostile documents turned up, not a proof that nothing else
// exists — the sweeps named above are where one that came from a key would show itself.
//
// Leniency is never silence. A circle that loses geometry this way is still counted in the
// report, and an ARC that loses any of it lands in `arcs_unresolved` rather than being
// resolved from a fallback; a circle that loses its `repeat` is counted in
// `circle_repeats_unreadable`, because the tally for that element is then low.
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

/// A number this module treats as cosmetic geometry, exactly as the document wrote it.
///
/// The three cases are kept apart because leniency degrades a field to its DEFAULT, and a
/// default is only an honest answer when the document did not state one. A `circle` with no
/// `x` really is centered; a `circle` whose `x` is `"left"` is centered only because this
/// module gave up reading it, and an arc resolved about that center would hand the caller an
/// apex nobody wrote. So [`Self::Unreadable`] is remembered rather than folded into
/// [`Self::Absent`], and the report says a shape was lost instead of inventing one.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
enum Cosmetic {
    /// The document did not write this key — or wrote `null`, which is how JSON spells the
    /// absence of a value.
    ///
    /// The other place this model accepts a `null` is a mark's optional `repeat`, where
    /// serde's `Option` reads it the same way. It is NOT how the reticle-level fields behave:
    /// `#[serde(default)]` fills in a key the document omitted and still refuses a present
    /// `null`, which is 0.32.0's behaviour and is deliberately left alone. Both halves are
    /// pinned by `tests::null_is_read_as_absence_not_as_a_bad_value`.
    #[default]
    Absent,
    /// The key is there, holding something this module cannot read as a number.
    Unreadable,
    /// A usable number. Finite by construction — see [`Cosmetic::from_json`].
    Value(f64),
}

impl Cosmetic {
    /// The number, when there is a usable one.
    fn value(self) -> Option<f64> {
        match self {
            Cosmetic::Value(number) => Some(number),
            Cosmetic::Absent | Cosmetic::Unreadable => None,
        }
    }

    /// Whether the document wrote this key at all, readable or not.
    fn is_declared(self) -> bool {
        !matches!(self, Cosmetic::Absent)
    }

    /// Whether the document wrote this key and this module could not read it.
    fn is_unreadable(self) -> bool {
        matches!(self, Cosmetic::Unreadable)
    }

    /// Classify a cosmetic number the document wrote, never failing the document over one
    /// this module cannot use.
    ///
    /// A typed field — `Option<f64>` included — absorbs an ABSENT or `null` value but still
    /// *rejects* a present one of the wrong type, and one rejected field fails the entire
    /// document, not the element that carried it. That is the wrong trade for a key this
    /// module only ever reads as a nicety: a Ventum tool is free to write `"r": "2mil"`,
    /// `"r": {"v": 2, "unit": "mil"}` or anything else its UI finds convenient, and none of
    /// that is a reason to refuse a reticle whose dots are perfectly good. So the value
    /// arrives as arbitrary JSON and is classified here, never rejected.
    ///
    /// What the classification can and cannot see: `Value::as_f64` answers `None` for every
    /// JSON value that is not a number, which is the case this exists for. The `is_finite`
    /// guard behind it cannot fire on JSON at all — `NaN` and `Infinity` are not JSON
    /// literals, and a numeric literal outside `f64`'s range (`1e400`) is refused by
    /// serde_json's PARSER, which fails the whole document before any [`serde_json::Value`]
    /// exists for this function to see. The guard is kept only so a non-JSON deserializer
    /// could not slip a non-finite number past; it catches nothing a Ventum file can contain.
    fn from_json(value: serde_json::Value) -> Cosmetic {
        if value.is_null() {
            return Cosmetic::Absent;
        }
        match value.as_f64().filter(|v| v.is_finite()) {
            Some(number) => Cosmetic::Value(number),
            None => Cosmetic::Unreadable,
        }
    }
}

/// A `circle`'s `repeat`, exactly as the document wrote it. See [`CosmeticRepeat::from_json`].
#[derive(Debug, Default)]
enum CosmeticRepeat {
    /// The document did not write a `repeat` (or wrote `null`).
    #[default]
    Absent,
    /// It wrote one this module cannot read — precisely, one [`Repeat`] itself refuses:
    /// missing `axis`, `step` or `n`, any of them of the wrong type or out of range, an
    /// `axis` that is neither `x` nor `y`, a `repeat` that is not even an object.
    Unreadable,
    /// A usable repeat.
    Present(Repeat),
}

impl CosmeticRepeat {
    /// The repeat to expand with, `None` when there is none to use.
    fn as_repeat(&self) -> Option<&Repeat> {
        match self {
            CosmeticRepeat::Present(repeat) => Some(repeat),
            CosmeticRepeat::Absent | CosmeticRepeat::Unreadable => None,
        }
    }

    /// Whether the document declared a repeat this module could not read.
    fn is_unreadable(&self) -> bool {
        matches!(self, CosmeticRepeat::Unreadable)
    }

    /// Classify a `circle`'s `repeat`, degrading an unusable one to "no repeat" instead of
    /// failing the document over it.
    ///
    /// This leniency belongs to the CIRCLE, not to [`Repeat`], and the two must not be
    /// confused. On a `dot`, `tick` or `text` a repeat is not decoration: it stamps the hold
    /// points the shooter aims with, so a malformed one that quietly collapsed to a single
    /// copy would drop marks without saying so — the exact silence MBA-1441 exists to end,
    /// and loud is the only safe answer there. On a `circle` it stamps copies of something
    /// that is never a mark under any reading, so an unreadable repeat costs the report a
    /// count (declared in [`VentumImportReport::circle_repeats_unreadable`]) and costs the
    /// reticle nothing. [`Repeat`] itself therefore keeps its strict derive; this is the one
    /// place a value is offered to it and a failure turned into a shrug.
    ///
    /// The degradation is all-or-nothing on purpose. `axis`, `step` and `n` have no
    /// defensible defaults — `x` is not more plausible than `y`, and no copy count is more
    /// plausible than another — so a repeat with any of them unusable is not a repeat. "No
    /// repeat", one instance at `(x, y)`, is what an absent `repeat` gives, and is what
    /// 0.32.0 gave every `circle` in every document, since the variant read no fields at all.
    fn from_json(value: serde_json::Value) -> CosmeticRepeat {
        if value.is_null() {
            return CosmeticRepeat::Absent;
        }
        match serde_json::from_value::<Repeat>(value) {
            Ok(repeat) => CosmeticRepeat::Present(repeat),
            Err(_) => CosmeticRepeat::Unreadable,
        }
    }
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
    /// either way: still dropped, but no longer in silence (MBA-1441). Its keys are read by
    /// [`CircleFields`], which refuses the document over none of them.
    Circle(CircleFields),
    /// Decoration — dropped.
    Rect {},
    /// Decoration — dropped.
    Grid {},
    /// Any unknown future element type — dropped.
    #[serde(other)]
    Unknown,
}

/// The keys of a `circle`, read so that none of them can refuse the document.
///
/// A derived `Deserialize` gives a struct three ways to reject an element, and a `circle` is
/// decoration this module must be able to skip rather than die on, so all three are answered
/// here instead:
///
/// * **an unknown key** — already harmless, since no variant in this model sets
///   `deny_unknown_fields`; the `_` arm below keeps it that way,
/// * **a key whose value is the wrong type** — [`Cosmetic::from_json`] and
///   [`CosmeticRepeat::from_json`] classify a value rather than typing it, so `"r": "2mil"`
///   degrades to "unreadable" instead of failing the parse,
/// * **a key written twice** — which is what a derived struct calls `duplicate field \`x\``
///   and refuses the whole document over. That is the one this hand-written visitor exists
///   for. It bites two ways a Ventum file really can be written: the same key repeated
///   (`{"r": 1, "r": 2}`, which JSON permits and `serde_json::Value` itself resolves by
///   keeping the last), and the schema's two spellings of one center used together
///   (`{"x": 1, "cx": 2}`, where the derive's `alias` collapses both onto `x` and then calls
///   the second one a duplicate). 0.32.0 accepted all of those, because its `Circle {}` read
///   no keys at all; a reticle whose dots are perfectly good must not start failing over how
///   its rings spell a center.
///
/// The rule is last-one-wins, in document order — the same answer `serde_json::Value` gives a
/// repeated key, applied to the alias pairs as well. A mark's fields get none of this: see
/// the serde section above for why a `dot`'s `x` stays strict.
#[derive(Debug, Default)]
struct CircleFields {
    /// Center, `+x` right. `cx` is the schema's other spelling of the same field.
    x: Cosmetic,
    /// Center, `+y` down. `cy` is the other spelling.
    y: Cosmetic,
    /// Radius, in the reticle's own unit.
    r: Cosmetic,
    /// Sweep start in degrees from 3 o'clock, clockwise (see the module documentation — 270
    /// is the TOP of the reticle). Present only on an arc.
    start: Cosmetic,
    /// Sweep end, same convention. An arc needs BOTH: one angle alone describes no sweep.
    end: Cosmetic,
    /// Copies of the ring or arc.
    repeat: CosmeticRepeat,
}

impl<'de> Deserialize<'de> for CircleFields {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct CircleFieldsVisitor;

        impl<'de> Visitor<'de> for CircleFieldsVisitor {
            type Value = CircleFields;

            fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str("a circle element")
            }

            fn visit_map<A>(self, mut map: A) -> Result<CircleFields, A::Error>
            where
                A: de::MapAccess<'de>,
            {
                let mut fields = CircleFields::default();
                while let Some(key) = map.next_key::<String>()? {
                    match key.as_str() {
                        "x" | "cx" => fields.x = Cosmetic::from_json(map.next_value()?),
                        "y" | "cy" => fields.y = Cosmetic::from_json(map.next_value()?),
                        "r" => fields.r = Cosmetic::from_json(map.next_value()?),
                        "start" => fields.start = Cosmetic::from_json(map.next_value()?),
                        "end" => fields.end = Cosmetic::from_json(map.next_value()?),
                        "repeat" => fields.repeat = CosmeticRepeat::from_json(map.next_value()?),
                        // Every other key — `color`, `width`, `max_extent`, the next thing
                        // the authoring tool invents — is decoration on decoration. Skipped
                        // without materializing it.
                        _ => {
                            map.next_value::<de::IgnoredAny>()?;
                        }
                    }
                }
                Ok(fields)
            }
        }

        deserializer.deserialize_map(CircleFieldsVisitor)
    }
}

/// The `repeat` operator: stamp `n` copies of an element along one axis, optionally
/// mirrored, optionally auto-numbering a text ladder's labels.
///
/// Strictly typed, and deliberately so: on a `dot`, `tick` or `text` this is what stamps the
/// hold points, and a malformed repeat that degraded to one copy would drop marks in silence.
/// A `circle`'s repeat stamps nothing holdable, so that one field — and only that one — is
/// offered to this type through [`CosmeticRepeat::from_json`], which keeps the failure
/// instead of propagating it.
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
        let mut circle_repeats_unreadable = 0usize;
        let instances =
            expand_elements(&elements, 1.0, &mut dropped, &mut circle_repeats_unreadable).unwrap();
        assert!(dropped.is_empty(), "a text ladder drops nothing");
        assert_eq!(circle_repeats_unreadable, 0);
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

    /// The same leniency on the angles — with the same SAYING SO. An arc whose `start`/`end`
    /// are not numbers has lost its sweep, and that is the same kind of loss as an unreadable
    /// radius: the document declared an arc this module cannot resolve. Calling it a plain
    /// ring would drop the declaration and report nothing, which is exactly the silence
    /// MBA-1441 is named for.
    #[test]
    fn a_declared_sweep_that_cannot_be_read_is_reported_not_quietly_a_ring() {
        for (start, end) in [
            (r#""200deg""#, r#"{"deg":340}"#),
            ("200", r#""340deg""#),
            (r#""200deg""#, "340"),
        ] {
            let json = format!(
                r#"{{"name":"A","unit":"mil","spec":[
                    {{"type":"dot","x":0,"y":1}},
                    {{"type":"circle","x":0,"y":0,"r":2,"start":{start},"end":{end}}}
                ]}}"#
            );
            let (desc, report) = import_ventum_reticle_with_report(&json)
                .unwrap_or_else(|e| panic!("{start}/{end} must not fail the import, got {e:?}"));
            assert_eq!(desc.marks.len(), 1, "{start}/{end} — the dot must survive");
            assert_eq!(
                report.dropped_element_types,
                vec![("arc".to_string(), 1)],
                "{start}/{end} — a declared sweep keeps the arc tag even when unreadable"
            );
            assert!(report.arcs.is_empty(), "{start}/{end} — nothing to resolve");
            assert_eq!(
                report.arcs_unresolved, 1,
                "{start}/{end} — an unreadable sweep is reported exactly like an unreadable \
                 radius, not swallowed"
            );
            assert_eq!(
                report.arcs.len() + report.arcs_unresolved,
                report.dropped_elements,
                "{start}/{end} — the stated arcs-vs-tally invariant still holds"
            );
        }
    }

    /// One angle alone still describes no sweep, whatever type it is, so a `circle` carrying
    /// only `start` is a ring — and remains one when that lone angle is unreadable. This is
    /// the boundary of the rule above: an arc needs BOTH angles declared before there is a
    /// declaration to lose.
    #[test]
    fn a_lone_angle_is_a_ring_whether_or_not_it_can_be_read() {
        for spec in [
            r#"{"type":"circle","x":0,"y":0,"r":2,"start":200}"#,
            r#"{"type":"circle","x":0,"y":0,"r":2,"start":"200deg"}"#,
            r#"{"type":"circle","x":0,"y":0,"r":2,"end":340}"#,
        ] {
            let json = format!(r#"{{"name":"R","unit":"mil","spec":[{spec}]}}"#);
            let (_, report) = import_ventum_reticle_with_report(&json).unwrap();
            assert_eq!(
                report.dropped_element_types,
                vec![("circle".to_string(), 1)],
                "{spec} — one angle is no sweep"
            );
            assert_eq!(report.arcs_unresolved, 0, "{spec}");
        }
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

    /// An ABSENT center — omitted, or `null`, which is how JSON spells absent — is the
    /// reticle center, and the arc resolves around it. That is the format's own default, not
    /// a guess.
    #[test]
    fn an_absent_or_null_center_is_the_reticle_center() {
        let (_, report) = import_ventum_reticle_with_report(
            r#"{"name":"C","unit":"mil","spec":[
                {"type":"circle","r":2,"start":200,"end":340},
                {"type":"circle","x":null,"y":null,"r":2,"start":200,"end":340}
            ]}"#,
        )
        .unwrap();
        assert_eq!(report.arcs.len(), 2);
        assert_eq!(report.arcs_unresolved, 0);
        assert_eq!(report.arcs[0].center, VentumArcPoint::default());
        assert_eq!(report.arcs[1].center, VentumArcPoint::default());
    }

    /// A center the document WROTE and this module cannot read is a different thing entirely.
    /// Resolving the arc about the origin anyway would print an apex and two tips at
    /// coordinates nobody wrote, under a notice that invites the reader to add them as marks
    /// — inventing a hold, which is the other half of MBA-1441's complaint. So it is
    /// unresolved, reported beside an unreadable radius and an unreadable sweep.
    #[test]
    fn a_center_that_cannot_be_read_leaves_the_arc_unresolved() {
        for (x, y) in [
            (r#""left""#, "0"),
            ("0", r#""up""#),
            (r#"{"v":1}"#, r#"{"v":1}"#),
        ] {
            let json = format!(
                r#"{{"name":"C","unit":"mil","spec":[
                    {{"type":"circle","x":{x},"y":{y},"r":2,"start":200,"end":340}}
                ]}}"#
            );
            let (_, report) = import_ventum_reticle_with_report(&json)
                .unwrap_or_else(|e| panic!("x:{x} y:{y} must not fail the import, got {e:?}"));
            assert_eq!(
                report.dropped_element_types,
                vec![("arc".to_string(), 1)],
                "x:{x} y:{y} — still an arc"
            );
            assert!(
                report.arcs.is_empty(),
                "x:{x} y:{y} — no points may be reported from a center nobody wrote"
            );
            assert_eq!(report.arcs_unresolved, 1, "x:{x} y:{y}");
        }
    }

    /// A ring is never resolved, so an unreadable center costs it nothing: it is counted as
    /// the plain `circle` it is, and no arc bucket moves.
    #[test]
    fn an_unreadable_center_on_a_plain_ring_is_just_a_counted_ring() {
        let (_, report) = import_ventum_reticle_with_report(
            r#"{"name":"C","unit":"mil","spec":[{"type":"circle","x":"left","r":2}]}"#,
        )
        .unwrap();
        assert_eq!(report.dropped_element_types, vec![("circle".to_string(), 1)]);
        assert_eq!(report.arcs_unresolved, 0);
        assert!(report.arcs.is_empty());

        // It is not entirely without consequence, though, and the consequence is a count.
        // Mirroring a ring across the axis it sits on normally dedupes the twin away; an
        // unreadable center makes the stepped coordinate a fallback, which proves nothing, so
        // the twin is kept and the document is credited with the two shapes it drew.
        let (_, mirrored) = import_ventum_reticle_with_report(
            r#"{"name":"C","unit":"mil","spec":[
                {"type":"circle","x":"left","y":0,"r":2,
                 "repeat":{"axis":"x","step":0,"n":1,"mirror":true}}
            ]}"#,
        )
        .unwrap();
        assert_eq!(
            mirrored.dropped_element_types,
            vec![("circle".to_string(), 2)],
            "an unprovable dedupe keeps the twin, ring or not"
        );

        // The same ring with a center this module CAN read is one shape, as it always was.
        let (_, readable) = import_ventum_reticle_with_report(
            r#"{"name":"C","unit":"mil","spec":[
                {"type":"circle","x":0,"y":0,"r":2,
                 "repeat":{"axis":"x","step":0,"n":1,"mirror":true}}
            ]}"#,
        )
        .unwrap();
        assert_eq!(readable.dropped_element_types, vec![("circle".to_string(), 1)]);
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

    /// A `circle`'s `repeat` is the sixth cosmetic field and is read like the other five: a
    /// repeat this module cannot use is the same answer as no repeat, which is what 0.32.0
    /// gave every circle in every document. Each of these spellings imported under 0.32.0 —
    /// `Circle {}` read no fields at all — so each must import now.
    #[test]
    fn an_unusable_repeat_on_a_circle_is_ignored_rather_than_rejected() {
        for repeat in [
            r#"{"axis":"x","step":"1mil","n":1}"#,
            r#"{"axis":"x","step":1,"n":"two"}"#,
            r#"{"axis":"diagonal","step":1,"n":1}"#,
            r#"{"step":1,"n":1}"#,
            r#"{"axis":"x","step":1,"n":1,"mirror":"yes"}"#,
            r#""every 1 mil""#,
            "[1]",
            "7",
        ] {
            let json = format!(
                r#"{{"name":"L","unit":"mil","spec":[
                    {{"type":"dot","x":0,"y":1}},
                    {{"type":"circle","x":0,"y":0,"r":2,"repeat":{repeat}}}
                ]}}"#
            );
            let (desc, report) = import_ventum_reticle_with_report(&json).unwrap_or_else(|e| {
                panic!("repeat {repeat} must not fail the import, got {e:?}")
            });

            // The document imported: the sibling dot is still a hold point...
            assert_eq!(desc.marks.len(), 1, "repeat {repeat} — the dot must survive");
            assert_eq!(desc.marks[0].down_mil, 1.0);
            // ...the ring is counted once, as an unrepeated circle...
            assert_eq!(
                report.dropped_element_types,
                vec![("circle".to_string(), 1)],
                "repeat {repeat} — expanded as though it had no repeat"
            );
            // ...and the tally says it may be low, rather than passing itself off as exact.
            assert_eq!(
                report.circle_repeats_unreadable, 1,
                "repeat {repeat} — an unreadable repeat must be declared, not swallowed"
            );
        }
    }

    /// A `repeat` a mark carries stays STRICT, because it stamps hold points: degrading it to
    /// one copy would drop marks in silence. Leniency is the circle's, not [`Repeat`]'s.
    #[test]
    fn an_unusable_repeat_on_a_hold_bearing_mark_still_fails_loudly() {
        for kind in [
            r#"{"type":"dot","x":0,"y":1"#,
            r#"{"type":"tick","x":0,"y":1"#,
            r#"{"type":"text","x":0,"y":1,"text":"5""#,
        ] {
            let json = format!(
                r#"{{"name":"S","unit":"mil","spec":[
                    {kind},"repeat":{{"axis":"x","step":"1mil","n":3}}}}
                ]}}"#
            );
            let error = import_ventum_reticle(&json)
                .expect_err("a mark's malformed repeat must not be swallowed");
            assert!(
                matches!(error, ReticleError::InvalidSpec(_)),
                "{kind} — expected InvalidSpec, got {error:?}"
            );
        }
    }

    /// A readable `repeat` on a circle is untouched by that leniency: it still stamps its
    /// copies, and the report still counts every one of them.
    #[test]
    fn a_readable_repeat_on_a_circle_still_stamps_its_copies() {
        let (_, report) = import_ventum_reticle_with_report(
            r#"{"name":"L","unit":"mil","spec":[
                {"type":"circle","x":1,"y":0,"r":0.5,"repeat":{"axis":"x","step":1,"n":4}}
            ]}"#,
        )
        .unwrap();
        assert_eq!(report.dropped_element_types, vec![("circle".to_string(), 4)]);
        assert_eq!(report.circle_repeats_unreadable, 0);
    }

    /// The mirror dedupe drops a twin only when it can PROVE the twin is a duplicate, and
    /// geometry the document wrote but this module cannot read proves nothing. A centered arc
    /// whose sweep or center is unreadable therefore keeps its twin: both are counted, both
    /// are unresolved, and the caller is told two shapes were lost. Guessing "duplicate" would
    /// undercount in exactly the way this report exists to prevent.
    #[test]
    fn an_unreadable_centered_arc_keeps_its_mirrored_twin() {
        // start + end = 540 = 180 (mod 360) reads as symmetric about the vertical axis and
        // would be deduped to one — but only when the angles can be read at all.
        for spec in [
            r#""start":"200deg","end":340"#,
            r#""start":200,"end":"340deg""#,
        ] {
            let json = format!(
                r#"{{"name":"U","unit":"mil","spec":[
                    {{"type":"circle","x":0,"y":0,"r":2,{spec},
                     "repeat":{{"axis":"x","step":1,"n":1,"mirror":true}}}}
                ]}}"#
            );
            let (_, report) = import_ventum_reticle_with_report(&json).unwrap();
            assert_eq!(
                report.dropped_element_types,
                vec![("arc".to_string(), 2)],
                "{spec} — an unprovable dedupe keeps the twin"
            );
            assert_eq!(report.arcs_unresolved, 2, "{spec}");
            assert!(report.arcs.is_empty(), "{spec}");
        }

        // Same for a center that cannot be read: the stepped coordinate this dedupe tests is
        // itself the fallback, so it proves nothing either.
        let (_, report) = import_ventum_reticle_with_report(
            r#"{"name":"U","unit":"mil","spec":[
                {"type":"circle","x":"left","y":0,"r":2,"start":200,"end":340,
                 "repeat":{"axis":"x","step":1,"n":1,"mirror":true}}
            ]}"#,
        )
        .unwrap();
        assert_eq!(report.dropped_element_types, vec![("arc".to_string(), 2)]);
        assert_eq!(report.arcs_unresolved, 2);
    }

    /// KNOWN DIVERGENCE from 0.32.0, recorded so it is a fact rather than a surprise.
    ///
    /// Resolving arcs put `circle` instances through repeat expansion, and expansion enforces
    /// [`MAX_RETICLE_MARKS`] over everything it materializes. A `circle` therefore consumes
    /// mark budget it never consumed before — under 0.32.0 it was skipped before expansion
    /// began — so a document already at the cap, or one whose circle carries a huge
    /// `repeat.n`, is refused where 0.32.0 imported it. Leniency cannot reach this: the
    /// document is well-formed, it is the cap that refuses it. Both CHANGELOG and CLI_USAGE
    /// state the exception; this pins it.
    #[test]
    fn a_circle_counts_against_the_mark_cap_unlike_0_32_0() {
        // One dot short of the cap, plus one ring: fine, because the ring is the last slot.
        let at_cap = format!(
            r#"{{"name":"C","unit":"mil","spec":[
                {{"type":"dot","x":1,"y":1,"repeat":{{"axis":"x","step":0.001,"n":{}}}}},
                {{"type":"circle","x":0,"y":0,"r":1}}
            ]}}"#,
            MAX_RETICLE_MARKS - 1
        );
        let (desc, report) = import_ventum_reticle_with_report(&at_cap).unwrap();
        assert_eq!(desc.marks.len(), MAX_RETICLE_MARKS - 1);
        assert_eq!(report.dropped_element_types, vec![("circle".to_string(), 1)]);

        // A full cap of marks plus one ring: refused, where 0.32.0 imported the marks and
        // dropped the ring unexamined.
        let over_cap = format!(
            r#"{{"name":"C","unit":"mil","spec":[
                {{"type":"dot","x":1,"y":1,"repeat":{{"axis":"x","step":0.001,"n":{MAX_RETICLE_MARKS}}}}},
                {{"type":"circle","x":0,"y":0,"r":1}}
            ]}}"#
        );
        assert!(
            matches!(
                import_ventum_reticle(&over_cap),
                Err(ReticleError::TooManyMarks { .. })
            ),
            "a circle past the cap is refused — the documented 0.32.0 divergence"
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

    // -----------------------------------------------------------------------------------
    // The strictness claims, pinned.
    //
    // The module's serde section, the CHANGELOG and CLI_USAGE all say the same thing about
    // this importer in prose — which fields can refuse a document, which absorb anything,
    // and what an arc needs before its points are reported. Prose drifts; three separate
    // rounds of review caught a sentence here that had gone one member out of date. These
    // sweeps are what those sentences point at: each states its claim as code, then drives
    // every case it covers through the real importer, so a change to the model that
    // contradicts the prose turns a test red instead of leaving a lie in a document.
    // -----------------------------------------------------------------------------------

    /// Values the Ventum schema declares no field as. `null` is deliberately absent: it is
    /// how JSON spells "no value", every optional field in this model reads it that way, and
    /// it has its own sweep in [`null_is_read_as_absence_not_as_a_bad_value`].
    const WRONG_TYPED: [&str; 4] = [r#""a string""#, r#"{"nested":1}"#, "[1]", "true"];

    /// [`WRONG_TYPED`] minus anything that is in fact the right type for `key`. `text` is the
    /// schema's only string-valued key, so a string is a perfectly good value there and is
    /// not evidence of anything; every other key is a number or the `repeat` object. (An
    /// object IS the right shape for a `repeat`, but `{"nested":1}` is missing every field
    /// one needs, so it is still a value a strict `repeat` refuses.)
    fn wrong_typed_for(key: &str) -> Vec<&'static str> {
        WRONG_TYPED
            .into_iter()
            .filter(|value| !(key == "text" && value.starts_with('"')))
            .collect()
    }

    /// Every key this importer reads on any element — the union across the variants — plus
    /// two cosmetic ones it never reads, so the sweeps cover "a key this element does not
    /// have" as well as "a key it does". [`strictness_is_exactly_what_a_hold_is_built_from`]
    /// asserts that nothing [`is_strict`] names is missing from this list.
    const EVERY_ELEMENT_KEY: [&str; 11] = [
        "x",
        "y",
        "cx",
        "cy",
        "r",
        "start",
        "end",
        "text",
        "repeat",
        "color",
        "max_extent",
    ];

    /// The element types this model can produce, taken FROM the model instead of from
    /// memory: the match below is exhaustive, so adding a variant to [`Element`] stops this
    /// file compiling until the new type is named — and once named it is swept with every
    /// key above.
    ///
    /// What the guard does not do, stated so nobody trusts it further than it goes: the
    /// sample list it maps over is hand-written, so a new variant needs adding in two places
    /// rather than one. The compile error lands on the match, three lines from the list.
    fn every_element_type() -> Vec<&'static str> {
        let samples = [
            Element::Dot {
                x: 0.0,
                y: 0.0,
                repeat: None,
            },
            Element::Tick {
                x: 0.0,
                y: 0.0,
                repeat: None,
            },
            Element::Text {
                x: 0.0,
                y: 0.0,
                text: String::new(),
                repeat: None,
            },
            Element::Line {},
            Element::Circle(CircleFields::default()),
            Element::Rect {},
            Element::Grid {},
            Element::Unknown,
        ];
        samples
            .iter()
            .map(|element| match element {
                Element::Dot { .. } => "dot",
                Element::Tick { .. } => "tick",
                Element::Text { .. } => "text",
                Element::Line {} => "line",
                Element::Circle(_) => "circle",
                Element::Rect {} => "rect",
                Element::Grid {} => "grid",
                // The schema defines no such type; `Unknown` is this model's answer for
                // anything it has never heard of, so the sweeps drive it with a name no
                // schema will define.
                Element::Unknown => "a_type_this_module_has_never_heard_of",
            })
            .collect()
    }

    /// The keys that make an element of `element_type` well-formed, so a sweep changing one
    /// key at a time is testing that key rather than a missing sibling.
    fn valid_base(element_type: &str) -> &'static [(&'static str, &'static str)] {
        match element_type {
            "dot" | "tick" => &[("x", "0"), ("y", "0")],
            "text" => &[("x", "0"), ("y", "0"), ("text", "\"t\"")],
            "circle" => &[("r", "1")],
            _ => &[],
        }
    }

    /// THE CLAIM, written where a test can disagree with it: the (element type, key) pairs
    /// that refuse a whole document over a value of the wrong type. Every other pair absorbs
    /// one.
    ///
    /// This is what the module's serde section says in words — strict is what a hold is
    /// built from: a `dot`/`tick` position, a `text`'s own position and string, and the
    /// `repeat` that stamps copies of any of those.
    fn is_strict(element_type: &str, key: &str) -> bool {
        matches!(
            (element_type, key),
            ("dot" | "tick", "x" | "y" | "repeat") | ("text", "x" | "y" | "text" | "repeat")
        )
    }

    /// A one-element document of `element_type`, well-formed except that `key` holds `value`
    /// (raw JSON). A base key of the same name is replaced rather than repeated, so the
    /// document tests the value's type and not a duplicate key.
    fn one_element_document(element_type: &str, key: &str, value: &str) -> String {
        let mut keys: Vec<String> = vec![format!(r#""type":"{element_type}""#)];
        for (base_key, base_value) in valid_base(element_type) {
            if *base_key != key {
                keys.push(format!(r#""{base_key}":{base_value}"#));
            }
        }
        keys.push(format!(r#""{key}":{value}"#));
        format!(
            r#"{{"name":"S","unit":"mil","spec":[{{{}}}]}}"#,
            keys.join(",")
        )
    }

    /// Every element type this model knows, crossed with every key the schema puts on one,
    /// crossed with values the schema never allows: the document must be refused for exactly
    /// the pairs [`is_strict`] names and imported for all the rest.
    ///
    /// Both directions matter and both are asserted by the same line. Loosening a mark's
    /// field turns this red (a refusal that stopped happening); tightening any field of a
    /// `circle`, or adding a typed one to it, turns it red too (a refusal that started).
    #[test]
    fn strictness_is_exactly_what_a_hold_is_built_from() {
        let mut refusals = 0usize;
        let mut acceptances = 0usize;
        for element_type in every_element_type() {
            for key in EVERY_ELEMENT_KEY {
                for value in wrong_typed_for(key) {
                    let json = one_element_document(element_type, key, value);
                    let refused = import_ventum_reticle(&json).is_err();
                    let expected = is_strict(element_type, key);
                    assert_eq!(
                        refused, expected,
                        "{element_type}.{key} = {value}: the importer refused={refused}, the \
                         documented boundary says strict={expected}\n  {json}"
                    );
                    if expected {
                        refusals += 1;
                    } else {
                        acceptances += 1;
                    }
                }
            }
        }
        // A sweep that swept nothing, or only one side of the boundary, would pass silently.
        // `dot` and `tick` each contribute x, y and repeat against four wrong values; `text`
        // those three plus its own string, which has only three wrong values because a
        // string is the right type for it.
        assert_eq!(refusals, 12 + 12 + 15, "the strict side must be covered");
        assert!(acceptances > 100, "the lenient side must be covered too");

        // And nothing the claim names may sit outside the keys the sweep drives, or it would
        // be asserted about nothing at all.
        for element_type in every_element_type() {
            for key in ["x", "y", "cx", "cy", "r", "start", "end", "text", "repeat", "label"] {
                assert!(
                    !is_strict(element_type, key) || EVERY_ELEMENT_KEY.contains(&key),
                    "{element_type}.{key} is claimed strict but is not in the swept key list"
                );
            }
        }
    }

    /// The other direction, at the one element the leniency is about, with the two shapes a
    /// type sweep cannot reach: a key written twice, and the schema's two spellings of one
    /// center used together.
    ///
    /// Both of those refuse a derived struct with `duplicate field \`x\``, and 0.32.0
    /// accepted them because its `circle` read no keys at all. A ring must not be able to
    /// fail a document over how it spells a center, so [`CircleFields`] reads the map itself.
    #[test]
    fn no_circle_key_refuses_a_document_however_it_is_spelled_or_repeated() {
        const CIRCLE_KEYS: [&str; 8] = ["x", "y", "cx", "cy", "r", "start", "end", "repeat"];
        let values: Vec<&str> = WRONG_TYPED
            .iter()
            .copied()
            .chain(["null", "0", "-1", "1e308"])
            .collect();

        for key in CIRCLE_KEYS {
            for value in &values {
                // Once, twice, and — for a center — beside its other spelling.
                let mut circles = vec![
                    format!(r#"{{"type":"circle","{key}":{value}}}"#),
                    format!(r#"{{"type":"circle","{key}":{value},"{key}":{value}}}"#),
                ];
                if let Some(other) = match key {
                    "x" => Some("cx"),
                    "cx" => Some("x"),
                    "y" => Some("cy"),
                    "cy" => Some("y"),
                    _ => None,
                } {
                    circles.push(format!(r#"{{"type":"circle","{other}":1,"{key}":{value}}}"#));
                    circles.push(format!(r#"{{"type":"circle","{key}":{value},"{other}":1}}"#));
                }

                for circle in circles {
                    let json = format!(
                        r#"{{"name":"L","unit":"mil","spec":[
                            {{"type":"dot","x":0,"y":1}},{circle}
                        ]}}"#
                    );
                    let (desc, report) = import_ventum_reticle_with_report(&json)
                        .unwrap_or_else(|e| panic!("{circle} must not fail the import: {e:?}"));
                    assert_eq!(desc.marks.len(), 1, "{circle} — the dot must survive");
                    assert_eq!(
                        report.dropped_elements, 1,
                        "{circle} — the ring must still be counted, not silently absorbed"
                    );
                }
            }
        }
    }

    /// And the rule that resolves a repeated or double-spelled key: the LAST one in document
    /// order wins, which is the answer `serde_json::Value` itself gives a repeated key.
    #[test]
    fn a_repeated_circle_key_takes_its_last_value() {
        // `x` and `cx` are one field, so the later spelling is the one that counts...
        for (spec, expected_right) in [
            (r#""x":5,"cx":1"#, 1.0),
            (r#""cx":1,"x":5"#, 5.0),
            (r#""x":5,"x":1"#, 1.0),
        ] {
            let json = format!(
                r#"{{"name":"O","unit":"mil","spec":[
                    {{"type":"circle",{spec},"y":0,"r":2,"start":200,"end":340}}
                ]}}"#
            );
            let (_, report) = import_ventum_reticle_with_report(&json).unwrap();
            assert_eq!(report.arcs.len(), 1, "{spec}");
            assert!(
                (report.arcs[0].center.right_mil - expected_right).abs() < 1e-12,
                "{spec}: center {} should be {expected_right}",
                report.arcs[0].center.right_mil
            );
        }

        // ...and the same for a repeat, whose last value here is `null` — no repeat at all,
        // so one ring is counted rather than three.
        let (_, report) = import_ventum_reticle_with_report(
            r#"{"name":"O","unit":"mil","spec":[
                {"type":"circle","r":2,"repeat":{"axis":"x","step":1,"n":3},"repeat":null}
            ]}"#,
        )
        .unwrap();
        assert_eq!(report.dropped_element_types, vec![("circle".to_string(), 1)]);
        assert_eq!(report.circle_repeats_unreadable, 0, "`null` is absence, not garbage");
    }

    /// What an arc's points are built from, and therefore what leaves one unresolved.
    ///
    /// [`CircleShape::resolve`] needs three things and computes nothing without all of them:
    /// a center it could read, a radius that is a positive length, and a sweep. This drives
    /// each of the three away in turn — including the case an earlier version of the
    /// `arcs_unresolved` doc comment left out, a radius that reads perfectly well as a
    /// number and is not a length (`0`, or negative).
    #[test]
    fn an_arc_resolves_only_from_a_center_a_positive_radius_and_a_sweep() {
        // The positive control: all three readable, so the arc is resolved rather than
        // counted as lost. Without this the sweep below could pass on an importer that
        // resolved nothing at all.
        let (_, resolved) = import_ventum_reticle_with_report(
            r#"{"name":"A","unit":"mil","spec":[
                {"type":"circle","x":0,"y":0,"r":2,"start":200,"end":340}
            ]}"#,
        )
        .unwrap();
        assert_eq!(resolved.arcs.len(), 1);
        assert_eq!(resolved.arcs_unresolved, 0);

        for (what, circle) in [
            ("no radius", r#""x":0,"y":0,"start":200,"end":340"#),
            ("unreadable radius", r#""x":0,"y":0,"r":"2mil","start":200,"end":340"#),
            ("zero radius", r#""x":0,"y":0,"r":0,"start":200,"end":340"#),
            ("negative radius", r#""x":0,"y":0,"r":-2,"start":200,"end":340"#),
            ("unreadable start", r#""x":0,"y":0,"r":2,"start":"200deg","end":340"#),
            ("unreadable end", r#""x":0,"y":0,"r":2,"start":200,"end":"340deg""#),
            ("unreadable center x", r#""x":"left","y":0,"r":2,"start":200,"end":340"#),
            ("unreadable center y", r#""x":0,"y":"up","r":2,"start":200,"end":340"#),
        ] {
            let json =
                format!(r#"{{"name":"A","unit":"mil","spec":[{{"type":"circle",{circle}}}]}}"#);
            let (_, report) = import_ventum_reticle_with_report(&json)
                .unwrap_or_else(|e| panic!("{what} must not fail the import: {e:?}"));
            assert!(report.arcs.is_empty(), "{what}: no points may be reported");
            assert_eq!(report.arcs_unresolved, 1, "{what}: must be declared unresolved");
            assert_eq!(
                report.dropped_element_types,
                vec![(ARC_TAG.to_string(), 1)],
                "{what}: keeps its arc tag — the loss is reported as an arc's"
            );
        }
    }

    /// `null` is how JSON writes "no value", and this model reads it as absence everywhere a
    /// value is optional — never as a value it could not understand.
    ///
    /// The consequence worth pinning is the one that looks like an inconsistency until the
    /// rule is stated: a `circle` carrying `"start":null` alongside an unreadable `end` is a
    /// RING, not an arc with a lost sweep. A sweep needs two angles, `null` supplies neither,
    /// and one angle alone describes no sweep however it is spelled.
    #[test]
    fn null_is_read_as_absence_not_as_a_bad_value() {
        for (what, circle) in [
            ("start null, end unreadable", r#""r":2,"start":null,"end":"340deg""#),
            ("start unreadable, end null", r#""r":2,"start":"200deg","end":null"#),
            ("both null", r#""r":2,"start":null,"end":null"#),
            ("every key null", r#""x":null,"y":null,"r":null,"start":null,"end":null"#),
        ] {
            let json =
                format!(r#"{{"name":"N","unit":"mil","spec":[{{"type":"circle",{circle}}}]}}"#);
            let (_, report) = import_ventum_reticle_with_report(&json).unwrap();
            assert_eq!(
                report.dropped_element_types,
                vec![("circle".to_string(), 1)],
                "{what}: a declared-but-null angle is no angle, so this is a ring"
            );
            assert_eq!(report.arcs_unresolved, 0, "{what}: a ring loses no sweep");
        }

        // The one other field in this model that takes a `null`: a mark's optional `repeat`,
        // where serde's `Option` reads it as absence too. One dot is emitted.
        let desc = import_ventum_reticle(
            r#"{"name":"N","unit":"mil","spec":[{"type":"dot","x":0,"y":1,"repeat":null}]}"#,
        )
        .unwrap();
        assert_eq!(desc.marks.len(), 1);

        // And the half that is NOT true, pinned so the sentence above cannot quietly grow
        // into "null is absence everywhere". A reticle-level `#[serde(default)]` field fills
        // in a key the document omitted and still refuses a present `null` — 0.32.0's
        // behaviour, left alone — as does a defaulted field inside a mark's `repeat`.
        for metadata in [
            r#""name":null"#,
            r#""plane":null"#,
            r#""unit":null"#,
            r#""ref_magnification":null"#,
            r#""spec":null"#,
        ] {
            // The null key on its own: every other reticle-level field has a default, so
            // nothing but the null under test can be what refuses this document.
            let json = format!(r#"{{{metadata}}}"#);
            assert!(
                import_ventum_reticle(&json).is_err(),
                "{metadata}: a present null is not an absent key at reticle level\n  {json}"
            );
        }
        assert!(
            import_ventum_reticle(
                r#"{"unit":"mil","spec":[{"type":"dot","x":0,"y":1,
                   "repeat":{"axis":"x","step":1,"n":2,"mirror":null}}]}"#
            )
            .is_err(),
            "a null inside a mark's repeat is not an absent key either"
        );
    }

    /// What a `repeat` does to the tally, which is not one answer for every type.
    ///
    /// A `text` is expanded before it is bound, so each copy that lands on no mark is its own
    /// lost label and its own count. A `line` is counted and skipped before expansion begins,
    /// so its `repeat` is never read and the element counts once however many copies the
    /// document draws. The report's "How the counts are taken" section says exactly this; the
    /// two numbers below are the difference it is talking about.
    #[test]
    fn a_repeat_counts_per_instance_for_text_and_once_for_a_dropped_shape() {
        let (_, text_report) = import_ventum_reticle_with_report(
            r#"{"name":"R","unit":"mil","spec":[
                {"type":"text","x":20,"y":20,"text":"t","repeat":{"axis":"x","step":1,"n":3}}
            ]}"#,
        )
        .unwrap();
        assert_eq!(
            text_report.dropped_element_types,
            vec![("text".to_string(), 3)],
            "three copies that bound to nothing are three lost labels"
        );

        let (_, line_report) = import_ventum_reticle_with_report(
            r#"{"name":"R","unit":"mil","spec":[
                {"type":"line","x1":0,"y1":0,"x2":1,"y2":1,
                 "repeat":{"axis":"x","step":1,"n":3}}
            ]}"#,
        )
        .unwrap();
        assert_eq!(
            line_report.dropped_element_types,
            vec![("line".to_string(), 1)],
            "a dropped shape's repeat is never expanded, so the element counts once"
        );
    }

    /// A refusal that is not about the `circle` at all, pinned as such.
    ///
    /// Three of these exist and none of them is this module's to give away. A numeric literal
    /// outside `f64`'s range and JSON nested past serde_json's recursion limit are refused by
    /// the PARSER before any element model is consulted; an element that writes its own
    /// `type` twice is refused by serde's tag reader before a variant is even chosen. Nothing
    /// in this module's leniency can see any of them, and 0.32.0 refused all three.
    ///
    /// The evidence is that they fall identically on element types that read no fields
    /// whatsoever: if one of these documents ever imports as a `line` but not as a `circle`,
    /// the difference IS about the circle and this turns red.
    #[test]
    fn a_refusal_that_is_not_about_a_key_falls_on_every_element_type_alike() {
        let deep = format!("{}1{}", "[".repeat(200), "]".repeat(200));
        for value in [String::from("1e400"), deep] {
            for element_type in every_element_type() {
                let json = one_element_document(element_type, "r", &value);
                assert!(
                    import_ventum_reticle(&json).is_err(),
                    "{element_type} with r={value} must be refused by the parser\n  {json}"
                );
            }
        }

        // The element's own tag, written twice. Not a key this module reads, and not a key
        // any leniency here could reach: the tag is consumed before the variant exists.
        let mut refusals = 0usize;
        for element_type in every_element_type() {
            let json = format!(
                r#"{{"name":"S","unit":"mil","spec":[
                    {{"type":"{element_type}","type":"{element_type}","x":0,"y":1}}
                ]}}"#
            );
            assert!(
                import_ventum_reticle(&json).is_err(),
                "{element_type} written with two `type` keys must be refused\n  {json}"
            );
            refusals += 1;
        }
        assert_eq!(refusals, every_element_type().len());
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
