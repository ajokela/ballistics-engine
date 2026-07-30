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
//! # Coordinate and unit conventions
//!
//! Ventum uses `+x = right`, `+y = down`, origin at the reticle center — identical to the
//! engine's [`ReticleMark::right_mil`] / [`ReticleMark::down_mil`], so no sign flips are
//! needed. All coordinates are expressed in the reticle's own `unit`; MOA input is
//! converted to milliradians with `1 MOA = 0.2908882 mil` ([`MOA_TO_MIL`]). Auto-numbered
//! ladder labels are NOT unit-converted — a label is the reticle-unit value it names.
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
/// # Errors
///
/// Returns [`ReticleError::InvalidSpec`] if `json` cannot be parsed as a Ventum reticle,
/// and [`ReticleError::TooManyMarks`] if repeat expansion would exceed the mark cap.
pub fn import_ventum_reticle(json: &str) -> Result<ReticleDescription, ReticleError> {
    let reticle: VentumReticle =
        serde_json::from_str(json).map_err(|e| ReticleError::InvalidSpec(e.to_string()))?;

    let scale = match reticle.unit {
        Unit::Mil => 1.0,
        Unit::Moa => MOA_TO_MIL,
    };

    // Expand every element's `repeat` into concrete, unit-scaled instances first, with the
    // mark cap enforced during expansion.
    let instances = expand_elements(&reticle.spec.0, scale)?;

    // Classify: dot/tick become marks; text is collected for nearest-mark binding.
    let mut marks: Vec<ReticleMark> = Vec::new();
    let mut texts: Vec<TextLabel> = Vec::new();
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
        }
    }

    bind_text_labels(&mut marks, &texts);

    // `ref_magnification` is meaningful only for SFP; an FFP reticle's subtensions do not
    // depend on magnification, so we normalize its reference to 1.0 (as the engine does).
    let reference_magnification = match reticle.plane {
        FocalPlane::Second => reticle.ref_magnification,
        FocalPlane::First => 1.0,
    };

    Ok(ReticleDescription {
        name: reticle.name,
        focal_plane: reticle.plane,
        reference_magnification,
        marks,
    })
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
fn bind_text_labels(marks: &mut [ReticleMark], texts: &[TextLabel]) {
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
        if let Some(index) = nearest {
            if nearest_distance <= DEFAULT_TEXT_BIND_MIL {
                marks[index].label = Some(text.label.clone());
            }
        }
    }
}

/// The role an expanded point instance plays once classified.
#[derive(Debug, Clone, PartialEq)]
enum Role {
    Dot,
    Tick,
    /// A text label carrying its (possibly auto-numbered) string.
    Text(String),
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
}

impl PointRole<'_> {
    /// The [`Role`] for one instance. `ladder_label` is `Some` only for an auto-numbered
    /// text ladder copy; otherwise a text keeps its element's own string.
    fn role_with_label(&self, ladder_label: Option<String>) -> Role {
        match self {
            PointRole::Dot => Role::Dot,
            PointRole::Tick => Role::Tick,
            PointRole::Text(base) => Role::Text(ladder_label.unwrap_or_else(|| (*base).to_string())),
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
/// and enforcing the mark cap along the way. Decoration and unknown element types are
/// dropped (not emitted, not counted).
fn expand_elements(elements: &[Element], scale: f64) -> Result<Vec<ExpandedInstance>, ReticleError> {
    let mut out: Vec<ExpandedInstance> = Vec::new();
    for element in elements {
        let (x, y, repeat, base) = match element {
            Element::Dot { x, y, repeat } => (*x, *y, repeat.as_ref(), PointRole::Dot),
            Element::Tick { x, y, repeat } => (*x, *y, repeat.as_ref(), PointRole::Tick),
            Element::Text {
                x, y, text, repeat,
            } => (*x, *y, repeat.as_ref(), PointRole::Text(text)),
            // Decoration and unknown future element types carry no hold and no label.
            Element::Line {}
            | Element::Circle {}
            | Element::Rect {}
            | Element::Grid {}
            | Element::Unknown => continue,
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

        if repeat.mirror && stepped != 0.0 {
            // Emit the axis-negated twin as well, but never a duplicate at the center. Only
            // the stepped axis is negated. A mirrored labeled-text copy keeps the SAME label
            // as its positive twin (magnitude convention: the "5" on the right and the "5"
            // on the left both read 5). The schema does not specify this; it is a deliberate
            // assumption.
            let (mirror_down, mirror_right) = match repeat.axis {
                Axis::X => (down, -right),
                Axis::Y => (-down, right),
            };
            push_capped(down, right, role.clone(), out)?;
            push_capped(mirror_down, mirror_right, role, out)?;
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
// notes, manufacturer, ...) are ignored rather than rejected.
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
    /// Decoration — dropped.
    Circle {},
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
    /// When true, also emit the axis-negated copy (a `±` symmetric ladder), except at the
    /// center.
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
        let instances = expand_elements(&elements, 1.0).unwrap();
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
