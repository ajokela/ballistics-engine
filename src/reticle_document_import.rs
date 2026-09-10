//! Cleanroom importer for the `.reticle` XML document format.
//!
//! A `.reticle` file is a third-party **drawing** format for scope reticles: a single
//! `<reticle>` root carrying an `<elements>` list of shapes (lines, circles, rectangles,
//! filled paths) plus a `<bdc>` section. This module is the one-way transform from that
//! authoring format into the engine's own [`ReticleDescription`], so a reticle drawn in
//! that format can be hold-solved by the existing
//! [`hold_point_in_reticle`](crate::reticle::hold_point_in_reticle) without any of the
//! CLI/FFI/WASM surfaces having to learn a second schema. It is a converter, not new
//! physics — every coordinate is carried straight through, only unit-normalized to
//! milliradians and sign-corrected for the engine's axis convention.
//!
//! It is a sibling of [`crate::reticle_import`] (Bero's "Ventum" JSON), and meets the same
//! hold-vs-decoration problem from the opposite direction: Ventum is an aiming-mark format
//! that happens to contain some decoration, while `.reticle` is a picture that happens to
//! contain some aiming marks.
//!
//! Everything here was derived from a sample `.reticle` document and the SVG it renders
//! to — the observable file structure, and the observable output that pins down what its
//! coordinates mean. **No third-party implementation of this format was consulted, and
//! none is vendored.** A file format is not itself copyrightable; this crate stays MIT OR
//! Apache-2.0. Every fixture in this module's tests and in
//! `tests/fixtures/reticle_document/` is invented content.
//!
//! Like [`crate::reticle_import`], this module is fs-free (it transforms an in-memory
//! `&str`), carries no feature gate, and compiles for `wasm32-unknown-unknown`.
//!
//! # What carries a hold, and what does not
//!
//! This is the whole design decision, and the format answers it itself rather than leaving
//! it to be inferred:
//!
//! * **`<bdc>` entries are the holds.** A `<bdc position-x="…" position-y="…">` element is
//!   a declared aiming point. Each one becomes a [`ReticleMark`].
//! * **Everything under `<elements>` is the picture, and is dropped.** `reticle-line`,
//!   `reticle-circle`, `reticle-rectangle`, `reticle-path` (and its
//!   `reticle-path-move-to` / `reticle-path-line-to` vertices), and any unknown future tag
//!   carry no hold. Their geometry attributes are never read as angles at all — a shape
//!   this module drops cannot fail its import.
//!
//! The evidence that this split is the format's own and not our invention: in a real
//! document the `<bdc>` positions are *duplicates* of shapes already drawn under
//! `<elements>`. A BDC reticle draws its three holdover marks as three short
//! `reticle-line` ticks, and then lists those same three coordinates again under `<bdc>`.
//! The drawing says how it looks; the `<bdc>` section says what it means. Importing the
//! drawing as well would double every holdover.
//!
//! ## The cost of that rule, stated plainly
//!
//! A `.reticle` document's windage marks are drawn as `reticle-line` ticks and are *not*
//! listed under `<bdc>`. They are therefore dropped, and an imported reticle usually has
//! marks on the elevation axis only. That has a visible consequence downstream: with every
//! mark at `right_mil == 0`, the windage axis of the mark bounding box is degenerate, so
//! [`ReticleHold::off_reticle`](crate::reticle::ReticleHold::off_reticle) reads `true` for
//! *any* nonzero wind hold — see that field's documentation. For a reticle that really does
//! have windage hashes etched on it, that answer comes from this import, not from the optic.
//!
//! Recovering those marks would mean classifying `reticle-line` elements as hash marks by
//! their geometry (short, axis-perpendicular, off-center) — inventing holds out of
//! decoration, and doing it with a heuristic that cannot distinguish a windage hash from a
//! crosshair arm. This module will not guess. Instead the drop is made countable rather
//! than silent: [`import_reticle_document_with_report`] returns a
//! [`ReticleImportReport`] tallying every drawing tag it walked past, so a caller can tell
//! the shooter how many lines and paths this document drew that are not holds, instead of
//! leaving them to infer it from a suspiciously empty reticle.
//!
//! Two related non-inventions, for the reader who wonders:
//!
//! * **No synthetic center mark.** The zero *is* the coordinate origin every imported mark
//!   is measured from, and a hold is reported in those same coordinates, so it survives as
//!   the frame rather than as a mark. If a document wants the center to be a selectable
//!   aiming point it can say so with a `<bdc position-x="0…" position-y="0…">` entry.
//! * **No labels.** A `<bdc>` entry carries no range and no text — only a position and the
//!   `text-offset` / `text-height` hints for where a renderer should draw a label. The range
//!   that belongs to a holdover is a property of the firing solution, not of the document,
//!   so every imported mark has `label: None`.
//!
//! # Coordinate and unit conventions
//!
//! **The document's `+y` is UP; the engine's is DOWN.** This is the one sign trap of the
//! format, and it is the opposite of [`crate::reticle_import`], where no flip was needed.
//! Confirmed from rendered output: an element at `position-y="8moa"` renders *above* the
//! zero, and a holdover BDC entry sits at a *negative* `position-y`. So
//!
//! ```text
//!   right_mil =  position-x           down_mil = -position-y
//! ```
//!
//! and a holdover's negative `position-y` becomes the positive `down_mil` the engine wants.
//!
//! Every value names its own unit inline as a suffix — `"96moa"`, `"-0.55moa"`,
//! `"1.5mil"`. There is no document-level unit declaration, so each attribute is converted
//! on its own: `moa` scales by [`MOA_TO_MIL`] (shared with [`crate::reticle_import`], so
//! the two importers cannot drift apart), `mil` / `mils` / `mrad` pass through. A value
//! with no suffix, or with a suffix this module does not recognize, is a hard error rather
//! than a guessed scale — silently mis-scaling a reticle by 3.44x is worse than refusing it.
//!
//! ## `size-x` / `size-y` / `zero-x` / `zero-y` are canvas placement, not a mark offset
//!
//! The root declares a canvas extent (`size-x`, `size-y`) and the position of the zero
//! within it (`zero-x`, `zero-y`, measured from the canvas top-left). They exist so a
//! renderer can map a reticle coordinate onto its output surface:
//!
//! ```text
//!   canvas_x = (zero-x + position-x) / size-x * width
//!   canvas_y = (zero-y - position-y) / size-y * height     (note the minus: +y is up)
//! ```
//!
//! The origin of the element coordinate system *is* the zero — that is what the attributes
//! are named after — so element coordinates are already zero-relative and `zero-x`/`zero-y`
//! do **not** enter this transform. Changing them re-frames the picture on its canvas; it
//! does not move a hold relative to the point of aim. They are reported (as
//! [`ReticleImportReport::zero_in_canvas_mil`]) rather than applied, and there is a test
//! that two documents differing only in their zero attributes import to identical marks.
//!
//! # What the format does not carry
//!
//! `.reticle` has no focal-plane and no magnification attribute. A document's own name
//! string may hint at one in prose (`"… (4-12x40, at 12x)"`), but that is text for a human,
//! not a field, and this module does not mine it.
//!
//! An imported description is therefore always [`FocalPlane::First`] with
//! `reference_magnification = 1.0` — the identity, the only choice that scales nothing.
//! **If the reticle is actually SFP, the caller must set
//! [`ReticleDescription::focal_plane`] and
//! [`ReticleDescription::reference_magnification`] from the optic before solving**, or
//! every hold will be computed as though the subtensions were magnification-independent.
//! Both fields are public precisely so that is a one-line fix at the call site.
//!
//! # Safety
//!
//! [`crate::reticle::MAX_RETICLE_MARKS`] is enforced *while* `<bdc>` entries are collected,
//! so a hostile document with a million entries returns [`ReticleError::TooManyMarks`]
//! promptly instead of allocating. Element nesting is bounded by
//! [`MAX_DOCUMENT_DEPTH`], and the scanner is iterative — no recursion, so a deeply nested
//! document cannot overflow the stack. This is the same defense as the 0.31.0 reticle
//! generator size guards and the Ventum importer's repeat-expansion cap.

use std::collections::BTreeMap;

use crate::reticle::{
    FocalPlane, MarkKind, ReticleDescription, ReticleError, ReticleMark, MAX_RETICLE_MARKS,
};
use crate::reticle_import::MOA_TO_MIL;

/// Deepest element nesting a `.reticle` document may use before it is rejected.
///
/// Real documents nest three deep (`reticle` > `elements` > `reticle-path` > `elements` >
/// `reticle-path-line-to` is five). The cap exists so a hostile document cannot grow the
/// open-element stack without bound; nothing legitimate approaches it.
pub const MAX_DOCUMENT_DEPTH: usize = 64;

/// The document root.
const TAG_RETICLE: &str = "reticle";
/// A structural container for child shapes. Appears both directly under the root and
/// inside a `reticle-path` (whose vertices are its `<elements>`), so it is treated as
/// structure wherever it occurs.
const TAG_ELEMENTS: &str = "elements";
/// Both the hold-bearing entry and the section that wraps a list of them; the two are told
/// apart by whether the element carries a position (see [`bdc_mark`]).
const TAG_BDC: &str = "bdc";

const ATTR_NAME: &str = "name";
const ATTR_SIZE_X: &str = "size-x";
const ATTR_SIZE_Y: &str = "size-y";
const ATTR_ZERO_X: &str = "zero-x";
const ATTR_ZERO_Y: &str = "zero-y";
const ATTR_POSITION_X: &str = "position-x";
const ATTR_POSITION_Y: &str = "position-y";

/// What an import left behind, so the hold/decoration rule is countable rather than
/// implicit. See the module documentation for why the drawing is dropped.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ReticleImportReport {
    /// How many drawing tags were walked past — the sum of [`Self::dropped_element_tags`].
    ///
    /// This counts *tags*, not shapes: a filled `reticle-path` contributes itself plus one
    /// per `reticle-path-move-to` / `reticle-path-line-to` vertex, because that is what the
    /// document literally contains.
    pub drawing_elements_dropped: usize,
    /// Per-tag counts of that drawing, sorted by tag name so the report is deterministic.
    pub dropped_element_tags: Vec<(String, usize)>,
    /// The document's declared canvas extent as `(size-x, size-y)` in milliradians, when
    /// the root declares both. Rendering metadata — see the module documentation.
    pub canvas_mil: Option<(f64, f64)>,
    /// Where the zero sits inside that canvas, as `(zero-x, zero-y)` in milliradians from
    /// the canvas top-left, when the root declares both. Rendering metadata: it does NOT
    /// offset the imported marks.
    pub zero_in_canvas_mil: Option<(f64, f64)>,
}

/// Import a `.reticle` XML document into a [`ReticleDescription`].
///
/// `xml` is one `.reticle` document (see the module documentation for the grammar). The
/// returned description is ready for
/// [`hold_point_in_reticle`](crate::reticle::hold_point_in_reticle): its marks are in
/// milliradians from the zero, positive `down_mil` below it.
///
/// The description is always [`FocalPlane::First`] with `reference_magnification = 1.0`,
/// because the format carries neither; a caller holding an SFP optic must overwrite both
/// fields before solving.
///
/// This is a pure transform: like [`crate::reticle_import::import_ventum_reticle`] it does
/// not itself reject an all-decoration document (the resulting description simply carries
/// no marks, which [`hold_point_in_reticle`](crate::reticle::hold_point_in_reticle) then
/// reports as [`ReticleError::NoMarks`]). It DOES enforce
/// [`crate::reticle::MAX_RETICLE_MARKS`] while collecting.
///
/// Use [`import_reticle_document_with_report`] when you want to know what was dropped.
///
/// # Errors
///
/// Returns [`ReticleError::InvalidSpec`] if `xml` is not well-formed, is not rooted at
/// `<reticle>`, or carries an attribute value this module cannot read as an angle (missing
/// unit suffix, unrecognized unit, non-numeric, non-finite, non-positive canvas), and
/// [`ReticleError::TooManyMarks`] if the document declares more `<bdc>` entries than the
/// mark cap.
pub fn import_reticle_document(xml: &str) -> Result<ReticleDescription, ReticleError> {
    import_reticle_document_with_report(xml).map(|(description, _)| description)
}

/// Import a `.reticle` document, also returning a [`ReticleImportReport`] of what the
/// hold/decoration rule left behind.
///
/// Identical to [`import_reticle_document`] in every other respect; that function is this
/// one with the report discarded.
///
/// # Errors
///
/// As [`import_reticle_document`].
pub fn import_reticle_document_with_report(
    xml: &str,
) -> Result<(ReticleDescription, ReticleImportReport), ReticleError> {
    let mut scanner = XmlScanner::new(xml);

    // Names of the currently open elements, innermost last. A Vec rather than recursion:
    // depth is data here, and a document must not be able to nest us into a stack overflow.
    let mut open: Vec<&str> = Vec::new();
    let mut root_seen = false;

    let mut name = String::new();
    let (mut size_x, mut size_y) = (None, None);
    let (mut zero_x, mut zero_y) = (None, None);
    let mut marks: Vec<ReticleMark> = Vec::new();
    // Keyed by a borrowed slice of the input, so a document with many distinct tag names
    // costs pointers rather than copies; owned only when the report is built.
    let mut dropped: BTreeMap<&str, usize> = BTreeMap::new();

    while let Some(event) = scanner.next_event()? {
        match event {
            XmlEvent::Start {
                tag,
                attributes,
                self_closing,
            } => {
                if open.is_empty() {
                    // Top level: this must be the one and only root, and it must be a
                    // `.reticle` document rather than some other XML that happens to parse.
                    if root_seen {
                        return Err(invalid(
                            "a .reticle document must have exactly one root element".to_string(),
                        ));
                    }
                    if tag != TAG_RETICLE {
                        return Err(invalid(format!(
                            "expected a <{TAG_RETICLE}> root element, found <{tag}>"
                        )));
                    }
                    root_seen = true;

                    if let Some(raw) = attribute(&attributes, ATTR_NAME) {
                        name = decode_entities(raw);
                    }
                    size_x = parse_canvas_extent(&attributes, ATTR_SIZE_X)?;
                    size_y = parse_canvas_extent(&attributes, ATTR_SIZE_Y)?;
                    zero_x = parse_optional_angle_mil(&attributes, ATTR_ZERO_X)?;
                    zero_y = parse_optional_angle_mil(&attributes, ATTR_ZERO_Y)?;
                } else {
                    match tag {
                        // Pure structure: neither a hold nor a drawn shape.
                        TAG_ELEMENTS => {}
                        // The hold-bearing section. A `<bdc>` that carries a position is an
                        // entry; one that carries none is the wrapper around entries, and
                        // contributes nothing itself.
                        TAG_BDC => {
                            if let Some(mark) = bdc_mark(&attributes)? {
                                push_capped(&mut marks, mark)?;
                            }
                        }
                        // Everything else is the picture. Its geometry attributes are never
                        // parsed — only tallied, so the drop is visible in the report.
                        drawing => *dropped.entry(drawing).or_insert(0) += 1,
                    }
                }

                if !self_closing {
                    if open.len() >= MAX_DOCUMENT_DEPTH {
                        return Err(invalid(format!(
                            "element nesting deeper than {MAX_DOCUMENT_DEPTH} at <{tag}>"
                        )));
                    }
                    open.push(tag);
                }
            }
            XmlEvent::End { tag } => match open.pop() {
                Some(expected) if expected == tag => {}
                Some(expected) => {
                    return Err(invalid(format!(
                        "</{tag}> closes <{expected}>, which is still open"
                    )))
                }
                None => return Err(invalid(format!("</{tag}> has no matching open element"))),
            },
        }
    }

    if let Some(unclosed) = open.last() {
        return Err(invalid(format!("<{unclosed}> is never closed")));
    }
    if !root_seen {
        return Err(invalid(format!(
            "the document has no <{TAG_RETICLE}> root element"
        )));
    }

    let report = ReticleImportReport {
        drawing_elements_dropped: dropped.values().sum(),
        dropped_element_tags: dropped
            .into_iter()
            .map(|(tag, count)| (tag.to_string(), count))
            .collect(),
        canvas_mil: size_x.zip(size_y),
        zero_in_canvas_mil: zero_x.zip(zero_y),
    };

    let description = ReticleDescription {
        name,
        // The format declares no focal plane, so the import applies the identity: FFP
        // subtensions do not scale with magnification. See the module documentation.
        focal_plane: FocalPlane::First,
        reference_magnification: 1.0,
        marks,
    };

    Ok((description, report))
}

/// Turn one `<bdc>` element into a mark, or `None` when it is the section wrapper rather
/// than an entry.
///
/// An entry is any `<bdc>` carrying at least one of `position-x` / `position-y`; the
/// wrapper carries neither. A missing axis defaults to `0` — an entry that names only its
/// `position-y` sits on the vertical stadia, which is where a holdover ladder lives.
fn bdc_mark(attributes: &[(&str, &str)]) -> Result<Option<ReticleMark>, ReticleError> {
    let raw_x = attribute(attributes, ATTR_POSITION_X);
    let raw_y = attribute(attributes, ATTR_POSITION_Y);
    if raw_x.is_none() && raw_y.is_none() {
        return Ok(None);
    }

    let right_mil = match raw_x {
        Some(raw) => parse_angle_mil(raw, ATTR_POSITION_X)?,
        None => 0.0,
    };
    let up_mil = match raw_y {
        Some(raw) => parse_angle_mil(raw, ATTR_POSITION_Y)?,
        None => 0.0,
    };

    // THE sign flip: the document measures +y UP from the zero, the engine measures
    // +down_mil DOWN from it. A holdover's negative position-y is a positive down_mil.
    // `text-offset` / `text-height` are label-rendering hints and are dropped; a `<bdc>`
    // entry carries no range text to become a mark label.
    Ok(Some(ReticleMark::new(
        -up_mil,
        right_mil,
        // The format gives a bdc point no shape (the shape lives in the drawing), so this
        // is the conventional drawing of a holdover point on a stadia line. `MarkKind` is
        // descriptive only and has no effect on the hold math.
        MarkKind::Hash,
    )))
}

/// Push one mark, rejecting the document the instant the running count would exceed
/// [`MAX_RETICLE_MARKS`]. Checking before every push bounds the allocation at the cap
/// regardless of how many `<bdc>` entries the document declares.
fn push_capped(marks: &mut Vec<ReticleMark>, mark: ReticleMark) -> Result<(), ReticleError> {
    if marks.len() >= MAX_RETICLE_MARKS {
        return Err(ReticleError::TooManyMarks {
            count: marks.len() + 1,
            max: MAX_RETICLE_MARKS,
        });
    }
    marks.push(mark);
    Ok(())
}

/// First value for `wanted`, or `None`. XML forbids duplicate attribute names; if a
/// document has them anyway the first wins rather than the parse failing.
fn attribute<'a>(attributes: &[(&'a str, &'a str)], wanted: &str) -> Option<&'a str> {
    attributes
        .iter()
        .find(|(key, _)| *key == wanted)
        .map(|(_, value)| *value)
}

/// Parse an optional angle-valued root attribute.
fn parse_optional_angle_mil(
    attributes: &[(&str, &str)],
    wanted: &'static str,
) -> Result<Option<f64>, ReticleError> {
    match attribute(attributes, wanted) {
        Some(raw) => parse_angle_mil(raw, wanted).map(Some),
        None => Ok(None),
    }
}

/// Parse a canvas extent (`size-x` / `size-y`), which must be strictly positive.
///
/// The extent never reaches a mark, so this check is not protecting the math — it is
/// protecting the reading. A canvas of zero or negative width means the document is not
/// shaped the way this module believes it is, and importing marks out of a file we have
/// evidently misread is worse than refusing it.
fn parse_canvas_extent(
    attributes: &[(&str, &str)],
    wanted: &'static str,
) -> Result<Option<f64>, ReticleError> {
    let Some(value) = parse_optional_angle_mil(attributes, wanted)? else {
        return Ok(None);
    };
    if value <= 0.0 {
        return Err(invalid(format!(
            "{wanted} must be greater than zero (got {value} mil)"
        )));
    }
    Ok(Some(value))
}

/// Convert one unit-suffixed attribute value to milliradians.
///
/// `.reticle` writes the unit inline on every value (`"96moa"`, `"-0.55moa"`, `"1.5mil"`),
/// so each attribute is converted on its own terms; there is no document-level unit to
/// fall back on. A value with no suffix or an unrecognized one is refused rather than
/// assigned a guessed scale.
fn parse_angle_mil(raw: &str, attribute: &str) -> Result<f64, ReticleError> {
    let text = raw.trim();

    // The unit is the trailing run of ASCII letters. Splitting from the right (rather than
    // scanning for the first letter) keeps scientific notation intact: "1e5moa" splits into
    // "1e5" and "moa", because the digit stops the run.
    let suffix_len = text
        .chars()
        .rev()
        .take_while(|c| c.is_ascii_alphabetic())
        .count();
    let (number, unit) = text.split_at(text.len() - suffix_len);

    let scale = if unit.eq_ignore_ascii_case("moa") {
        MOA_TO_MIL
    } else if unit.eq_ignore_ascii_case("mil")
        || unit.eq_ignore_ascii_case("mils")
        || unit.eq_ignore_ascii_case("mrad")
    {
        1.0
    } else if unit.is_empty() {
        return Err(invalid(format!(
            "{attribute}=\"{raw}\" has no unit suffix; every .reticle value names its own \
             unit (e.g. \"1.5moa\")"
        )));
    } else {
        return Err(invalid(format!(
            "{attribute}=\"{raw}\" uses the unrecognized unit \"{unit}\"; this importer \
             reads moa, mil and mrad"
        )));
    };

    let value: f64 = number.trim().parse().map_err(|_| {
        invalid(format!(
            "{attribute}=\"{raw}\" is not a number followed by a unit"
        ))
    })?;
    if !value.is_finite() {
        return Err(invalid(format!("{attribute}=\"{raw}\" is not finite")));
    }

    Ok(value * scale)
}

fn invalid(reason: String) -> ReticleError {
    ReticleError::InvalidSpec(reason)
}

/// Expand the five predefined XML entities and numeric character references.
///
/// Only ever applied to `name`, the one attribute that becomes a human-readable string.
/// An unrecognized entity is left verbatim rather than failing the import — a mangled
/// display name is not worth refusing a reticle over.
fn decode_entities(raw: &str) -> String {
    if !raw.contains('&') {
        return raw.to_string();
    }

    let mut out = String::with_capacity(raw.len());
    let mut rest = raw;
    while let Some(start) = rest.find('&') {
        out.push_str(&rest[..start]);
        let tail = &rest[start..];
        let Some(end) = tail.find(';') else {
            // No terminator at all: the remainder is literal text.
            out.push_str(tail);
            return out;
        };
        let entity = &tail[1..end];
        match entity {
            "amp" => out.push('&'),
            "lt" => out.push('<'),
            "gt" => out.push('>'),
            "quot" => out.push('"'),
            "apos" => out.push('\''),
            _ => {
                let decoded = entity.strip_prefix('#').and_then(|digits| {
                    let code = match digits.strip_prefix(['x', 'X']) {
                        Some(hex) => u32::from_str_radix(hex, 16).ok()?,
                        None => digits.parse::<u32>().ok()?,
                    };
                    char::from_u32(code)
                });
                match decoded {
                    Some(c) => out.push(c),
                    None => out.push_str(&tail[..=end]),
                }
            }
        }
        rest = &tail[end + 1..];
    }
    out.push_str(rest);
    out
}

// ---------------------------------------------------------------------------------------
// Minimal XML scanner.
//
// Hand-rolled rather than pulled from a crate, for the same reason `profile_import` hand-
// rolls its protobuf wire decoder: this crate ships to thirteen platforms plus wasm32, and
// a whole XML library is a large dependency to carry for one importer that needs elements
// and attributes and nothing else. Text nodes, CDATA, namespaces, DTDs and entity
// declarations are all skipped — a `.reticle` document puts every value in an attribute,
// so no element content is ever read.
//
// It is a scanner, not a validator: it enforces the well-formedness a misread would hide
// (balanced tags, quoted attribute values, one root) and ignores the rest.
// ---------------------------------------------------------------------------------------

/// One scanned markup event. Comments, processing instructions, doctypes, CDATA and text
/// never surface.
#[derive(Debug)]
enum XmlEvent<'a> {
    Start {
        tag: &'a str,
        attributes: Vec<(&'a str, &'a str)>,
        self_closing: bool,
    },
    End {
        tag: &'a str,
    },
}

struct XmlScanner<'a> {
    src: &'a str,
    pos: usize,
}

impl<'a> XmlScanner<'a> {
    fn new(src: &'a str) -> Self {
        Self { src, pos: 0 }
    }

    fn bytes(&self) -> &'a [u8] {
        self.src.as_bytes()
    }

    fn rest(&self) -> &'a str {
        &self.src[self.pos..]
    }

    /// Advance past `open`-prefixed markup terminated by `close`, or fail naming it.
    fn skip_delimited(&mut self, open: &str, close: &str, what: &str) -> Result<(), ReticleError> {
        let after_open = self.pos + open.len();
        match self.src[after_open..].find(close) {
            Some(offset) => {
                self.pos = after_open + offset + close.len();
                Ok(())
            }
            None => Err(invalid(format!(
                "unterminated {what} at byte {} (expected {close:?})",
                self.pos
            ))),
        }
    }

    /// The next start or end tag, or `None` at end of input.
    fn next_event(&mut self) -> Result<Option<XmlEvent<'a>>, ReticleError> {
        loop {
            // Everything between tags is text content, which this format never uses.
            match self.rest().find('<') {
                Some(offset) => self.pos += offset,
                None => {
                    self.pos = self.src.len();
                    return Ok(None);
                }
            }

            let rest = self.rest();
            if rest.starts_with("<!--") {
                self.skip_delimited("<!--", "-->", "comment")?;
            } else if rest.starts_with("<![CDATA[") {
                self.skip_delimited("<![CDATA[", "]]>", "CDATA section")?;
            } else if rest.starts_with("<?") {
                self.skip_delimited("<?", "?>", "processing instruction")?;
            } else if rest.starts_with("<!") {
                // DOCTYPE and friends. A document type declaration with an internal subset
                // would end at the wrong '>', which is precisely why nothing here reads
                // element content: the worst case is a spurious parse error, never a
                // silently different reticle.
                self.skip_delimited("<!", ">", "declaration")?;
            } else if rest.starts_with("</") {
                return self.read_end_tag().map(Some);
            } else {
                return self.read_start_tag().map(Some);
            }
        }
    }

    fn read_end_tag(&mut self) -> Result<XmlEvent<'a>, ReticleError> {
        self.pos += 2; // past "</"
        let start = self.pos;
        let Some(tag) = self.read_name() else {
            return Err(invalid(format!(
                "an end tag has an empty name at byte {start}"
            )));
        };
        self.skip_whitespace();
        if self.bytes().get(self.pos) != Some(&b'>') {
            return Err(invalid(format!(
                "</{tag}> is malformed at byte {}",
                self.pos
            )));
        }
        self.pos += 1;
        Ok(XmlEvent::End { tag })
    }

    fn read_start_tag(&mut self) -> Result<XmlEvent<'a>, ReticleError> {
        self.pos += 1; // past '<'
        let start = self.pos;
        let Some(tag) = self.read_name() else {
            return Err(invalid(format!(
                "a start tag has an empty name at byte {start}"
            )));
        };
        let mut attributes: Vec<(&'a str, &'a str)> = Vec::new();

        loop {
            self.skip_whitespace();
            let rest = self.rest();
            if rest.starts_with("/>") {
                self.pos += 2;
                return Ok(XmlEvent::Start {
                    tag,
                    attributes,
                    self_closing: true,
                });
            }
            if rest.starts_with('>') {
                self.pos += 1;
                return Ok(XmlEvent::Start {
                    tag,
                    attributes,
                    self_closing: false,
                });
            }
            if rest.is_empty() {
                return Err(invalid(format!("<{tag}> is never terminated")));
            }

            let Some(key) = self.read_name() else {
                return Err(invalid(format!(
                    "<{tag}> has a malformed attribute at byte {}",
                    self.pos
                )));
            };
            self.skip_whitespace();
            if self.bytes().get(self.pos) != Some(&b'=') {
                return Err(invalid(format!(
                    "attribute {key:?} of <{tag}> has no value"
                )));
            }
            self.pos += 1;
            self.skip_whitespace();

            let quote = match self.bytes().get(self.pos) {
                Some(&q @ (b'"' | b'\'')) => q,
                _ => {
                    return Err(invalid(format!(
                        "attribute {key:?} of <{tag}> must be quoted"
                    )))
                }
            };
            self.pos += 1;
            let start = self.pos;
            let Some(offset) = self.rest().find(quote as char) else {
                return Err(invalid(format!(
                    "attribute {key:?} of <{tag}> has an unterminated value"
                )));
            };
            self.pos = start + offset;
            attributes.push((key, &self.src[start..self.pos]));
            self.pos += 1; // past the closing quote
        }
    }

    /// Read an element or attribute name: everything up to whitespace or one of `/>=`.
    /// `None` when there is no name there at all, which every caller reports as malformed.
    fn read_name(&mut self) -> Option<&'a str> {
        let start = self.pos;
        let bytes = self.bytes();
        while let Some(&byte) = bytes.get(self.pos) {
            if byte.is_ascii_whitespace() || matches!(byte, b'/' | b'>' | b'=') {
                break;
            }
            self.pos += 1;
        }
        if self.pos == start {
            return None;
        }
        Some(&self.src[start..self.pos])
    }

    fn skip_whitespace(&mut self) {
        let bytes = self.bytes();
        while matches!(bytes.get(self.pos), Some(byte) if byte.is_ascii_whitespace()) {
            self.pos += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reticle::hold_point_in_reticle;

    /// A miniature but complete document in the shape a real one takes: a filled background
    /// rectangle, a field-of-view circle, crosshair arms, one windage hash, a filled post
    /// path, and a `<bdc>` section listing three holdovers that the drawing also draws.
    /// Invented content — the geometry is not any product's reticle.
    const LADDER_MOA: &str = r#"<?xml version="1.0" encoding="utf-8"?>
<reticle name="Test Ladder MOA" size-x="40moa" size-y="40moa" zero-x="20moa" zero-y="20moa">
  <elements>
    <reticle-rectangle position-x="-20moa" position-y="20moa" size-x="40moa" size-y="40moa" fill="true" color="black" />
    <reticle-circle center-x="0moa" center-y="0moa" radius="19moa" fill="true" line-width="0moa" color="white" />
    <!-- crosshair arms -->
    <reticle-line start-x="0moa" start-y="0moa" end-x="0moa" end-y="6moa" line-width="0.2moa" line-color="black" />
    <reticle-line start-x="0moa" start-y="0moa" end-x="0moa" end-y="-9moa" line-width="0.2moa" line-color="black" />
    <reticle-line start-x="-3moa" start-y="-0.5moa" end-x="-3moa" end-y="0.5moa" line-width="0.2moa" line-color="black" />
    <!-- the holdover ticks the bdc section names again -->
    <reticle-line start-x="-0.5moa" start-y="-2moa" end-x="0.5moa" end-y="-2moa" line-width="0.2moa" line-color="black" />
    <reticle-line start-x="-0.5moa" start-y="-5moa" end-x="0.5moa" end-y="-5moa" line-width="0.2moa" line-color="black" />
    <reticle-line start-x="-0.5moa" start-y="-8moa" end-x="0.5moa" end-y="-8moa" line-width="0.2moa" line-color="black" />
    <reticle-path fill="true" color="black">
      <elements>
        <reticle-path-move-to position-x="-0.3moa" position-y="19moa" />
        <reticle-path-line-to position-x="0.3moa" position-y="19moa" />
        <reticle-path-line-to position-x="0.3moa" position-y="6moa" />
        <reticle-path-line-to position-x="-0.3moa" position-y="6moa" />
      </elements>
    </reticle-path>
  </elements>
  <bdc>
    <bdc position-x="0moa" position-y="-2moa" text-offset="1moa" text-height="0.8moa" />
    <bdc position-x="0moa" position-y="-5moa" text-offset="1moa" text-height="0.8moa" />
    <bdc position-x="0moa" position-y="-8moa" text-offset="1moa" text-height="0.8moa" />
  </bdc>
</reticle>"#;

    /// The imported `down_mil` values, sorted.
    fn downs(description: &ReticleDescription) -> Vec<f64> {
        let mut downs: Vec<f64> = description.marks.iter().map(|m| m.down_mil).collect();
        downs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        downs
    }

    #[test]
    fn bdc_entries_are_the_only_holds_and_the_drawing_is_dropped() {
        let (description, report) = import_reticle_document_with_report(LADDER_MOA).unwrap();

        assert_eq!(description.name, "Test Ladder MOA");
        assert_eq!(
            description.marks.len(),
            3,
            "only the three <bdc> entries carry holds; the drawing draws six lines, a \
             circle, a rectangle and a path and none of them is a mark"
        );
        assert!(
            description.marks.iter().all(|m| m.kind == MarkKind::Hash),
            "a bdc point has no declared shape; every imported mark is a hash"
        );
        assert!(
            description.marks.iter().all(|m| m.label.is_none()),
            "a bdc entry carries no range text, so no mark is labeled"
        );

        // The three holdover ticks are drawn AND listed; importing both would double them.
        assert_eq!(
            report
                .dropped_element_tags
                .iter()
                .find(|(tag, _)| tag == "reticle-line")
                .map(|(_, count)| *count),
            Some(6),
            "the three drawn holdover ticks are dropped as drawing, not imported twice"
        );
        assert_eq!(
            report.dropped_element_tags,
            vec![
                ("reticle-circle".to_string(), 1),
                ("reticle-line".to_string(), 6),
                ("reticle-path".to_string(), 1),
                ("reticle-path-line-to".to_string(), 3),
                ("reticle-path-move-to".to_string(), 1),
                ("reticle-rectangle".to_string(), 1),
            ],
            "the report tallies every drawing tag it walked past, sorted by tag"
        );
        assert_eq!(report.drawing_elements_dropped, 13);
    }

    #[test]
    fn positive_y_in_the_document_is_down_negated_in_the_engine() {
        let description = import_reticle_document(LADDER_MOA).unwrap();
        // -2/-5/-8 MOA in the document (below the zero, +y up) become POSITIVE down_mil.
        let expected: Vec<f64> = vec![2.0, 5.0, 8.0]
            .into_iter()
            .map(|moa| moa * MOA_TO_MIL)
            .collect();
        for (got, want) in downs(&description).iter().zip(expected.iter()) {
            assert!(
                (got - want).abs() < 1e-12,
                "holdover {got} mil should be {want} mil below the zero"
            );
        }
        assert!(
            description.marks.iter().all(|m| m.down_mil > 0.0),
            "a holdover is BELOW the zero after the sign flip"
        );
    }

    #[test]
    fn a_mark_above_the_zero_imports_as_negative_down() {
        // A stadia mark above the zero (a positive position-y) is a negative down_mil.
        let description = import_reticle_document(
            r#"<reticle name="Up" size-x="10mil" size-y="10mil" zero-x="5mil" zero-y="5mil">
                 <bdc><bdc position-x="0mil" position-y="3mil" /></bdc>
               </reticle>"#,
        )
        .unwrap();
        assert_eq!(description.marks.len(), 1);
        assert_eq!(description.marks[0].down_mil, -3.0);
    }

    #[test]
    fn each_value_converts_on_its_own_unit_suffix() {
        // No document-level unit: mixed suffixes in one document all normalize to mil.
        let description = import_reticle_document(
            r#"<reticle name="Mixed" size-x="40moa" size-y="10mrad" zero-x="20moa" zero-y="5mrad">
                 <bdc>
                   <bdc position-x="0mil" position-y="-1.5mil" />
                   <bdc position-x="0mrad" position-y="-3mrad" />
                   <bdc position-x="0moa" position-y="-2MOA" />
                 </bdc>
               </reticle>"#,
        )
        .unwrap();
        // Document order, not sorted: each entry is converted on its own suffix.
        let got: Vec<f64> = description.marks.iter().map(|m| m.down_mil).collect();
        assert_eq!(got[0], 1.5, "mil passes through");
        assert_eq!(got[1], 3.0, "mrad passes through");
        assert!(
            (got[2] - 2.0 * MOA_TO_MIL).abs() < 1e-12,
            "MOA is case-insensitive and scales by MOA_TO_MIL (got {})",
            got[2]
        );
    }

    #[test]
    fn a_value_this_module_cannot_scale_is_refused_rather_than_guessed() {
        for (xml, why) in [
            (
                r#"<reticle><bdc><bdc position-y="-2" /></bdc></reticle>"#,
                "no unit suffix",
            ),
            (
                r#"<reticle><bdc><bdc position-y="-2cm100m" /></bdc></reticle>"#,
                "unrecognized unit",
            ),
            (
                r#"<reticle><bdc><bdc position-y="twomoa" /></bdc></reticle>"#,
                "not a number",
            ),
            (
                r#"<reticle><bdc><bdc position-y="1e400moa" /></bdc></reticle>"#,
                "overflows to infinity",
            ),
        ] {
            assert!(
                matches!(
                    import_reticle_document(xml),
                    Err(ReticleError::InvalidSpec(_))
                ),
                "{why} must be an InvalidSpec, not a guessed scale"
            );
        }
    }

    #[test]
    fn a_non_positive_canvas_is_refused() {
        assert!(matches!(
            import_reticle_document(r#"<reticle size-x="0moa" size-y="10moa"></reticle>"#),
            Err(ReticleError::InvalidSpec(_))
        ));
        assert!(matches!(
            import_reticle_document(r#"<reticle size-x="-4moa" size-y="10moa"></reticle>"#),
            Err(ReticleError::InvalidSpec(_))
        ));
    }

    #[test]
    fn the_zero_attributes_place_the_canvas_and_never_move_a_mark() {
        // Two documents identical but for zero-x/zero-y: the zero is the origin of the
        // element coordinate system, so re-framing the canvas cannot move a hold.
        let centered = import_reticle_document(
            r#"<reticle name="Z" size-x="40moa" size-y="40moa" zero-x="20moa" zero-y="20moa">
                 <bdc><bdc position-x="1moa" position-y="-5moa" /></bdc>
               </reticle>"#,
        )
        .unwrap();
        let offset = import_reticle_document(
            r#"<reticle name="Z" size-x="40moa" size-y="40moa" zero-x="12moa" zero-y="31moa">
                 <bdc><bdc position-x="1moa" position-y="-5moa" /></bdc>
               </reticle>"#,
        )
        .unwrap();
        assert_eq!(
            centered, offset,
            "zero-x/zero-y are canvas placement; they must not offset the marks"
        );
    }

    #[test]
    fn the_canvas_and_zero_are_reported_in_milliradians() {
        let (_, report) = import_reticle_document_with_report(LADDER_MOA).unwrap();
        let (size_x, size_y) = report.canvas_mil.expect("canvas declared");
        let (zero_x, zero_y) = report.zero_in_canvas_mil.expect("zero declared");
        assert!((size_x - 40.0 * MOA_TO_MIL).abs() < 1e-12);
        assert!((size_y - 40.0 * MOA_TO_MIL).abs() < 1e-12);
        assert!((zero_x - 20.0 * MOA_TO_MIL).abs() < 1e-12);
        assert!((zero_y - 20.0 * MOA_TO_MIL).abs() < 1e-12);

        // Undeclared canvas metadata is absent, not defaulted.
        let (_, bare) = import_reticle_document_with_report(
            r#"<reticle name="Bare"><bdc><bdc position-y="-1mil" /></bdc></reticle>"#,
        )
        .unwrap();
        assert_eq!(bare.canvas_mil, None);
        assert_eq!(bare.zero_in_canvas_mil, None);
    }

    #[test]
    fn a_bdc_wrapper_is_not_itself_a_mark() {
        // The section wrapper carries no position; only its children are entries. An entry
        // naming one axis sits at 0 on the other.
        let description = import_reticle_document(
            r#"<reticle name="W" size-x="10mil" size-y="10mil">
                 <bdc>
                   <bdc position-y="-1mil" />
                   <bdc position-x="2mil" />
                 </bdc>
               </reticle>"#,
        )
        .unwrap();
        assert_eq!(description.marks.len(), 2, "the wrapper is not a mark");
        assert_eq!(description.marks[0].down_mil, 1.0);
        assert_eq!(description.marks[0].right_mil, 0.0);
        assert_eq!(description.marks[1].down_mil, 0.0);
        assert_eq!(description.marks[1].right_mil, 2.0);

        // An empty section contributes nothing at all.
        let empty = import_reticle_document(r#"<reticle name="E"><bdc /></reticle>"#).unwrap();
        assert!(empty.marks.is_empty());
    }

    #[test]
    fn an_all_decoration_document_imports_with_no_marks() {
        // A pure transform: a document that draws a reticle but declares no bdc entry is
        // not an error here, it is simply a description with nothing to hold on. The
        // rejection belongs to the solver.
        let description = import_reticle_document(
            r#"<reticle name="Picture" size-x="20moa" size-y="20moa">
                 <elements>
                   <reticle-line start-x="-8moa" start-y="0moa" end-x="8moa" end-y="0moa" line-width="0.2moa" />
                   <reticle-circle center-x="0moa" center-y="0moa" radius="2moa" />
                   <future-shape whatever="7pizzas" />
                 </elements>
               </reticle>"#,
        )
        .unwrap();
        assert!(description.marks.is_empty());
        assert_eq!(description.validate(), Err(ReticleError::NoMarks));
        assert!(
            matches!(
                hold_point_in_reticle(1.0, 0.0, 1.0, &description),
                Err(ReticleError::NoMarks)
            ),
            "an all-decoration reticle has nothing to hold on"
        );
    }

    #[test]
    fn dropped_geometry_is_never_read_as_an_angle() {
        // Decoration geometry is tallied, not read: an unreadable unit inside a shape we
        // drop must not fail the import, exactly as the Ventum importer ignores the fields
        // of its decoration variants.
        let (description, report) = import_reticle_document_with_report(
            r#"<reticle name="Odd" size-x="20moa" size-y="20moa">
                 <elements>
                   <reticle-line start-x="nonsense" end-y="12furlongs" />
                 </elements>
                 <bdc><bdc position-y="-1moa" /></bdc>
               </reticle>"#,
        )
        .unwrap();
        assert_eq!(description.marks.len(), 1);
        assert_eq!(report.drawing_elements_dropped, 1);
    }

    #[test]
    fn the_import_is_always_ffp_identity_because_the_format_says_nothing() {
        // The name may talk about magnification in prose; the importer does not mine it.
        let description = import_reticle_document(
            r#"<reticle name="Scope 4-12x40, at 12x" size-x="20moa" size-y="20moa">
                 <bdc><bdc position-y="-2moa" /></bdc>
               </reticle>"#,
        )
        .unwrap();
        assert_eq!(description.focal_plane, FocalPlane::First);
        assert_eq!(description.reference_magnification, 1.0);
        assert_eq!(
            description.name, "Scope 4-12x40, at 12x",
            "the magnification hint stays prose in the name"
        );
    }

    #[test]
    fn malformed_documents_report_invalid_spec() {
        for (xml, why) in [
            ("<reticle>", "unclosed root"),
            ("<reticle><bdc></reticle>", "mismatched close tag"),
            ("</reticle>", "close with nothing open"),
            ("", "no root at all"),
            ("   \n  ", "whitespace only"),
            ("<scope name=\"x\" />", "root is not <reticle>"),
            ("<reticle /><reticle />", "two roots"),
            (r#"<reticle name=x />"#, "unquoted attribute value"),
            (r#"<reticle name />"#, "attribute with no value"),
            (r#"<reticle name="unterminated />"#, "unterminated value"),
            ("<reticle", "unterminated start tag"),
            ("<!-- never ends", "unterminated comment"),
            ("<>", "empty tag name"),
        ] {
            assert!(
                matches!(
                    import_reticle_document(xml),
                    Err(ReticleError::InvalidSpec(_))
                ),
                "{why} must be reported as an InvalidSpec, not accepted"
            );
        }
    }

    #[test]
    fn prologue_comments_and_cdata_are_skipped() {
        let description = import_reticle_document(
            r#"<?xml version="1.0"?>
               <!DOCTYPE reticle>
               <!-- a comment mentioning <reticle-line /> and </reticle> -->
               <reticle name="Skip" size-x="10mil" size-y="10mil">
                 <![CDATA[ <bdc position-y="-99mil" /> ]]>
                 <bdc><bdc position-y="-1mil" /></bdc>
               </reticle>"#,
        )
        .unwrap();
        assert_eq!(
            description.marks.len(),
            1,
            "markup inside a comment or CDATA is text, not an element"
        );
        assert_eq!(description.marks[0].down_mil, 1.0);
    }

    #[test]
    fn too_many_bdc_entries_are_capped_during_collection() {
        let mut xml = String::from(r#"<reticle name="Flood"><bdc>"#);
        for _ in 0..(MAX_RETICLE_MARKS + 10) {
            xml.push_str(r#"<bdc position-y="-1mil" />"#);
        }
        xml.push_str("</bdc></reticle>");
        assert!(matches!(
            import_reticle_document(&xml),
            Err(ReticleError::TooManyMarks { max, .. }) if max == MAX_RETICLE_MARKS
        ));
    }

    #[test]
    fn nesting_deeper_than_the_cap_is_refused_without_recursing() {
        let mut xml = String::from("<reticle>");
        for _ in 0..(MAX_DOCUMENT_DEPTH + 5) {
            xml.push_str("<elements>");
        }
        assert!(matches!(
            import_reticle_document(&xml),
            Err(ReticleError::InvalidSpec(_))
        ));
    }

    #[test]
    fn entities_in_the_name_are_decoded() {
        let description = import_reticle_document(
            r#"<reticle name="Duplex &amp; BDC &quot;Pro&quot; &#8212; &#x33;x" />"#,
        )
        .unwrap();
        assert_eq!(description.name, r#"Duplex & BDC "Pro" — 3x"#);

        // An entity this module does not know stays verbatim rather than failing the import.
        let odd = import_reticle_document(r#"<reticle name="A &nbsp; B &broken C" />"#).unwrap();
        assert_eq!(odd.name, "A &nbsp; B &broken C");
    }

    #[test]
    fn an_imported_ladder_hold_solves() {
        let description = import_reticle_document(LADDER_MOA).unwrap();
        // A 5 MOA drop lands exactly on the middle holdover.
        let hold = hold_point_in_reticle(5.0 * MOA_TO_MIL, 0.0, 1.0, &description).unwrap();
        assert!(hold.nearest_mark.is_some());
        assert!(
            hold.nearest_mark_distance_mil < 1e-9,
            "the hold sits on the mark (distance {})",
            hold.nearest_mark_distance_mil
        );
        assert!(!hold.off_reticle);

        // And the documented cost of the hold/decoration rule, asserted rather than
        // described: the drawn windage hash was dropped, so the windage axis is degenerate
        // and any wind hold reads off-reticle.
        let windy = hold_point_in_reticle(5.0 * MOA_TO_MIL, 0.5, 1.0, &description).unwrap();
        assert!(
            windy.off_reticle,
            "with elevation-only marks, any wind hold is off-reticle by construction"
        );
    }
}
