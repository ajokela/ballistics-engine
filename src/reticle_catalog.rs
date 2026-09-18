//! Named reticles, by id — the tracked follow-up [`crate::reticle`]'s header defers to.
//!
//! `reticle` deliberately has no catalog: it ships three generic constructors and says so,
//! because "manufacturer subtension sheets are published facts and are a legally viable
//! catalog source, but curating one is a separate, per-vendor IP-reviewed data project".
//! This module is that project, kept separate so that statement stays true of the module
//! it was made about.
//!
//! # Why a catalog is the whole feature (MBA-1545)
//!
//! A holdover display is only worth anything if it shows the shooter THEIR reticle. A hold
//! against a generic mil grid reads "hold 4.2 mil down", which is what a range card
//! already prints — the same number, worse. The value is "hold on the third dot down, two
//! right", and that needs the marks actually etched in their glass.
//!
//! Nothing public supplies those. Surveyed 2026-09-17: the `.reticle` XML format
//! (gehtsoft-usa/BallisticCalculator1, LGPL-2.1) ships a designer and worked examples, not
//! commercial designs; `o-murphy/bcrg` (LGPL-3.0) generates reticles procedurally and emits
//! 1-bit BMPs, which carry no coordinates and cannot drive hold math; and no
//! machine-readable dataset of real reticles appears to exist. The real source is
//! manufacturer subtension sheets, transcribed by hand.
//!
//! # The two constraints, which are different constraints
//!
//! **Licensing.** Both open projects in this space are LGPL and this crate is MIT OR
//! Apache-2.0. Every entry here must be OUR transcription of published subtension facts —
//! measurements, which are not copyrightable — and never a file lifted from an LGPL
//! project. This is the same reason `reticle_import` and `profile_import` are cleanroom
//! implementations written from format specifications.
//!
//! **Patents.** Horus grid reticles and Time-of-Flight Wind Dots are actively patented.
//! That bars reproducing those LAYOUTS however the numbers were obtained, and it is not a
//! data-licensing question — a patented layout stays barred even where its subtensions are
//! published. No TREMOR-family or Horus geometry enters this module, on the same terms
//! [`crate::reticle`] states for itself.
//!
//! # Provenance is mandatory, because a wrong subtension is a silently wrong hold
//!
//! Every [`CatalogEntry`] carries a [`CatalogEntry::source`] naming where its geometry came
//! from. That is not decoration. A mistyped subtension produces a hold that is wrong by a
//! mark and looks entirely normal — there is no downstream check that can catch it, the
//! same property that makes a grams-for-grains bullet weight invisible (MBA-1521). Naming
//! the source per entry is what makes a bad number traceable to a document someone can
//! re-read, rather than anonymous.
//!
//! So this module starts with what can be stated EXACTLY from a public standard, and each
//! vendor reticle added later is its own change citing its own sheet.

use crate::reticle::{FocalPlane, MarkKind, ReticleDescription, ReticleError, ReticleMark};

/// One named reticle this build can produce, and where its geometry came from.
#[derive(Debug, Clone, PartialEq)]
pub struct CatalogEntry {
    /// Stable wire id — lowercase, hyphenated. Callers persist this, so it does not change.
    pub id: &'static str,
    /// What a picker shows.
    pub display_name: &'static str,
    /// Where the subtensions came from, specifically enough to re-check. Mandatory; see
    /// the module header.
    pub source: &'static str,
    /// What a shooter should know before trusting it — variations between scopes that
    /// carry this reticle, and anything the model does not capture.
    pub notes: &'static str,
}

/// The classic mil-dot, at its standard subtensions.
///
/// Dots centred every **1 mil** along both stadia, each dot subtending **0.2 mil**. That
/// spacing is the standard and is what a hold is read against; it is why the reticle is
/// called what it is, and it does not vary between the scopes that carry it.
///
/// `dots_per_side` is how many dots are etched out from centre on each arm, and it DOES
/// vary — four and five are both common, and some scopes suppress the innermost dot. It is
/// a parameter rather than a constant for that reason, and it changes nothing about a hold
/// that falls inside the ladder: the marks are at the same angles either way, there are
/// just more or fewer of them.
///
/// The dot DIAMETER is carried in the docs rather than the model because
/// [`ReticleDescription`] positions marks and does not size them — a hold is measured to a
/// mark's centre. A 0.2 mil dot is worth knowing when reading a hold off glass (it is its
/// own coarse ranging reference) and changes no arithmetic here.
///
/// FFP, with `reference_magnification` 1.0, which is unused on that plane. A mil-dot etched
/// in the second focal plane is the same geometry at its own reference magnification;
/// callers set those two fields, exactly as they do for a generated reticle.
pub fn mil_dot(dots_per_side: usize) -> Result<ReticleDescription, ReticleError> {
    if dots_per_side == 0 {
        return Err(ReticleError::InvalidGeneratorParameter {
            parameter: "dots_per_side",
            value: 0.0,
            rule: "at least one dot on each arm",
        });
    }

    let mut marks = Vec::with_capacity(1 + dots_per_side * 4);
    marks.push(ReticleMark::new(0.0, 0.0, MarkKind::Center));
    for n in 1..=dots_per_side {
        let mil = n as f64;
        // Down and up the vertical stadium, then right and left the horizontal one. Labels
        // name the angle, which is the only thing a shooter reads a mil-dot for.
        marks.push(ReticleMark::labeled(
            mil,
            0.0,
            MarkKind::Dot,
            format!("{n} mil down"),
        ));
        marks.push(ReticleMark::labeled(
            -mil,
            0.0,
            MarkKind::Dot,
            format!("{n} mil up"),
        ));
        marks.push(ReticleMark::labeled(
            0.0,
            mil,
            MarkKind::Dot,
            format!("{n} mil right"),
        ));
        marks.push(ReticleMark::labeled(
            0.0,
            -mil,
            MarkKind::Dot,
            format!("{n} mil left"),
        ));
    }

    Ok(ReticleDescription {
        name: format!("Mil-Dot ({dots_per_side} per side)"),
        focal_plane: FocalPlane::First,
        reference_magnification: 1.0,
        marks,
    })
}

/// How many dots per arm [`MIL_DOT_ID`] builds when no count is asked for.
///
/// Five, which reaches 5 mil — far enough for the holds this is used for, and a count real
/// scopes carry. See [`mil_dot`] on why this is a convention and the 1 mil spacing is not.
pub const MIL_DOT_DEFAULT_DOTS_PER_SIDE: usize = 5;

/// Wire id of the mil-dot entry.
pub const MIL_DOT_ID: &str = "mil-dot";

/// A plain hash ladder: marks every `spacing` along both stadia out to `extent`.
///
/// GENERIC GEOMETRY, and that is the whole point of these entries. A hash ladder
/// is not anybody's design — it is what a reticle looks like when it is just a
/// ruler, and a great many scopes carry exactly that under a vendor's name. So
/// these can be stated exactly, from the spacing alone, without transcribing
/// anyone's subtension sheet and without reproducing anyone's layout.
///
/// `unit_label` names the angular unit the spacing is quoted in, for the display
/// name only. The geometry is milliradians either way, because
/// [`ReticleMark`] is milliradians — an MOA ladder is a ladder whose marks
/// happen to fall on MOA multiples, not a different kind of object.
fn hash_ladder(
    name: String,
    spacing_mil: f64,
    extent_mil: f64,
) -> Result<ReticleDescription, ReticleError> {
    let mut description = ReticleDescription::mil_grid(spacing_mil, extent_mil)?;
    description.name = name;
    Ok(description)
}

/// One MOA in milliradians, on the engine's locked printed-table constant.
///
/// `adjustment.rs` pins MOA at 3438 rather than the exact-angle 3437.7467 and
/// prints every MOA column on it (MBA-724), so an MOA ladder built here has to
/// use the same number or its marks sit a hair off the dial column beside them.
/// 1 MOA = 1000/3438 mil.
const MIL_PER_MOA: f64 = 1000.0 / 3438.0;

/// Wire id of the half-mil hash ladder.
pub const MIL_HASH_HALF_ID: &str = "mil-hash-0.5";

/// Wire id of the fine two-tenth mil hash ladder.
pub const MIL_HASH_FINE_ID: &str = "mil-hash-0.2";

/// Wire id of the one-MOA hash ladder.
pub const MOA_HASH_1_ID: &str = "moa-hash-1";

/// Wire id of the two-MOA hash ladder.
pub const MOA_HASH_2_ID: &str = "moa-hash-2";

/// Every named reticle in this build.
///
/// Deliberately short. It holds what can be stated EXACTLY from a public standard; vendor
/// reticles arrive one change at a time, each citing the sheet it was transcribed from, so
/// that a bad subtension is traceable. An entry whose `source` would have to read "from
/// memory" does not belong here at all — see the module header.
pub fn catalog() -> Vec<CatalogEntry> {
    vec![
        CatalogEntry {
            id: MIL_DOT_ID,
            display_name: "Mil-Dot",
            source: "US military mil-dot standard: dots on 1 mil centres, 0.2 mil dot \
                 subtension. Long-published and carried unchanged by every scope that \
                 names the reticle; no vendor sheet is involved.",
            notes: "The 1 mil spacing is the standard and does not vary. How many dots are \
                etched on each arm DOES vary between scopes (four and five are both \
                common); this builds five, and a hold inside the ladder reads the same \
                either way. Dots subtend 0.2 mil, which this model does not draw — marks \
                are positions, and a hold is measured to a mark's centre.",
        },
        CatalogEntry {
            id: MIL_HASH_HALF_ID,
            display_name: "Mil hash, 0.5 mil",
            source: "Generic geometry, not a vendor design: hash marks on 0.5 mil \
                 centres along both stadia. Stated from the spacing alone — \
                 nothing here is transcribed from anyone's subtension sheet.",
            notes: "Many scopes carry a plain half-mil ladder under a vendor name. \
                Use this when yours is simply a ruler; if it has a tree, \
                irregular windage marks or labelled BDC marks, it is not this \
                and the holds will name the wrong mark.",
        },
        CatalogEntry {
            id: MIL_HASH_FINE_ID,
            display_name: "Mil hash, 0.2 mil",
            source: "Generic geometry, not a vendor design: hash marks on 0.2 mil \
                 centres along both stadia.",
            notes: "A fine ladder, typical of a first-focal-plane target reticle at \
                high magnification. Same caveat as the half-mil ladder: this is \
                a ruler, not a tree.",
        },
        CatalogEntry {
            id: MOA_HASH_1_ID,
            display_name: "MOA hash, 1 MOA",
            source: "Generic geometry, not a vendor design: hash marks on 1 MOA \
                 centres along both stadia, converted at the engine's locked \
                 printed-table 3438 MOA per radian (adjustment.rs, MBA-724) so \
                 the marks agree with the MOA column beside them.",
            notes: "An MOA ladder is a ladder whose marks fall on MOA multiples; \
                the model is milliradians either way. 1 MOA is about 0.29 mil, \
                so these sit closer together than a mil ladder's.",
        },
        CatalogEntry {
            id: MOA_HASH_2_ID,
            display_name: "MOA hash, 2 MOA",
            source: "Generic geometry, not a vendor design: hash marks on 2 MOA \
                 centres along both stadia, on the same locked 3438 constant.",
            notes: "The coarser of the two MOA ladders, and the commoner one on a \
                hunting scope.",
        },
    ]
}

/// Build the reticle named by `id`, or `None` when this build has no such entry.
///
/// `None` rather than a default: a caller asking for a specific reticle and silently
/// getting a different one would be shown holds for glass they are not looking through.
pub fn by_id(id: &str) -> Option<Result<ReticleDescription, ReticleError>> {
    match id {
        MIL_DOT_ID => Some(mil_dot(MIL_DOT_DEFAULT_DOTS_PER_SIDE)),
        // Ten mils of ladder, which reaches past any hold a shooter takes on a
        // reticle rather than a dial, and matches what mil_grid's own extent
        // means: marks out TO this, not beyond it.
        MIL_HASH_HALF_ID => Some(hash_ladder("Mil hash 0.5".into(), 0.5, 10.0)),
        MIL_HASH_FINE_ID => Some(hash_ladder("Mil hash 0.2".into(), 0.2, 5.0)),
        MOA_HASH_1_ID => Some(hash_ladder(
            "MOA hash 1".into(),
            MIL_PER_MOA,
            MIL_PER_MOA * 30.0,
        )),
        MOA_HASH_2_ID => Some(hash_ladder(
            "MOA hash 2".into(),
            MIL_PER_MOA * 2.0,
            MIL_PER_MOA * 30.0,
        )),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reticle::hold_point_in_reticle;

    #[test]
    fn mil_dot_puts_its_dots_on_one_mil_centres() {
        // The one thing about this reticle that is a standard rather than a convention.
        let reticle = mil_dot(5).unwrap();
        let mut down: Vec<f64> = reticle
            .marks
            .iter()
            .filter(|m| m.right_mil == 0.0 && m.down_mil > 0.0)
            .map(|m| m.down_mil)
            .collect();
        down.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(down, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    }

    #[test]
    fn it_is_a_cross_not_a_grid() {
        // The IP line `reticle`'s header draws, restated where a new module could cross it:
        // marks live on the two stadia and nowhere else. A filled two-dimensional grid is
        // patented geometry.
        let reticle = mil_dot(5).unwrap();
        for mark in &reticle.marks {
            assert!(
                mark.down_mil == 0.0 || mark.right_mil == 0.0,
                "mark at ({}, {}) is off both stadia — that is a grid",
                mark.down_mil,
                mark.right_mil
            );
        }
        assert_eq!(reticle.marks.len(), 1 + 5 * 4);
    }

    #[test]
    fn a_four_dot_scope_holds_the_same_as_a_five_dot_one() {
        // Why dots_per_side is a parameter and not a correctness question: the marks are
        // at the same angles, there are just fewer of them.
        let four = mil_dot(4).unwrap();
        let five = mil_dot(5).unwrap();
        for reticle in [&four, &five] {
            let hold = hold_point_in_reticle(3.0, 0.0, 10.0, reticle).unwrap();
            let mark = &reticle.marks[hold.nearest_mark.unwrap()];
            assert_eq!(mark.down_mil, 3.0);
            assert_eq!(hold.nearest_mark_distance_mil, 0.0);
        }
    }

    #[test]
    fn every_entry_names_a_source_and_builds() {
        // Provenance is the module's rule; a test is what makes it one. An entry that
        // cannot say where its numbers came from is exactly the entry that will be wrong
        // and untraceable.
        for entry in catalog() {
            assert!(
                !entry.source.trim().is_empty(),
                "{} names no source",
                entry.id
            );
            assert!(
                !entry.notes.trim().is_empty(),
                "{} carries no notes",
                entry.id
            );
            assert!(!entry.display_name.trim().is_empty());
            assert_eq!(entry.id, entry.id.to_lowercase(), "ids are lowercase");
            let built = by_id(entry.id)
                .unwrap_or_else(|| panic!("{} is listed but by_id does not build it", entry.id))
                .unwrap_or_else(|e| panic!("{} does not build: {e}", entry.id));
            built
                .validate()
                .unwrap_or_else(|e| panic!("{} builds an invalid reticle: {e}", entry.id));
        }
    }

    #[test]
    fn a_moa_ladder_falls_on_moa_multiples() {
        // An MOA ladder is a ladder whose marks fall on MOA multiples; the model
        // is milliradians either way. The constant is the engine's LOCKED 3438,
        // not the exact-angle 3437.7467, so these marks agree with the MOA
        // column printed beside them (MBA-724).
        let reticle = by_id(MOA_HASH_1_ID).unwrap().unwrap();
        let one_moa = 1000.0 / 3438.0;
        let mut down: Vec<f64> = reticle
            .marks
            .iter()
            .filter(|m| m.right_mil == 0.0 && m.down_mil > 0.0)
            .map(|m| m.down_mil)
            .collect();
        down.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (index, mil) in down.iter().enumerate() {
            let expected = one_moa * (index + 1) as f64;
            assert!(
                (mil - expected).abs() < 1e-9,
                "mark {index} at {mil} mil is not {} MOA",
                index + 1
            );
        }
        assert!(
            down.len() >= 20,
            "a 30 MOA ladder should carry 30 marks a side"
        );
    }

    #[test]
    fn every_ladder_is_a_cross_and_carries_its_centre() {
        // The IP line, restated for every entry rather than only the mil-dot: a
        // filled two-dimensional grid is patented geometry, and nothing here
        // may drift into one.
        for entry in catalog() {
            let reticle = by_id(entry.id).unwrap().unwrap();
            assert!(
                reticle
                    .marks
                    .iter()
                    .any(|m| m.kind == MarkKind::Center && m.down_mil == 0.0 && m.right_mil == 0.0),
                "{} has no centre",
                entry.id
            );
            for mark in &reticle.marks {
                assert!(
                    mark.down_mil == 0.0 || mark.right_mil == 0.0,
                    "{}: mark at ({}, {}) is off both stadia — that is a grid",
                    entry.id,
                    mark.down_mil,
                    mark.right_mil
                );
            }
        }
    }

    #[test]
    fn the_generic_ladders_say_they_are_generic() {
        // These are stated from a spacing, not transcribed from a sheet, and
        // the entry has to say so — a shooter comparing this against their own
        // glass needs to know it is a ruler and not their vendor's design.
        for id in [
            MIL_HASH_HALF_ID,
            MIL_HASH_FINE_ID,
            MOA_HASH_1_ID,
            MOA_HASH_2_ID,
        ] {
            let entry = catalog().into_iter().find(|e| e.id == id).unwrap();
            assert!(
                entry
                    .source
                    .contains("Generic geometry, not a vendor design"),
                "{id} does not say it is generic"
            );
        }
    }

    #[test]
    fn an_unknown_id_is_none_rather_than_a_substitute() {
        // Silently handing back a different reticle would show a shooter holds for glass
        // they are not looking through.
        assert!(by_id("no-such-reticle").is_none());
        assert!(by_id("").is_none());
        assert!(
            by_id("MIL-DOT").is_none(),
            "ids are case-sensitive on the wire"
        );
    }

    #[test]
    fn zero_dots_is_refused() {
        assert!(mil_dot(0).is_err());
    }
}
