//! End-to-end cover for the `.reticle` document importer against real files on disk.
//!
//! The unit tests in `src/reticle_document_import.rs` pin the transform against inline
//! strings; this harness runs the same importer over hand-written `.reticle` documents in
//! `tests/fixtures/reticle_document/`, so the format is exercised as a *file* — prologue,
//! comments, tabs, single-quoted attributes and all — and the hold/decoration rule is
//! asserted against a document a shooter would recognize.
//!
//! Every fixture is invented content. No third-party `.reticle` document is vendored here.

use ballistics_engine::reticle::{
    hold_point_in_reticle, FocalPlane, MarkKind, ReticleError, MAX_RETICLE_MARKS,
};
use ballistics_engine::reticle_document_import::{
    import_reticle_document, import_reticle_document_with_report,
};
use ballistics_engine::reticle_import::MOA_TO_MIL;

const DUPLEX_BDC_MOA: &str = include_str!("fixtures/reticle_document/duplex_bdc_moa.reticle");
const OFFSET_ZERO_MOA: &str = include_str!("fixtures/reticle_document/offset_zero_moa.reticle");
const WIND_GRID_MIL: &str = include_str!("fixtures/reticle_document/wind_grid_mil.reticle");
const DECORATION_ONLY: &str = include_str!("fixtures/reticle_document/decoration_only.reticle");
const MALFORMED_UNCLOSED: &str =
    include_str!("fixtures/reticle_document/malformed_unclosed.reticle");
const UNITLESS_VALUE: &str = include_str!("fixtures/reticle_document/unitless_value.reticle");

/// A duplex BDC document imports its three declared holdovers and nothing else, even
/// though the drawing contains fifteen shapes — four of which sit at aiming coordinates.
#[test]
fn a_duplex_bdc_document_imports_only_its_declared_holdovers() {
    let (description, report) = import_reticle_document_with_report(DUPLEX_BDC_MOA).unwrap();

    assert_eq!(description.name, "Fixture Duplex BDC MOA");
    assert_eq!(description.marks.len(), 3);
    assert!(description.marks.iter().all(|m| m.kind == MarkKind::Hash));
    assert!(description.marks.iter().all(|m| m.label.is_none()));

    let downs: Vec<f64> = description.marks.iter().map(|m| m.down_mil).collect();
    for (got, moa) in downs.iter().zip([2.0_f64, 5.0, 9.0]) {
        assert!(
            (got - moa * MOA_TO_MIL).abs() < 1e-12,
            "{moa} MOA below the zero should import as {} mil, got {got}",
            moa * MOA_TO_MIL
        );
    }
    assert!(
        description.marks.iter().all(|m| m.right_mil == 0.0),
        "this document's bdc ladder is on the vertical stadia"
    );

    // The drawing is dropped, and countably so. Four of the dropped lines are the four
    // windage hashes; three more are the very ticks the bdc section names again.
    assert_eq!(report.drawing_elements_dropped, 18);
    assert_eq!(
        report.dropped_element_tags,
        vec![
            ("reticle-circle".to_string(), 1),
            ("reticle-line".to_string(), 11),
            ("reticle-path".to_string(), 1),
            ("reticle-path-line-to".to_string(), 3),
            ("reticle-path-move-to".to_string(), 1),
            ("reticle-rectangle".to_string(), 1),
        ]
    );
    assert!(
        !report
            .dropped_element_tags
            .iter()
            .any(|(tag, _)| tag == "bdc" || tag == "elements"),
        "the bdc section and the structural containers are not drawing"
    );

    // Canvas metadata is reported in mil, never applied to a mark.
    let (size_x, size_y) = report.canvas_mil.expect("canvas declared");
    assert!((size_x - 60.0 * MOA_TO_MIL).abs() < 1e-12);
    assert!((size_y - 60.0 * MOA_TO_MIL).abs() < 1e-12);
    assert_eq!(
        report.zero_in_canvas_mil,
        Some((30.0 * MOA_TO_MIL, 30.0 * MOA_TO_MIL))
    );
}

/// The format carries no focal plane, so an import is always the identity. A caller with an
/// SFP optic has to say so; the fixture's own name only hints at it in prose.
#[test]
fn an_imported_document_is_ffp_identity_until_a_caller_says_otherwise() {
    let mut description = import_reticle_document(DUPLEX_BDC_MOA).unwrap();
    assert_eq!(description.focal_plane, FocalPlane::First);
    assert_eq!(description.reference_magnification, 1.0);

    // The documented one-line fix at the call site, and proof it takes effect: at half the
    // reference magnification an SFP mark subtends twice its nominal value.
    description.focal_plane = FocalPlane::Second;
    description.reference_magnification = 12.0;
    let hold = hold_point_in_reticle(2.0 * 2.0 * MOA_TO_MIL, 0.0, 6.0, &description).unwrap();
    assert_eq!(hold.mark_scale, 2.0);
    assert!(hold.nearest_mark_distance_mil < 1e-9);
}

/// Re-framing the canvas cannot move a hold: `zero-x`/`zero-y` place the drawing origin on
/// the output surface, and element coordinates are already measured from that origin.
#[test]
fn the_zero_attributes_do_not_offset_the_marks() {
    let centered = import_reticle_document(DUPLEX_BDC_MOA).unwrap();
    let (offset, report) = import_reticle_document_with_report(OFFSET_ZERO_MOA).unwrap();

    assert_eq!(
        centered.marks, offset.marks,
        "the two documents declare the same bdc entries, so they import the same marks"
    );
    assert_eq!(
        report.zero_in_canvas_mil,
        Some((19.0 * MOA_TO_MIL, 41.0 * MOA_TO_MIL)),
        "the off-centre zero is reported, not applied"
    );
}

/// A document whose `<bdc>` entries carry a `position-x` gets marks on both axes: the
/// importer does not assume a BDC ladder lives on the vertical stadia.
#[test]
fn bdc_entries_off_the_vertical_stadia_import_as_windage() {
    let description = import_reticle_document(WIND_GRID_MIL).unwrap();
    assert_eq!(description.marks.len(), 9, "three rows of three");

    let mut rights: Vec<f64> = description.marks.iter().map(|m| m.right_mil).collect();
    rights.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(rights, vec![-2.0, -1.5, -1.0, 0.0, 0.0, 0.0, 1.0, 1.5, 2.0]);

    let mut downs: Vec<f64> = description.marks.iter().map(|m| m.down_mil).collect();
    downs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(downs, vec![2.0, 2.0, 2.0, 4.0, 4.0, 4.0, 6.0, 6.0, 6.0]);

    // mil and mrad are the same unit spelled two ways; the row authored in mrad lands on
    // the same grid as the rows authored in mil.
    assert!(
        description
            .marks
            .iter()
            .any(|m| m.down_mil == 4.0 && m.right_mil == 1.5),
        "the mrad-authored row imports onto the mil grid"
    );

    // With marks on both axes the windage span is real, so a wind hold inside it is on the
    // reticle — the opposite of the elevation-only case below.
    let hold = hold_point_in_reticle(4.0, 1.5, 1.0, &description).unwrap();
    assert!(!hold.off_reticle);
    assert!(hold.nearest_mark_distance_mil < 1e-9);
}

/// The documented cost of the hold/decoration rule, asserted rather than described: the
/// duplex fixture etches four windage hashes, none of them is a `<bdc>` entry, so the
/// imported reticle offers nothing to hold on for wind.
#[test]
fn dropping_the_drawn_windage_hashes_makes_every_wind_hold_off_reticle() {
    let description = import_reticle_document(DUPLEX_BDC_MOA).unwrap();

    // A pure elevation hold lands on the middle holdover.
    let on = hold_point_in_reticle(5.0 * MOA_TO_MIL, 0.0, 1.0, &description).unwrap();
    assert!(on.nearest_mark.is_some());
    assert!(on.nearest_mark_distance_mil < 1e-9);
    assert!(!on.off_reticle);

    // The same hold with a breath of wind — well inside the ±6 MOA hashes the document
    // actually draws — reads as off-reticle, because those hashes are not imported.
    let windy =
        hold_point_in_reticle(5.0 * MOA_TO_MIL, 1.0 * MOA_TO_MIL, 1.0, &description).unwrap();
    assert!(
        windy.off_reticle,
        "an elevation-only mark set has a degenerate windage axis by construction"
    );
}

/// An all-decoration document is not an error here. The import is a pure transform; the
/// "there is nothing to aim with" judgement belongs to the solver.
#[test]
fn a_decoration_only_document_imports_with_no_marks() {
    let (description, report) = import_reticle_document_with_report(DECORATION_ONLY).unwrap();

    assert_eq!(description.name, "Fixture Decoration Only");
    assert!(description.marks.is_empty());
    assert_eq!(description.validate(), Err(ReticleError::NoMarks));
    assert!(matches!(
        hold_point_in_reticle(1.0, 0.0, 1.0, &description),
        Err(ReticleError::NoMarks)
    ));

    // Including the tag this importer has never heard of: unknown shapes are dropped like
    // any other drawing, not rejected.
    assert!(report
        .dropped_element_tags
        .iter()
        .any(|(tag, count)| tag == "reticle-illumination-dot" && *count == 1));
    assert_eq!(report.drawing_elements_dropped, 9);
}

/// Two documents the format allows but this importer cannot read: one structurally
/// malformed, one carrying a value with no unit. Both are refused rather than
/// half-imported.
#[test]
fn unreadable_documents_are_refused() {
    match import_reticle_document(MALFORMED_UNCLOSED) {
        Err(ReticleError::InvalidSpec(reason)) => assert!(
            reason.contains("elements"),
            "the error should name the element left open, got {reason:?}"
        ),
        other => panic!("an unclosed element must be an InvalidSpec, got {other:?}"),
    }

    match import_reticle_document(UNITLESS_VALUE) {
        Err(ReticleError::InvalidSpec(reason)) => assert!(
            reason.contains("position-y") && reason.contains("unit"),
            "the error should name the attribute and the missing unit, got {reason:?}"
        ),
        other => panic!("a unitless value must be an InvalidSpec, got {other:?}"),
    }
}

/// A hostile document cannot make the importer allocate without bound: the mark cap is
/// enforced while `<bdc>` entries are collected, not after.
#[test]
fn a_flood_of_bdc_entries_is_capped() {
    let mut xml = String::from(r#"<reticle name="Flood" size-x="10mil" size-y="10mil"><bdc>"#);
    for _ in 0..(MAX_RETICLE_MARKS * 4) {
        xml.push_str(r#"<bdc position-x="0mil" position-y="-1mil" />"#);
    }
    xml.push_str("</bdc></reticle>");

    assert!(matches!(
        import_reticle_document(&xml),
        Err(ReticleError::TooManyMarks { max, .. }) if max == MAX_RETICLE_MARKS
    ));
}
