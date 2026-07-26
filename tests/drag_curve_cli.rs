//! MBA-1424 / MBA-1426: the `drag-curve` command, and the shared-formatter contract behind it.
//!
//! `drag::format_reference_drag_curve` is the single formatter both the native CLI and the
//! browser terminal emit — that sharing is what prevents the two surfaces drifting the way the
//! `recoil` CSV header did (MBA-1418). The WASM half of the proof is browser-gated and therefore
//! compile-verified only, so this file carries the executable half: the native binary's stdout
//! must equal the shared formatter's string byte for byte, for every model and every form.
//!
//! When `run_drag_curve` was rewired onto the formatter, all 27 outputs (9 models x 3 formats)
//! were byte-diffed against the pre-refactor binary and matched; these tests keep that true.

use std::process::Command;

use ballistics_engine::drag::{format_reference_drag_curve, ReferenceDragCurveFormat};
use ballistics_engine::DragModel;

// The cargo-provided path to the binary built for THIS test run. The first draft hand-rolled
// target/debug with a release fallback, which can byte-compare against a STALE binary — under
// `cargo test --release`, a leftover debug build from any earlier commit would be the one
// tested, and a regression in the freshly built binary would pass. Caught in review; the newer
// sibling tests (dsf_profile_cli, compare_cli) already use this form.
const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

const ALL_MODELS: [(&str, DragModel); 9] = [
    ("g1", DragModel::G1),
    ("g2", DragModel::G2),
    ("g5", DragModel::G5),
    ("g6", DragModel::G6),
    ("g7", DragModel::G7),
    ("g8", DragModel::G8),
    ("gi", DragModel::GI),
    ("gs", DragModel::GS),
    ("ra4", DragModel::RA4),
];

/// The binary's stdout is the shared formatter's string, verbatim, for all 27 combinations.
/// Equality — not fragment spot-checks — because a fragment check is exactly what let the
/// recoil header divergence survive (the checked fragment happened to agree).
#[test]
fn native_output_is_the_shared_formatter_verbatim_for_every_model_and_form() {
    for (name, model) in ALL_MODELS {
        for (flag, format) in [
            ("table", ReferenceDragCurveFormat::Table),
            ("csv", ReferenceDragCurveFormat::Csv),
            ("json", ReferenceDragCurveFormat::Json),
        ] {
            let output = Command::new(BIN)
                .args(["drag-curve", "--drag-model", name, "-o", flag])
                .output()
                .expect("drag-curve");
            assert!(
                output.status.success(),
                "drag-curve {name} {flag} failed: {}",
                String::from_utf8_lossy(&output.stderr)
            );
            assert_eq!(
                String::from_utf8_lossy(&output.stdout),
                format_reference_drag_curve(&model, format),
                "native `drag-curve --drag-model {name} -o {flag}` diverged from the shared \
                 formatter"
            );
        }
    }
}

/// The JSON form parses, self-describes its per-table Mach domain, and its point count matches
/// the table the solver interpolates.
#[test]
fn json_form_is_parseable_and_self_describing()
{
    let output = Command::new(BIN)
        .args(["drag-curve", "--drag-model", "gs", "-o", "json"])
        .output()
        .expect("drag-curve");
    let doc: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("drag-curve json parses");

    let table = ballistics_engine::drag::reference_drag_table(&DragModel::GS);
    assert_eq!(doc["drag_model"], "GS");
    assert_eq!(doc["point_count"].as_u64().unwrap() as usize, table.mach_values.len());
    // GS genuinely stops at Mach 4 — the reason the domain is stated per-table.
    assert_eq!(doc["mach_max"].as_f64().unwrap(), 4.0);
    assert_eq!(
        doc["points"].as_array().unwrap().len(),
        table.mach_values.len()
    );
}

/// PDF stays refused with the documented message.
#[test]
fn pdf_is_refused_with_the_documented_message() {
    let output = Command::new(BIN)
        .args(["drag-curve", "-o", "pdf"])
        .output()
        .expect("drag-curve");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("drag-curve has no PDF form"),
        "unexpected refusal text: {stderr}"
    );
}
