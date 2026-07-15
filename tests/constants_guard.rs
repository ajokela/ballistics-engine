//! MBA-1327: guards against a duplicate, truncated, or rounded grain<->gram
//! conversion constant creeping back into `src/` outside `src/constants.rs`.
//!
//! Before this fix, three different approximations of the exact grain<->gram
//! conversion (0.45359237 kg / 7000 grains) coexisted in the crate:
//!   - `0.0647989`  (truncated grams-per-grain, in main.rs and api_client.rs)
//!   - `15.4323584` (rounded grains-per-gram, in wasm.rs)
//!   - `0.06479891` (the correct, exact value, duplicated locally in main.rs
//!     instead of coming from a shared constant)
//!
//! `src/constants.rs` now defines the single source of truth
//! (`GRAMS_PER_GRAIN` / `GRAINS_PER_GRAM`). This test scans the raw text of
//! every other file under `src/` (not just parsed code -- comments count
//! too, since a stray comment reintroducing the wrong digits is just as
//! misleading as code doing it) and fails if the truncated or rounded
//! literals reappear anywhere else.
//!
//! Deliberately scoped to `src/` only (not `tests/`): test fixtures are
//! allowed to hardcode the exact value independently as a cross-check
//! against the production constant.

use std::fs;
use std::path::{Path, PathBuf};

/// Recursively collect `.rs` file paths under `dir`.
fn collect_rs_files(dir: &Path, out: &mut Vec<PathBuf>) {
    let entries = match fs::read_dir(dir) {
        Ok(e) => e,
        Err(_) => return,
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            collect_rs_files(&path, out);
        } else if path.extension().and_then(|e| e.to_str()) == Some("rs") {
            out.push(path);
        }
    }
}

/// True if `text` contains the truncated grams-per-grain literal
/// `0.0647989` that is *not* immediately followed by a further `1`.
///
/// This deliberately does not flag the correct, exact `0.06479891` (which
/// is `0.0647989` followed by `1`) -- only the previously-shipped truncated
/// 7-significant-digit approximation.
fn contains_truncated_grams_per_grain(text: &str) -> bool {
    let needle = "0.0647989";
    text.match_indices(needle)
        .any(|(idx, _)| !text[idx + needle.len()..].starts_with('1'))
}

/// True if `text` contains the previously-shipped rounded grains-per-gram
/// literal (`15.4323584`) or any of its truncated prefixes.
fn contains_rounded_grains_per_gram(text: &str) -> bool {
    text.contains("15.4323") || text.contains("15.432358")
}

#[test]
fn no_duplicate_grain_gram_constants_outside_constants_rs() {
    let src_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
    let mut files = Vec::new();
    collect_rs_files(&src_dir, &mut files);
    assert!(
        !files.is_empty(),
        "sanity check failed: found no .rs files under {}",
        src_dir.display()
    );

    let mut offenders = Vec::new();
    for path in files {
        if path.file_name().and_then(|n| n.to_str()) == Some("constants.rs") {
            continue; // the single source of truth is allowed to define these
        }
        let Ok(text) = fs::read_to_string(&path) else {
            continue; // non-UTF8/unreadable file is not this test's concern
        };
        if contains_truncated_grams_per_grain(&text) || contains_rounded_grains_per_gram(&text) {
            offenders.push(path.display().to_string());
        }
    }

    assert!(
        offenders.is_empty(),
        "MBA-1327: found a duplicate/truncated/rounded grain<->gram conversion \
         literal outside src/constants.rs -- use \
         ballistics_engine::constants::{{GRAMS_PER_GRAIN, GRAINS_PER_GRAM}} instead:\n{}",
        offenders.join("\n")
    );
}
