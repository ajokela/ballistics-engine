//! Every test that runs the `ballistics` binary runs the one cargo built for this test run.
//!
//! Integration tests get that path from `env!("CARGO_BIN_EXE_ballistics")`. Sixteen of them used
//! to build it by hand instead, as `CARGO_MANIFEST_DIR/target/debug/ballistics`, some with a
//! `target/release` fallback, and `tests/test_monte_carlo.py` ran `./target/debug/ballistics`.
//! A hand-built path has two faults:
//!
//! - it ignores `CARGO_TARGET_DIR`, so a run with its own target directory fails with `NotFound`;
//! - it does not follow the build under test, so under `cargo test --release`, or in a checkout
//!   holding an old `./target/debug/ballistics`, the test runs that old binary and passes on it.
//!
//! This scans the text of every test file under `tests/` (`.rs`, `.py`, `.sh`) and fails if a
//! binary path is built under `target/` again. Whole-line comments are skipped, so a comment may
//! still explain the old pattern (as `tests/drag_curve_cli.rs` does).

use std::fs;
use std::path::{Path, PathBuf};

/// Recursively collect the test sources under `dir`.
fn collect_test_sources(dir: &Path, out: &mut Vec<PathBuf>) {
    let entries = fs::read_dir(dir).unwrap_or_else(|error| panic!("{}: {error}", dir.display()));
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            collect_test_sources(&path, out);
        } else if matches!(
            path.extension().and_then(|extension| extension.to_str()),
            Some("rs" | "py" | "sh")
        ) {
            out.push(path);
        }
    }
}

/// Whether a source line builds a path under the cargo target directory by hand. `comment` is
/// the file's line-comment marker; a line that is only a comment is never flagged.
fn builds_a_target_path(line: &str, comment: &str) -> bool {
    let code = line.trim_start();
    if code.starts_with(comment) {
        return false;
    }
    [
        "target/debug",
        "target/release",
        "push(\"target\")",
        "join(\"target",
        "\"./target",
        "'./target",
    ]
    .iter()
    .any(|pattern| code.contains(pattern))
}

#[test]
fn the_detector_flags_every_hand_built_binary_path_and_nothing_else() {
    // The shapes the sixteen helpers and the Python script used.
    for line in [
        "    PathBuf::from(env!(\"CARGO_MANIFEST_DIR\")).join(\"target/debug/ballistics\")",
        "    path.push(\"target\");",
        "        p.push(\"release\"); // target/release fallback",
        "    \"./target/debug/ballistics\", \"monte-carlo\",",
        "    let bin = root.join(\"target\").join(\"debug\");",
    ] {
        assert!(builds_a_target_path(line, "//"), "not flagged: {line}");
    }
    // `target` as a value, and a comment about the old pattern, are not binary paths.
    for line in [
        "    table_args.extend([\"--shooting-angle\", \"30\", \"--drops-reference\", \"target\"]);",
        "// target/debug with a release fallback, which can byte-compare against a STALE binary",
        "const BIN: &str = env!(\"CARGO_BIN_EXE_ballistics\");",
    ] {
        assert!(!builds_a_target_path(line, "//"), "flagged: {line}");
    }
    assert!(!builds_a_target_path(
        "# ./target/debug/ballistics used to be hard-coded here",
        "#"
    ));
}

#[test]
fn no_test_builds_the_binary_path_under_target() {
    let tests_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let mut files = Vec::new();
    collect_test_sources(&tests_dir, &mut files);
    files.sort();
    assert!(
        files.len() > 100,
        "sanity check failed: found only {} test sources under {}",
        files.len(),
        tests_dir.display()
    );

    let this_file = Path::new(file!())
        .file_name()
        .expect("this test's file name");
    let mut offenders = Vec::new();
    for path in &files {
        // The detector's own patterns live here.
        if path.file_name() == Some(this_file) {
            continue;
        }
        let comment = if path.extension().and_then(|extension| extension.to_str()) == Some("rs") {
            "//"
        } else {
            "#"
        };
        let text = fs::read_to_string(path)
            .unwrap_or_else(|error| panic!("read {}: {error}", path.display()));
        for (index, line) in text.lines().enumerate() {
            if builds_a_target_path(line, comment) {
                let relative = path
                    .strip_prefix(env!("CARGO_MANIFEST_DIR"))
                    .unwrap_or(path);
                offenders.push(format!(
                    "{}:{}: {}",
                    relative.display(),
                    index + 1,
                    line.trim()
                ));
            }
        }
    }
    assert!(
        offenders.is_empty(),
        "these tests build the binary's path under target/ by hand; use \
         env!(\"CARGO_BIN_EXE_ballistics\") (or `cargo run --bin ballistics` from a script):\n{}",
        offenders.join("\n")
    );
}
