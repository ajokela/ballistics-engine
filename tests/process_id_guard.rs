//! Public source and docs explain themselves: no comment, doc or release note cites the
//! humidity work's internal plan, steps, gates or review findings by an ID no reader of this
//! repository can look up (`plan C4`, `C7`, `H8`, `S5b`, `(R2)`, `gate G17`, `review HR-2`,
//! `humidity train`), or cites that unpublished plan at all ("the plan's grid", "the plan
//! allowed"). Each says the rule it stands for instead, or points to the CHANGELOG.
//!
//! Scanned, under the repository root (not `.git/`, a `target*` build directory, the generated
//! `doc/` or a wasm-pack `pkg/`; data such as JSON, CSV and fuzz seeds is not read): the comments
//! and prose string literals of every Rust and C source, read by a small lexer so a string
//! continued onto the next line, a raw string and a `"` in a char literal are each read as what
//! they are; every line of the Markdown docs, the scripts, the TOML and YAML files and the git
//! hooks; and the CHANGELOG's 0.50.0 section, which becomes the GitHub release notes. A phrase is
//! also looked for across a line wrap, so "the humidity" at one line's end and "train" at the
//! next line's start is found.
//!
//! Earlier work cites its own review IDs and uses the same spellings for other things
//! (`Tier 2 review C1`, `S8, 0.33.0 final-review fix wave`, a test's `S1`..`S8` scenario labels,
//! a truing or dial plan); those this vocabulary also matches are listed in [`EARLIER_WORK`], so
//! they stay visible here.

use std::fs;
use std::path::{Path, PathBuf};

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// A token that reads as one of the humidity work's IDs, by its form rather than a fixed list:
/// a step (`S` and digits, with an optional ASCII letter: `S1`, `S5b`, `S100B`),
/// a fix round (`R` and digits: `R2`, `R10`, `R700`), a plan item or an earlier review
/// round's finding (`C0a`, `C0b`, `C1`-`C19`), a finding (`H1`-`H24`, `T1`-`T9`), a gate
/// (`G9`-`G19`, `G8b`), or a review finding (`HR-`, `GAP-`, `DOC-`, `API-`, `PHYS-`, `UX-` or
/// `VAL-` and a number: `HR-4`, `UX-6`).
fn is_process_id(token: &str) -> bool {
    let number = |digits: &str, low: u32, high: u32| {
        !digits.is_empty()
            && digits.bytes().all(|b| b.is_ascii_digit())
            && digits
                .parse::<u32>()
                .is_ok_and(|n| (low..=high).contains(&n))
    };
    let digits = |text: &str| !text.is_empty() && text.bytes().all(|b| b.is_ascii_digit());
    if let Some((prefix, suffix)) = token.split_once('-') {
        return ["HR", "GAP", "DOC", "API", "PHYS", "UX", "VAL"].contains(&prefix)
            && digits(suffix);
    }
    let Some(first) = token.chars().next() else {
        return false;
    };
    let rest = &token[first.len_utf8()..];
    match first {
        'S' => {
            let suffix = rest
                .strip_suffix(|c: char| c.is_ascii_alphabetic())
                .unwrap_or(rest);
            digits(suffix)
        }
        'R' => digits(rest),
        'C' => matches!(rest, "0a" | "0b") || number(rest, 1, 19),
        'H' => number(rest, 1, 24),
        'T' => number(rest, 1, 9),
        'G' => match rest.strip_suffix('b') {
            Some(gate) => number(gate, 1, 19),
            None => number(rest, 9, 19),
        },
        _ => false,
    }
}

/// Phrases that cite the plan or the review rather than say what they mean.
const PHRASES: [&str; 7] = [
    "plan step",
    "humidity train",
    "Humidity train",
    "humidity-consistency train",
    "Humidity-consistency train",
    "gate G",
    "gate R",
];

/// The plan items and review findings cited in the forms `plan C4`, `plan's C7`, `plan item
/// C7`, `review C8` and `review T3`: a capital letter and a digit right after the lead. Each is
/// returned as cited (`C4`), so an [`EARLIER_WORK`] entry for that item excuses it.
fn review_items(text: &str) -> Vec<String> {
    let mut items = Vec::new();
    let mut cites = |lead: &str, letter: &dyn Fn(char) -> bool| {
        for (at, word) in text.match_indices(lead) {
            let item: String = text[at + word.len()..]
                .chars()
                .take_while(|c| c.is_ascii_alphanumeric())
                .collect();
            let mut chars = item.chars();
            if chars.next().is_some_and(letter) && chars.next().is_some_and(|c| c.is_ascii_digit())
            {
                items.push(item);
            }
        }
    };
    let any_capital = |c: char| c.is_ascii_uppercase();
    cites("plan ", &any_capital);
    cites("plan's ", &any_capital);
    cites("plan item ", &any_capital);
    cites("review ", &|c| matches!(c, 'C' | 'T'));
    items
}

/// Whether `text` cites a plan item or review finding in one of [`review_items`]'s forms.
fn cites_a_review_finding(text: &str) -> bool {
    !review_items(text).is_empty()
}

/// "the plan" or "this plan" as a word (not `plan-truing`, `planned` or `true.plan`): a
/// reference to the humidity work's unpublished plan, unless [`EARLIER_WORK`] says it names a
/// truing or dial plan.
fn cites_the_plan(text: &str) -> bool {
    let lower = text.to_lowercase();
    ["the plan", "this plan"].iter().any(|lead| {
        lower.match_indices(lead).any(|(at, word)| {
            let before = lower[..at].chars().next_back();
            let mut after = lower[at + word.len()..].chars();
            let next = after.next();
            let word_starts = before.is_none_or(|c| !(c.is_alphanumeric() || c == '_'));
            let word_ends = match next {
                None => true,
                Some('.') => after.next().is_none_or(|c| !c.is_alphanumeric()),
                Some(c) => !(c.is_alphanumeric() || c == '_' || c == '-'),
            };
            word_starts && word_ends
        })
    })
}

/// `(path, token or substring)` pairs that are not the humidity work's IDs: earlier work's own
/// labels and review IDs, and legitimate uses of the same spelling. Every entry must still be
/// needed ([`public_source_and_docs_cite_no_internal_plan_or_review_ids`] checks), so the list
/// cannot outlive what it excuses.
const EARLIER_WORK: [(&str, &str); 43] = [
    // Example zero-set and profile-row names, and a test CSV row's name: data, not IDs.
    ("CLI_USAGE.md", "R1"),
    ("CLI_USAGE.md", "R700"),
    ("src/main.rs", "R700"),
    // The initial target range in the moving-target intercept equation.
    ("docs/SOLVE_JSON_V1.md", "R0"),
    // A dispersion statistic's name: the radius holding 95 % of the shots.
    ("examples/02_monte_carlo.rs", "R95"),
    ("examples/README.md", "R95"),
    ("src/main.rs", "available zero sets"),
    ("tests/zero_sets.rs", "available zero sets"),
    ("tests/csv_numeric_cells.rs", "S1"),
    // Earlier work's review IDs (the 0.33.0 final review, the Tier 2 review, the MBA-1345,
    // MBA-1348 and MBA-1357 reviews, the JSON bridge's C2 item, the station sentinel's scenario
    // labels).
    ("src/adjustment.rs", "review C1"),
    ("src/explain.rs", "C1"),
    ("src/explain.rs", "S8"),
    ("src/main.rs", "C2"),
    ("src/main.rs", "MBA-1357 Task 2 review"),
    ("src/main.rs", "design R9"),
    ("src/main.rs", "the Plan B final review"),
    ("src/cli_api.rs", "C2"),
    ("src/optic.rs", "C1"),
    ("src/perturbation/taxonomy.rs", "review C1"),
    ("src/solve_v1.rs", "S5 (0.33.0"),
    ("src/wasm.rs", "(R6)"),
    ("src/wasm.rs", "Tier 2 review (C2)"),
    ("tests/adjustment_units_extended.rs", "C1"),
    ("tests/profile_bridge.rs", "(C2)"),
    ("tests/zero_range_from_angle.rs", "C2"),
    ("tests/station_presence_cli.rs", "S1"),
    ("tests/station_presence_cli.rs", "S1b"),
    ("tests/station_presence_cli.rs", "S2"),
    ("tests/station_presence_cli.rs", "S2b"),
    ("tests/station_presence_cli.rs", "S3"),
    ("tests/station_presence_cli.rs", "S4"),
    ("tests/station_presence_cli.rs", "S6"),
    ("tests/station_presence_cli.rs", "S8"),
    // A C1-continuous spline: the mathematical term.
    ("src/drag.rs", "remains C1 across"),
    // Earlier work's own plans (the Monte Carlo statistics, the DSF profile command, the fuzz
    // harness and the no_std spike).
    ("src/mc_stats.rs", "the plan"),
    ("tests/dsf_profile_cli.rs", "the plan: \"say so on stdout\""),
    ("scripts/fuzz-smoke.sh", "the plan's Global Constraints"),
    (
        "docs/embedded-no-std-spike-2026-07.md",
        "This is Phase 0 in the plan",
    ),
    // A truing plan, a dial or hold plan: the commands' own noun.
    ("src/truing_plan.rs", "the plan's drop unit"),
    ("src/optic.rs", "the plan's `limits_hit`"),
    ("src/main.rs", "limit violations this plan hit"),
    (
        "tests/station_presence_cli.rs",
        "the QNH must move the plan",
    ),
    ("CLI_USAGE.md", "the plan is still returned"),
];

/// This file: its own word list and examples.
const THIS_FILE: &str = "tests/process_id_guard.rs";

fn allowed_by(path: &str, line: &str, hit: &str) -> Option<usize> {
    EARLIER_WORK.iter().position(|&(file, what)| {
        file == path
            && (what == hit || (what.len() > 2 && !is_process_id(what) && line.contains(what)))
    })
}

thread_local! {
    static USED: std::cell::RefCell<Vec<bool>> = std::cell::RefCell::new(vec![false; EARLIER_WORK.len()]);
}

fn allowed(path: &str, line: &str, hit: &str) -> bool {
    if path == THIS_FILE {
        return true;
    }
    match allowed_by(path, line, hit) {
        Some(index) => {
            USED.with(|used| used.borrow_mut()[index] = true);
            true
        }
        None => false,
    }
}

/// The citations in `text`, one comment, string or document line. `carried` is the last word of
/// the line before it, when that line was prose too, so a phrase split by a line wrap is found.
fn citations_after(path: &str, carried: &str, text: &str) -> Vec<String> {
    let joined = if carried.is_empty() {
        text.to_string()
    } else {
        format!("{carried} {}", body(text))
    };
    let mut found = Vec::new();
    for phrase in PHRASES {
        if joined.contains(phrase) && !allowed(path, &joined, phrase) {
            found.push(phrase.to_string());
        }
    }
    for item in review_items(&joined) {
        if !allowed(path, &joined, &item) {
            found.push(format!("plan or review {item}"));
        }
    }
    if cites_the_plan(&joined) && !allowed(path, &joined, "the plan") {
        found.push("the plan".into());
    }
    for token in text.split(|c: char| !(c.is_ascii_alphanumeric() || c == '-')) {
        if is_process_id(token) && !allowed(path, &joined, token) {
            found.push(token.to_string());
        }
    }
    found
}

/// The citations in one line of prose on its own.
fn citations(path: &str, text: &str) -> Vec<String> {
    citations_after(path, "", text)
}

/// A comment or document line without its leading markers (`//`, `///`, `//!`, `/*`, `*`, `#`)
/// and indentation.
fn body(text: &str) -> &str {
    text.trim_start_matches(|c: char| c.is_whitespace() || matches!(c, '/' | '!' | '*' | '#'))
}

/// The last word of a line, for [`citations_after`].
fn last_word(text: &str) -> &str {
    text.split_whitespace().next_back().unwrap_or("")
}

/// One line's worth of prose: its 1-based line number and its text.
type Prose = (usize, String);

/// The prose of a Rust or C source, line by line: each comment, and each string literal that
/// reads as prose (it has a space somewhere, on any of its lines), cut at its line breaks.
/// A small lexer rather than a per-line split, so a string continued onto the next line, a `//`
/// inside a string, a `"` inside a char literal and a raw string are each read as what they are.
fn source_prose(src: &str) -> Vec<Prose> {
    let chars: Vec<char> = src.chars().collect();
    let at = |i: usize| chars.get(i).copied();
    let is_ident = |c: char| c.is_alphanumeric() || c == '_';
    let mut out: Vec<Prose> = Vec::new();
    let mut line = 1;
    let mut i = 0;
    while let Some(c) = at(i) {
        match c {
            '\n' => {
                line += 1;
                i += 1;
            }
            '/' if at(i + 1) == Some('/') => {
                let start = i;
                while at(i).is_some_and(|c| c != '\n') {
                    i += 1;
                }
                out.push((line, chars[start..i].iter().collect()));
            }
            '/' if at(i + 1) == Some('*') => {
                let mut depth = 0;
                let mut text = String::new();
                loop {
                    match (at(i), at(i + 1)) {
                        (None, _) => break,
                        (Some('/'), Some('*')) => {
                            depth += 1;
                            text.push_str("/*");
                            i += 2;
                        }
                        (Some('*'), Some('/')) => {
                            depth -= 1;
                            text.push_str("*/");
                            i += 2;
                            if depth == 0 {
                                break;
                            }
                        }
                        (Some('\n'), _) => {
                            out.push((line, std::mem::take(&mut text)));
                            line += 1;
                            i += 1;
                        }
                        (Some(other), _) => {
                            text.push(other);
                            i += 1;
                        }
                    }
                }
                out.push((line, text));
            }
            '"' => {
                let mut segments = vec![(line, String::new())];
                i += 1;
                loop {
                    match at(i) {
                        None => break,
                        Some('"') => {
                            i += 1;
                            break;
                        }
                        // A continuation: the next line's indentation is not in the string.
                        Some('\\') if at(i + 1) == Some('\n') => {
                            i += 1;
                            while let Some(space) = at(i).filter(|c| c.is_whitespace()) {
                                if space == '\n' {
                                    line += 1;
                                    segments.push((line, String::new()));
                                }
                                i += 1;
                            }
                        }
                        Some('\\') => {
                            let segment = &mut segments.last_mut().unwrap().1;
                            segment.push('\\');
                            segment.extend(at(i + 1));
                            i += 2;
                        }
                        Some('\n') => {
                            line += 1;
                            segments.push((line, String::new()));
                            i += 1;
                        }
                        Some(other) => {
                            segments.last_mut().unwrap().1.push(other);
                            i += 1;
                        }
                    }
                }
                push_literal(&mut out, segments);
            }
            // A raw string, `r"…"`, `r#"…"#` or `br"…"`; anything else that starts with `r`
            // or `b` (an identifier, `r#type`, a byte string) falls through to the next
            // character.
            'r' | 'b' if i == 0 || !is_ident(chars[i - 1]) => {
                let mut j = if c == 'b' && at(i + 1) == Some('r') {
                    i + 2
                } else if c == 'r' {
                    i + 1
                } else {
                    i += 1;
                    continue;
                };
                let mut hashes = 0;
                while at(j) == Some('#') {
                    hashes += 1;
                    j += 1;
                }
                if at(j) != Some('"') {
                    i += 1;
                    continue;
                }
                let mut segments = vec![(line, String::new())];
                i = j + 1;
                loop {
                    match at(i) {
                        None => break,
                        Some('"') if (1..=hashes).all(|h| at(i + h) == Some('#')) => {
                            i += 1 + hashes;
                            break;
                        }
                        Some('\n') => {
                            line += 1;
                            segments.push((line, String::new()));
                            i += 1;
                        }
                        Some(other) => {
                            segments.last_mut().unwrap().1.push(other);
                            i += 1;
                        }
                    }
                }
                push_literal(&mut out, segments);
            }
            // A char literal (`'x'`, `'"'`, `'\''`, `'\u{b0}'`) is skipped whole, so its quote
            // opens no string; a lifetime or label (`'a`, `'static`) is one character.
            '\'' => {
                if at(i + 1) == Some('\\') {
                    let mut k = i + 3;
                    while at(k).is_some_and(|c| c != '\'' && c != '\n') {
                        k += 1;
                    }
                    i = if at(k) == Some('\'') { k + 1 } else { k };
                } else if at(i + 1).is_some_and(|c| c != '\n') && at(i + 2) == Some('\'') {
                    i += 3;
                } else {
                    i += 1;
                }
            }
            _ => i += 1,
        }
    }
    out
}

/// A string literal's text, one segment per source line it covers, kept when it reads as prose:
/// a value such as `"S1"` or a CSV row has no space.
fn push_literal(out: &mut Vec<Prose>, segments: Vec<Prose>) {
    if segments.iter().any(|(_, text)| text.contains(' ')) {
        out.extend(
            segments
                .into_iter()
                .filter(|(_, text)| !text.trim().is_empty()),
        );
    }
}

/// The prose of a file, line by line: a source's comments and prose strings, or every line of a
/// document, script or configuration file.
fn file_prose(relative: &str, text: &str) -> Vec<Prose> {
    if [".rs", ".c", ".h"]
        .iter()
        .any(|ext| relative.ends_with(ext))
    {
        source_prose(text)
    } else {
        text.lines()
            .enumerate()
            .map(|(index, line)| (index + 1, line.to_string()))
            .collect()
    }
}

/// The citations in a run of prose lines, each found with the last word of the line before
/// carried over when the two lines are adjacent.
fn citations_in(path: &str, prose: &[Prose]) -> Vec<(usize, String, String)> {
    let mut found = Vec::new();
    let mut previous: Option<&Prose> = None;
    for item in prose {
        let (number, text) = item;
        let carried = match previous {
            Some((before, before_text)) if before + 1 == *number => last_word(before_text),
            _ => "",
        };
        for hit in citations_after(path, carried, text) {
            found.push((*number, hit, text.trim().to_string()));
        }
        previous = Some(item);
    }
    found
}

const SCANNED: [&str; 10] = [
    ".rs", ".c", ".h", ".md", ".sh", ".py", ".mjs", ".toml", ".yml", ".yaml",
];

fn walk(dir: &Path, top: bool, out: &mut Vec<PathBuf>) {
    let Ok(entries) = fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        let name = entry.file_name();
        let name = name.to_string_lossy();
        if path.is_dir() {
            let skipped = matches!(&*name, ".git" | "node_modules" | "pkg")
                || name.starts_with("target")
                || (top && name == "doc");
            if !skipped {
                walk(&path, false, out);
            }
        } else if (SCANNED.iter().any(|ext| name.ends_with(ext))
            && !(top && name == "CHANGELOG.md"))
            || in_git_hooks(&path)
        {
            out.push(path);
        }
    }
}

/// Whether `path` is one of the versioned git hooks: every file in `.githooks/`, whatever its
/// name (`pre-commit`, `commit-msg`, and any hook added later).
fn in_git_hooks(path: &Path) -> bool {
    path.parent()
        .and_then(Path::file_name)
        .is_some_and(|directory| directory == ".githooks")
}

/// The CHANGELOG's 0.50.0 section: under `## [0.50.0]` once cut, `## [Unreleased]` before.
fn changelog_section(changelog: &str) -> &str {
    let start = changelog
        .find("\n## [0.50.0]")
        .or_else(|| changelog.find("\n## [Unreleased]"))
        .expect("CHANGELOG.md has the 0.50.0 section");
    let body = &changelog[start + 1..];
    let end = body[3..].find("\n## [").map_or(body.len(), |at| at + 3);
    &body[..end]
}

#[test]
fn public_source_and_docs_cite_no_internal_plan_or_review_ids() {
    let root = root();
    let mut files = Vec::new();
    walk(&root, true, &mut files);
    let relative = |file: &Path| {
        file.strip_prefix(&root)
            .unwrap()
            .to_string_lossy()
            .replace('\\', "/")
    };
    let names: Vec<String> = files.iter().map(|file| relative(file)).collect();
    assert!(
        names.len() > 250,
        "the walk found only {} files",
        names.len()
    );
    for expected in [
        "src/lib.rs",
        "src/atmosphere/aloft_vapour_tests.rs",
        "tests/humidity_byte_identity.rs",
        "benches/humidity_air_bench.rs",
        "fuzz/src/domain.rs",
        "include/ballistics_bridge.h",
        "scripts/build-npm-postprocess.mjs",
        "README.md",
        "README-npm.md",
        "CLI_USAGE.md",
        "MCCOY_MIGRATION.md",
        "Cargo.toml",
        "fuzz/Cargo.toml",
        ".github/workflows/ci.yml",
        ".githooks/pre-commit",
        ".githooks/commit-msg",
        "scripts/fuzz-smoke.sh",
        "tests/fixtures/analyze_courtney_cd.py",
    ] {
        // A source export may leave the dotfiles out; a file that is here must be read.
        assert!(
            !root.join(expected).exists() || names.iter().any(|name| name == expected),
            "the walk missed {expected}"
        );
    }
    assert!(
        !names.iter().any(|name| name == "CHANGELOG.md"
            || name.starts_with("target")
            || name.starts_with("doc/")),
        "the walk read the whole CHANGELOG, a build directory or the generated doc/"
    );

    let mut offenders = Vec::new();
    for (file, name) in files.iter().zip(&names) {
        let Ok(text) = fs::read_to_string(file) else {
            continue;
        };
        for (number, hit, text) in citations_in(name, &file_prose(name, &text)) {
            offenders.push(format!("{name}:{number}: `{hit}` in {text}"));
        }
    }
    let changelog = fs::read_to_string(root.join("CHANGELOG.md")).unwrap();
    let section = file_prose("CHANGELOG.md", changelog_section(&changelog));
    for (number, hit, text) in citations_in("CHANGELOG.md (0.50.0 section)", &section) {
        offenders.push(format!(
            "CHANGELOG.md (0.50.0 section):{number}: `{hit}` in {text}"
        ));
    }
    assert!(
        offenders.is_empty(),
        "{} lines cite an internal plan, step, gate or review ID; say the rule instead:\n{}",
        offenders.len(),
        offenders.join("\n")
    );
    let unused: Vec<_> = USED.with(|used| {
        EARLIER_WORK
            .iter()
            .zip(used.borrow().iter())
            .filter(|(_, &used)| !used)
            .map(|(entry, _)| *entry)
            .collect()
    });
    assert!(
        unused.is_empty(),
        "EARLIER_WORK entries no line needs any more: {unused:?}"
    );
}

/// The detector itself: every ID form it exists to catch, and the spellings it must leave alone.
#[test]
fn the_detector_knows_an_id_from_a_name() {
    // A step or fix round of any number, not a fixed list: the work outgrew `R1`-`R7` and
    // `S1`-`S8` once, and `R8`, `R10` and `S9` went unseen.
    for id in [
        "S1",
        "S5b",
        "S8",
        "S9",
        "S12",
        "S12a",
        "R2",
        "R7",
        "R8",
        "R10",
        "R42",
        "C0a",
        "C0b",
        "C4",
        "C7",
        "C11",
        "H5",
        "H24",
        "T3",
        "G17",
        "G8b",
        "HR-4",
        "HR-12",
        "GAP-3",
        "UX-6",
        "VAL-1",
        "PHYS-5",
        "DOC-17",
        "API-1000",
        "R0",
        "S0",
        "HR-0",
        "R700",
        "S100",
        "S100B",
        "R999999999999999999999",
        "S999999999999999999999b",
        "API-999999999999999999999",
    ] {
        assert!(is_process_id(id), "{id} is an ID");
    }
    for name in [
        "G1", "G7", "G8", "R1a", "S5bb", "S5BB", "C0", "C20", "C99", "H25", "T0", "MBA-1605",
        "H2O", "S", "R", "C", "HR-", "API", "x_v", "1e-15",
    ] {
        assert!(!is_process_id(name), "{name} is not an ID");
    }
    assert!(cites_a_review_finding("see review C8 for the rule"));
    assert!(cites_a_review_finding("the one validator (plan C4)"));
    assert!(cites_a_review_finding("// The plan's C7 grid: T, P"));
    assert!(cites_a_review_finding("the plan item H5 says"));
    assert!(!cites_a_review_finding(
        "the stable dial-plan CLI command shape"
    ));
    assert!(cites_a_review_finding("(whole-branch review T3)"));
    assert!(!cites_a_review_finding("a review of the CHANGELOG"));
    for cite in [
        "spanning the plan's -40..45 C grid",
        "The two cases the plan names",
        "the plan allowed last-ulp movement",
        "as the plan said.",
        "(see the plan)",
        "THE PLAN: one validator",
        "this plan's order",
    ] {
        assert!(cites_the_plan(cite), "{cite:?} cites the plan");
    }
    for name in [
        "`plan-truing` reported a shortage",
        "`true.plan` reports any failure",
        "the planned shot",
        "the planner",
        "a request it cannot plan.",
        "the plan.json file",
    ] {
        assert!(!cites_the_plan(name), "{name:?} does not cite the plan");
    }
    assert_eq!(
        citations("src/x.rs", "// the station air (plan C4): one validator"),
        ["plan or review C4", "C4"]
    );
    assert_eq!(review_items("review C11 and plan's C7"), ["C7", "C11"]);
    assert!(citations("src/x.rs", "// the solver's G7 drag, at 15 °C").is_empty());
    // A legitimate row name or statistic excuses only its own token, not another ID
    // that someone appends to the same line.
    assert_eq!(
        citations("CLI_USAGE.md", "--profile-row R700_308 // R10"),
        ["R10"]
    );
    assert_eq!(
        citations("examples/README.md", "R95 dispersion, measured in S100B"),
        ["S100B"]
    );
}

/// The lexer reads each kind of Rust and C prose, and each look-alike, as what it is.
#[test]
fn the_lexer_reads_comments_and_prose_strings_as_a_reader_does() {
    let prose = |src: &str| source_prose(src);
    let lines = |src: &str| -> Vec<(usize, String)> {
        source_prose(src)
            .into_iter()
            .map(|(number, text)| (number, text.trim().to_string()))
            .collect()
    };
    // A prose string on its own line, and a string that is a value, not prose.
    assert_eq!(
        lines("    assert!(ok, \"S8 may change it\");"),
        [(1, "S8 may change it".into())]
    );
    assert!(prose("    let site = \"S1\";").is_empty());
    // A string continued onto a second line: the continuation is prose too, on its own line.
    let continued = "            assert_eq!(\n                a,\n                \"{solver_name}: must match its pin (the pre-change flat-Vec \\\n                 capture, re-captured in 0.50.0's S8)\"\n            );\n";
    assert_eq!(
        lines(continued),
        [
            (
                3,
                "{solver_name}: must match its pin (the pre-change flat-Vec".into()
            ),
            (4, "capture, re-captured in 0.50.0's S8)".into()),
        ]
    );
    assert_eq!(
        citations_in("src/x.rs", &prose(continued))
            .into_iter()
            .map(|(number, hit, _)| (number, hit))
            .collect::<Vec<_>>(),
        [(4, "S8".to_string())]
    );
    // A CSV or other multi-line value with no space is data, not prose.
    assert!(prose("let csv = \"A,B\\n\\\n         R1,2\\n\\\n         R2,3\\n\";").is_empty());
    // A string literal that spans lines with no continuation backslash.
    assert_eq!(
        lines("let s = \"one\nS5b two\";"),
        [(1, "one".into()), (2, "S5b two".into())]
    );
    // A `//` inside a string opens no comment; a `"` in a char literal opens no string.
    assert_eq!(
        lines("let url = \"https://example.org/a b\"; let q = '\"'; // S8 here"),
        [
            (1, "https://example.org/a b".into()),
            (1, "// S8 here".into())
        ]
    );
    assert_eq!(
        lines("let e = '\\''; let d = '\\u{b0}'; // after (R2)"),
        [(1, "// after (R2)".into())]
    );
    // Lifetimes and labels are not char literals.
    assert_eq!(
        lines("fn f<'a>(x: &'a str) -> &'a str { 'outer: loop { break \"a b S8\"; } }"),
        [(1, "a b S8".into())]
    );
    // Raw strings, with and without hashes, and a raw identifier that is not one.
    assert_eq!(
        lines("let r = r#\"the \"S8\" step\"#; let s = br\"x y\"; let r#type = 1; // end"),
        [
            (1, "the \"S8\" step".into()),
            (1, "x y".into()),
            (1, "// end".into())
        ]
    );
    // A block comment over lines, nested as Rust allows, and a C comment's `*` lines.
    assert_eq!(
        lines("/* one /* two */\n * S8 three */ code"),
        [(1, "/* one /* two */".into()), (2, "* S8 three */".into())]
    );
    // A phrase split by a line wrap is found on the line that completes it.
    let wrapped = "    // the rule of the humidity\n    // train, and\n    /// as the\n    /// plan allowed\n";
    assert_eq!(
        citations_in("src/x.rs", &prose(wrapped))
            .into_iter()
            .map(|(number, hit, _)| (number, hit))
            .collect::<Vec<_>>(),
        [
            (2, "humidity train".to_string()),
            (4, "the plan".to_string())
        ]
    );
    // Lines that are not adjacent carry nothing over.
    let apart = "    // the humidity\n    let x = 1;\n    // train";
    assert!(citations_in("src/x.rs", &prose(apart)).is_empty());
}
