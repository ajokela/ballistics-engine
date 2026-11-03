//! The browser terminal's help lists every flag its parser takes.
//!
//! `src/wasm.rs` parses each command's flags by hand and writes its help by hand, and the two
//! drifted: `trajectory` parsed `--zero-humidity`, `--zero-temperature`, `--zero-pressure`,
//! `--zero-pressure-type`, `--zero-altitude`, `--pressure-type` and `--density-altitude` (and
//! more), which CLI_USAGE.md documents, but its help listed none of them as an entry, so a
//! browser user could not find the zero day's air. The terminal only compiles for wasm32, so
//! this host test reads its source: for every `handle_<command>_command`, each `"--flag" =>`
//! arm must appear as an entry (a help line that starts with the flag) in that command's help
//! section, or be one of the listed exceptions below.

use std::collections::BTreeSet;

const SOURCE: &str = include_str!("../src/wasm.rs");

/// The deepest indent an entry line has: options sit at 2 to 4 spaces, descriptions at 33.
const ENTRY_INDENT_MAX: usize = 8;

/// `(command, flag, why it needs no entry of its own)`.
const EXCEPTIONS: &[(&str, &str, &str)] = &[
    (
        "*",
        "--units",
        "a global option, listed once under Global Options",
    ),
    (
        "trajectory",
        "--cant-angle",
        "an alias the --cant entry names",
    ),
    (
        "trajectory",
        "--bore-height",
        "an alias the --muzzle-height entry names",
    ),
    (
        "monte_carlo",
        "--wind-std",
        "an alias the --wind-speed-std entry names",
    ),
    (
        "monte_carlo",
        "--wind-dir-std",
        "an alias the --wind-direction-std entry names",
    ),
    (
        "reticle",
        "--profile",
        "native-only; the help says so and the arm refuses it",
    ),
    (
        "reticle",
        "--range",
        "native-only; the help says so and the arm refuses it",
    ),
];

/// The help title of each handler.
const TITLES: &[(&str, &str)] = &[
    ("trajectory", "Trajectory"),
    ("zero", "Zero"),
    ("lead", "Lead"),
    ("monte_carlo", "Monte Carlo"),
    ("true_velocity", "True Velocity"),
    ("true_wind", "True Wind"),
    ("estimate_bc", "Estimate BC"),
    ("bc_convert", "BC Convert"),
    ("drag_curve", "Drag Curve"),
    ("reticle", "Reticle"),
    ("powder", "Powder"),
    ("recoil", "Recoil"),
    ("power_factor", "Power Factor"),
];

/// The flags of every `"--flag" =>` (or `"--a" | "--b" =>`) match arm in `body`.
fn parsed_flags(body: &str) -> BTreeSet<String> {
    let mut flags = BTreeSet::new();
    for line in body.lines() {
        let trimmed = line.trim_start();
        if !trimmed.starts_with("\"-") || !trimmed.contains("=>") {
            continue;
        }
        let pattern = trimmed.split("=>").next().unwrap_or_default();
        for alternative in pattern.split('|') {
            let alternative = alternative.trim().trim_matches('"');
            if alternative.starts_with("--") && !alternative.contains(' ') {
                flags.insert(alternative.to_string());
            }
        }
    }
    flags
}

/// `(command, body)` for every `fn handle_<command>_command`.
fn handlers() -> Vec<(String, String)> {
    let starts: Vec<(usize, String)> = SOURCE
        .match_indices("    fn handle_")
        .filter_map(|(at, _)| {
            let rest = &SOURCE[at + "    fn handle_".len()..];
            let name = rest.split('(').next()?;
            name.strip_suffix("_command")
                .map(|command| (at, command.to_string()))
        })
        .collect();
    starts
        .iter()
        .enumerate()
        .map(|(i, (at, command))| {
            let end = starts.get(i + 1).map_or(SOURCE.len(), |(next, _)| *next);
            (command.clone(), SOURCE[*at..end].to_string())
        })
        .collect()
}

/// The help section titled `"<title> Command:"`, up to the next title or the end of its text.
fn help_section(title: &str) -> String {
    let header = format!("\n{title} Command:\n");
    let start = SOURCE
        .find(&header)
        .unwrap_or_else(|| panic!("no help section {header:?}"))
        + header.len();
    let mut section = String::new();
    for line in SOURCE[start..].lines() {
        let is_title = line.ends_with(" Command:")
            && line.chars().next().is_some_and(|c| c.is_ascii_uppercase());
        if is_title {
            break;
        }
        if let Some((text, _)) = line.split_once("\"#") {
            section.push_str(text);
            break;
        }
        section.push_str(line);
        section.push('\n');
    }
    section
}

/// The flags a help section lists as entries: a line that starts, in the entry column, with
/// the flag, or with a short option and then the flag (`-o, --output <FORMAT>`). A flag at the
/// start of a description's continuation line (indented to the description column) is no
/// entry: that is how `--zero-pressure-type` hid inside `--zero-density-altitude`'s text.
fn help_entries(section: &str) -> BTreeSet<String> {
    let mut entries = BTreeSet::new();
    for line in section.lines() {
        let mut rest = line.trim_start();
        if line.len() - rest.len() > ENTRY_INDENT_MAX {
            continue;
        }
        if rest.len() > 3 && rest.starts_with('-') && rest.as_bytes()[1] != b'-' {
            match rest.split_once(", ") {
                Some((_, after)) => rest = after,
                None => continue,
            }
        }
        if let Some(flag) = rest
            .strip_prefix("--")
            .map(|tail| tail.split([' ', '<', '[', ',']).next().unwrap_or_default())
        {
            if !flag.is_empty() {
                entries.insert(format!("--{flag}"));
            }
        }
    }
    entries
}

#[test]
fn every_terminal_flag_has_a_help_entry() {
    let mut missing = Vec::new();
    let handlers = handlers();
    for &(command, title) in TITLES {
        let (_, body) = handlers
            .iter()
            .find(|(name, _)| name == command)
            .unwrap_or_else(|| panic!("no handle_{command}_command"));
        let entries = help_entries(&help_section(title));
        for flag in parsed_flags(body) {
            let excepted = EXCEPTIONS
                .iter()
                .any(|&(c, f, _)| (c == command || c == "*") && f == flag);
            if !excepted && !entries.contains(&flag) {
                missing.push(format!("{command} {flag}"));
            }
        }
    }
    assert!(
        missing.is_empty(),
        "parsed by the browser terminal but missing from its help: {missing:?}"
    );
}

#[test]
fn every_excepted_flag_is_named_in_its_commands_help() {
    for &(command, flag, why) in EXCEPTIONS {
        if command == "*" {
            continue;
        }
        let title = TITLES
            .iter()
            .find(|(c, _)| *c == command)
            .map(|(_, t)| *t)
            .unwrap_or_else(|| panic!("no title for {command}"));
        assert!(
            help_section(title).contains(flag),
            "{command} {flag} is excepted as {why:?}, but its help never names it"
        );
    }
}

#[test]
fn the_guard_reads_the_terminal_it_guards() {
    // Its own detectors, on the source as it stands, so a reformatted wasm.rs cannot make the
    // guard pass by finding nothing.
    let handlers = handlers();
    assert_eq!(handlers.len(), TITLES.len(), "one help title per handler");
    let trajectory = &handlers.iter().find(|(c, _)| c == "trajectory").unwrap().1;
    let parsed = parsed_flags(trajectory);
    for flag in [
        "--zero-humidity",
        "--zero-temperature",
        "--humidity",
        "--velocity",
    ] {
        assert!(parsed.contains(flag), "{flag} is parsed: {parsed:?}");
    }
    let entries = help_entries(&help_section("Trajectory"));
    for flag in [
        "--humidity",
        "--velocity",
        "--auto-zero",
        "--zero-density-altitude",
    ] {
        assert!(entries.contains(flag), "{flag} is an entry: {entries:?}");
    }
    assert!(help_entries(&help_section("True Velocity")).contains("--output"));
    // A flag named only inside another entry's description is no entry, even at the start of
    // a continuation line.
    assert_eq!(
        help_entries("    --a <X>   see --b\n                                 --c, the old name\n"),
        ["--a".to_string()].into_iter().collect()
    );
}
