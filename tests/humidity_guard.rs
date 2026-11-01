//! Source guard for the humidity-consistency train (plan C2): no in-crate production path may
//! compute its own per-step air, fly dry sound, or pair a density with an RH of its own.
//!
//! `src/atmosphere.rs` owns humid air: `StationAir`, `KernelAtmosphere` and `air_at` are the
//! single source of truth, and its legacy wrappers stay there. This test reads the PRODUCTION
//! code of every other file under `src/` (comments, string contents and `#[cfg(test)]` items
//! blanked; test-only module files skipped) and fails on:
//!
//! - a call to the legacy air functions (`get_local_atmosphere*`, `moist_speed_of_sound`) or to
//!   the legacy-tuple escape hatches (`StationAir::with_density_ratio`,
//!   `KernelAtmosphere::from_legacy_params`), which exist only for bindings' 4-vectors;
//! - a second dry gas constant (`401.87...`, `1.4 * 287...`);
//! - the 0.49.0 patterns that put a density ratio into the humidity field, and the
//!   "humidity not available" note that went with the dry literal (that one is searched in
//!   comments too);
//! - a per-step acceleration file (`cli_api.rs`, `derivatives.rs`, `fast_trajectory.rs`) that
//!   does not call `air_at`;
//! - a fixed sound speed (`340.29`, `343.0`) outside its allowlisted, documented sites.
//!
//! Together with `-D warnings` on the deprecated dry `get_local_atmosphere`, no in-crate path can
//! silently go dry.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

fn collect_rs_files(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(entries) = fs::read_dir(dir) else {
        return;
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

/// `text` with every comment and every string or char literal's contents replaced by spaces,
/// newlines kept, so byte offsets and line numbers survive.
fn code_only(text: &str) -> String {
    let bytes = text.as_bytes();
    let mut out = bytes.to_vec();
    let blank = |out: &mut Vec<u8>, from: usize, to: usize| {
        for byte in &mut out[from..to] {
            if *byte != b'\n' {
                *byte = b' ';
            }
        }
    };
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'/' if bytes.get(i + 1) == Some(&b'/') => {
                let end = text[i..].find('\n').map_or(bytes.len(), |n| i + n);
                blank(&mut out, i, end);
                i = end;
            }
            b'/' if bytes.get(i + 1) == Some(&b'*') => {
                let (mut depth, mut j) = (1, i + 2);
                while j < bytes.len() && depth > 0 {
                    if bytes[j] == b'/' && bytes.get(j + 1) == Some(&b'*') {
                        depth += 1;
                        j += 2;
                    } else if bytes[j] == b'*' && bytes.get(j + 1) == Some(&b'/') {
                        depth -= 1;
                        j += 2;
                    } else {
                        j += 1;
                    }
                }
                blank(&mut out, i, j);
                i = j;
            }
            b'r' if (bytes.get(i + 1) == Some(&b'"') || bytes.get(i + 1) == Some(&b'#'))
                && (i == 0 || !(bytes[i - 1].is_ascii_alphanumeric() || bytes[i - 1] == b'_')) =>
            {
                // A raw string r"..." / r#"..."#, or an identifier starting with r#.
                let mut j = i + 1;
                let mut hashes = 0;
                while bytes.get(j) == Some(&b'#') {
                    hashes += 1;
                    j += 1;
                }
                if bytes.get(j) != Some(&b'"') {
                    i += 1;
                    continue;
                }
                let terminator: Vec<u8> = std::iter::once(b'"')
                    .chain(std::iter::repeat_n(b'#', hashes))
                    .collect();
                let body = j + 1;
                let end = bytes[body..]
                    .windows(terminator.len())
                    .position(|w| w == terminator.as_slice())
                    .map_or(bytes.len(), |n| body + n);
                blank(&mut out, body, end);
                i = end + terminator.len();
            }
            b'"' => {
                let mut j = i + 1;
                while j < bytes.len() && bytes[j] != b'"' {
                    j += if bytes[j] == b'\\' { 2 } else { 1 };
                }
                blank(&mut out, i + 1, j.min(bytes.len()));
                i = j + 1;
            }
            b'\'' => {
                // A char literal ('x', '\n', '\u{..}'), or a lifetime ('a).
                if bytes.get(i + 1) == Some(&b'\\') {
                    let end = text[i + 2..].find('\'').map_or(bytes.len(), |n| i + 2 + n);
                    blank(&mut out, i + 1, end);
                    i = end + 1;
                } else if let Some(c) = text[i + 1..].chars().next() {
                    let after = i + 1 + c.len_utf8();
                    if bytes.get(after) == Some(&b'\'') {
                        blank(&mut out, i + 1, after);
                        i = after + 1;
                    } else {
                        i += 1;
                    }
                } else {
                    i += 1;
                }
            }
            _ => i += 1,
        }
    }
    String::from_utf8(out).expect("blanking keeps UTF-8 (only ASCII spaces replace bytes)")
}

/// Whether an attribute's text marks a test-only item.
fn is_test_only_cfg(attribute: &str) -> bool {
    let compact: String = attribute.chars().filter(|c| !c.is_whitespace()).collect();
    compact == "#[cfg(test)]" || compact.starts_with("#[cfg(all(test,")
}

/// The byte ranges of every `#[cfg(test)]` item in `code` (already `code_only`): from the
/// attribute to the `;` or the matching `}` that ends the item.
fn test_item_ranges(code: &str) -> Vec<(usize, usize)> {
    let bytes = code.as_bytes();
    let mut ranges = Vec::new();
    let mut search = 0;
    while let Some(found) = code[search..].find("#[cfg(") {
        let start = search + found;
        let Some(close) = code[start..].find(']') else {
            break;
        };
        let attribute_end = start + close + 1;
        if !is_test_only_cfg(&code[start..attribute_end]) {
            search = attribute_end;
            continue;
        }
        // Walk the item: `;` at nesting depth 0 ends it, or the `}` matching its first `{`.
        let (mut parens, mut braces, mut j) = (0i32, 0i32, attribute_end);
        let mut end = bytes.len();
        while j < bytes.len() {
            match bytes[j] {
                b'(' | b'[' => parens += 1,
                b')' | b']' => parens -= 1,
                b'{' => braces += 1,
                b'}' => {
                    braces -= 1;
                    if braces == 0 {
                        end = j + 1;
                        break;
                    }
                }
                b';' if parens == 0 && braces == 0 => {
                    end = j + 1;
                    break;
                }
                _ => {}
            }
            j += 1;
        }
        ranges.push((start, end));
        search = end;
    }
    ranges
}

/// `code` (already `code_only`) with every `#[cfg(test)]` item blanked, newlines kept.
fn without_test_items(code: &str) -> String {
    let mut out = code.as_bytes().to_vec();
    for (start, end) in test_item_ranges(code) {
        for byte in &mut out[start..end] {
            if *byte != b'\n' {
                *byte = b' ';
            }
        }
    }
    String::from_utf8(out).expect("ASCII blanking")
}

/// Files declared `#[cfg(test)] mod name;`: test-only through and through.
fn test_only_module_files(files: &[PathBuf]) -> Vec<PathBuf> {
    let mut out = Vec::new();
    for path in files {
        let Ok(text) = fs::read_to_string(path) else {
            continue;
        };
        let code = code_only(&text);
        let mut search = 0;
        while let Some(found) = code[search..].find("#[cfg(test)]") {
            let after = search + found + "#[cfg(test)]".len();
            let rest = code[after..].trim_start();
            if let Some(name) = rest
                .strip_prefix("mod ")
                .and_then(|r| r.split_once(';'))
                .map(|(name, _)| name.trim())
                .filter(|name| name.chars().all(|c| c.is_alphanumeric() || c == '_'))
            {
                let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");
                let dir = path.parent().unwrap();
                let base = if matches!(stem, "lib" | "main" | "mod") {
                    dir.to_path_buf()
                } else {
                    dir.join(stem)
                };
                out.push(base.join(format!("{name}.rs")));
                out.push(base.join(name).join("mod.rs"));
            }
            search = after;
        }
    }
    out
}

/// Every file under `src/` except `src/atmosphere.rs`, its children and the test-only module
/// files, as `(path relative to the crate, raw text, production code)`; the production code has
/// the raw text's byte offsets and lines.
fn production_sources() -> Vec<(String, String, String)> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut files = Vec::new();
    collect_rs_files(&root.join("src"), &mut files);
    assert!(files.len() > 50, "sanity: found {} files", files.len());
    let test_only = test_only_module_files(&files);
    let mut out = Vec::new();
    for path in files {
        let relative = path
            .strip_prefix(root)
            .unwrap()
            .to_string_lossy()
            .replace('\\', "/");
        if relative == "src/atmosphere.rs"
            || relative.starts_with("src/atmosphere/")
            || test_only.contains(&path)
        {
            continue;
        }
        let Ok(text) = fs::read_to_string(&path) else {
            continue;
        };
        let production = without_test_items(&code_only(&text));
        out.push((relative, text, production));
    }
    out
}

fn line_of(text: &str, offset: usize) -> usize {
    text[..offset].matches('\n').count() + 1
}

/// Offsets of `needle` in `code` that are a call (`needle` then optional spaces then `(`) and not
/// the tail of a longer identifier.
fn calls(code: &str, needle: &str) -> Vec<usize> {
    code.match_indices(needle)
        .filter(|(at, _)| {
            let before_ok = *at == 0 || {
                let b = code.as_bytes()[at - 1];
                !(b.is_ascii_alphanumeric() || b == b'_')
            };
            let after = code[at + needle.len()..].trim_start();
            before_ok && after.starts_with('(')
        })
        .map(|(at, _)| at)
        .collect()
}

/// Offsets of a numeric literal `needle` not followed by another digit.
fn literal(code: &str, needle: &str) -> Vec<usize> {
    code.match_indices(needle)
        .filter(|(at, _)| {
            let next = code.as_bytes().get(at + needle.len());
            let prev = if *at == 0 {
                None
            } else {
                code.as_bytes().get(at - 1)
            };
            !next.is_some_and(u8::is_ascii_digit)
                && !prev.is_some_and(|b| b.is_ascii_digit() || *b == b'.')
        })
        .map(|(at, _)| at)
        .collect()
}

#[test]
fn the_guard_sees_production_code_only() {
    // The blanking must keep code and drop comments, strings and test items, or the guard below
    // checks the wrong text.
    let sample = "fn a() { get_local_atmosphere(1.0) } // get_local_atmosphere(\n\
                  let s = \"moist_speed_of_sound(\"; let c = '{';\n\
                  #[cfg(test)]\nmod tests { fn b() { moist_speed_of_sound(1.0); } }\n\
                  #[cfg(test)]\nfn sized() -> [u8; 3] { with_density_ratio(1) }\n\
                  fn after() { from_legacy_params(2) }\n";
    let production = without_test_items(&code_only(sample));
    assert_eq!(calls(&production, "get_local_atmosphere").len(), 1);
    assert!(calls(&production, "moist_speed_of_sound").is_empty());
    assert!(calls(&production, "with_density_ratio").is_empty());
    assert_eq!(calls(&production, "from_legacy_params").len(), 1);
    assert_eq!(production.lines().count(), sample.lines().count());
}

#[test]
fn no_production_path_outside_atmosphere_rs_computes_its_own_air() {
    let banned_calls = [
        "get_local_atmosphere",
        "get_local_atmosphere_humid",
        "get_local_atmosphere_humid_with_gravity",
        "moist_speed_of_sound",
        "with_density_ratio",
        "from_legacy_params",
    ];
    let banned_code = [
        "401.87",
        "1.4 * 287",
        "1.4*287",
        "humidity: params.atmos_params",
        ".humidity = params.atmo_params",
    ];
    let mut offenders = Vec::new();
    for (path, text, production) in production_sources() {
        for needle in banned_calls {
            for at in calls(&production, needle) {
                offenders.push(format!("{path}:{}: calls {needle}", line_of(&text, at)));
            }
        }
        for needle in banned_code {
            for (at, _) in production.match_indices(needle) {
                offenders.push(format!("{path}:{}: {needle}", line_of(&text, at)));
            }
        }
        // The note that went with the dry literal: searched in comments too.
        let test_items = test_item_ranges(&code_only(&text));
        for (at, _) in text.match_indices("humidity not available") {
            if !test_items.iter().any(|&(from, to)| from <= at && at < to) {
                offenders.push(format!(
                    "{path}:{}: humidity not available",
                    line_of(&text, at)
                ));
            }
        }
    }
    assert!(
        offenders.is_empty(),
        "humid air is computed only in src/atmosphere.rs (StationAir, KernelAtmosphere, \
         air_at); a production path here computes its own, goes dry, or pairs a density ratio \
         with a humidity of its own:\n{}",
        offenders.join("\n")
    );
}

#[test]
fn every_per_step_acceleration_file_calls_air_at() {
    let required = [
        "src/cli_api.rs",
        "src/derivatives.rs",
        "src/fast_trajectory.rs",
    ];
    let sources = production_sources();
    for file in required {
        let (_, _, production) = sources
            .iter()
            .find(|(path, _, _)| path == file)
            .unwrap_or_else(|| panic!("{file} is missing"));
        let uses =
            calls(production, "air_at").len() + calls(production, "air_at_with_gravity").len();
        assert!(
            uses > 0,
            "{file} computes a per-step acceleration and must take its air from air_at"
        );
    }
}

/// Offsets of the identifier `name` in `code` (not part of a longer identifier), skipping
/// `use` declarations.
fn identifier_uses(code: &str, name: &str) -> Vec<usize> {
    let is_ident = |b: u8| b.is_ascii_alphanumeric() || b == b'_';
    code.match_indices(name)
        .filter(|(at, _)| {
            let bytes = code.as_bytes();
            let before = *at == 0 || !is_ident(bytes[at - 1]);
            let after = bytes.get(at + name.len()).is_none_or(|b| !is_ident(*b));
            let line_start = code[..*at].rfind('\n').map_or(0, |n| n + 1);
            let in_use = code[line_start..].trim_start().starts_with("use ");
            before && after && !in_use
        })
        .map(|(at, _)| at)
        .collect()
}

#[test]
fn fixed_sound_speeds_appear_only_at_their_documented_sites() {
    // A fixed sound speed (the literals 340.29 and 343.0, or the ICAO constant
    // `SPEED_OF_SOUND_MPS`) has no air behind it, so it must never stand in for a flight's local
    // sound speed. These production sites are the documented exceptions, none of which flies a
    // trajectory: the constant's definition (constants.rs); the legacy
    // `aerodynamic_jump::calculate_aerodynamic_jump` estimator, which no solver calls (they fly
    // the Litz jump); `angle_calculations::quick_drop_estimate`; the bc-convert default sound
    // speed (main.rs, wasm.rs) and the WASM BC-segment table's Mach labels (wasm.rs, two); and
    // the deprecated MBA-1134 spin-drift helpers (spin_drift.rs, three).
    let allowed: BTreeMap<&str, usize> = BTreeMap::from([
        ("src/aerodynamic_jump.rs", 1),
        ("src/angle_calculations.rs", 1),
        ("src/constants.rs", 1),
        ("src/main.rs", 1),
        ("src/spin_drift.rs", 3),
        ("src/wasm.rs", 3),
    ]);
    let mut found: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (path, text, production) in production_sources() {
        let mut offsets: Vec<usize> = ["340.29", "343.0"]
            .iter()
            .flat_map(|needle| literal(&production, needle))
            .collect();
        offsets.extend(identifier_uses(&production, "SPEED_OF_SOUND_MPS"));
        for at in offsets {
            found
                .entry(path.clone())
                .or_default()
                .push(line_of(&text, at));
        }
    }
    let mut offenders = Vec::new();
    for (path, lines) in &mut found {
        // One site per line: the constant's definition is both the literal and the name.
        lines.sort_unstable();
        lines.dedup();
        if allowed.get(path.as_str()) != Some(&lines.len()) {
            offenders.push(format!(
                "{path}: {} fixed sound speed(s) at lines {lines:?}, {} allowed",
                lines.len(),
                allowed.get(path.as_str()).copied().unwrap_or(0)
            ));
        }
    }
    for (path, count) in &allowed {
        if !found.contains_key(*path) {
            offenders.push(format!(
                "{path}: allowlisted for {count} but has none (update the list)"
            ));
        }
    }
    assert!(offenders.is_empty(), "{}", offenders.join("\n"));
}
