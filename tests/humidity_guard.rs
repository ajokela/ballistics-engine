//! Source guard for the 0.50.0 humidity work: no in-crate production path may
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
//! - a fixed sound speed (`340.29`, `343.0`) outside its allowlisted, documented sites;
//! - a humidity scaled by 100 or clamped into a 0-1 fraction by hand: a Rust
//!   humidity is a percent, and a wire fraction converts through `RelativeHumidity`;
//! - a humidity bounded by hand: a range, an ordering comparison, `min`/`max` or a
//!   low-high argument pair against 0, 1, 100 or a sensor-band edge, on a line or in an
//!   expression that names a humidity (in a string too, such as a wire path). The one range is
//!   `RelativeHumidity`'s.
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
    blank_non_code(text, true)
}

/// `text` with every comment replaced by spaces but string and char literals KEPT, newlines
/// kept, so byte offsets and line numbers match [`code_only`]'s. For telling whether code NAMES
/// a humidity, which a wire path (`"$.atmosphere.relative_humidity"`) or a flag name does from
/// inside a string.
fn without_comments(text: &str) -> String {
    blank_non_code(text, false)
}

/// [`code_only`] and [`without_comments`]: comments always blanked, string and char literal
/// contents blanked when `blank_literals`.
fn blank_non_code(text: &str, blank_literals: bool) -> String {
    let bytes = text.as_bytes();
    let mut out = bytes.to_vec();
    let blank = |out: &mut Vec<u8>, from: usize, to: usize| {
        for byte in &mut out[from..to] {
            if *byte != b'\n' {
                *byte = b' ';
            }
        }
    };
    let blank_literal = |out: &mut Vec<u8>, from: usize, to: usize| {
        if blank_literals {
            blank(out, from, to);
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
                blank_literal(&mut out, body, end);
                i = end + terminator.len();
            }
            b'"' => {
                let mut j = i + 1;
                while j < bytes.len() && bytes[j] != b'"' {
                    j += if bytes[j] == b'\\' { 2 } else { 1 };
                }
                blank_literal(&mut out, i + 1, j.min(bytes.len()));
                i = j + 1;
            }
            b'\'' => {
                // A char literal ('x', '\n', '\u{..}'), or a lifetime ('a).
                if bytes.get(i + 1) == Some(&b'\\') {
                    let end = text[i + 2..].find('\'').map_or(bytes.len(), |n| i + 2 + n);
                    blank_literal(&mut out, i + 1, end);
                    i = end + 1;
                } else if let Some(c) = text[i + 1..].chars().next() {
                    let after = i + 1 + c.len_utf8();
                    if bytes.get(after) == Some(&b'\'') {
                        blank_literal(&mut out, i + 1, after);
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
    without_test_items_of(code, code)
}

/// `text` (the same file as `code`, `code_only`, byte for byte) with every `#[cfg(test)]` item
/// that `code` holds blanked, newlines kept. Each item runs from its ASCII `#` to its ASCII `;`
/// or `}`, so whole characters are blanked and `text` stays UTF-8.
fn without_test_items_of(text: &str, code: &str) -> String {
    let mut out = text.as_bytes().to_vec();
    for (start, end) in test_item_ranges(code) {
        for byte in &mut out[start..end] {
            if *byte != b'\n' {
                *byte = b' ';
            }
        }
    }
    String::from_utf8(out).expect("blanking whole characters keeps UTF-8")
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

#[test]
fn only_station_air_decides_whether_air_can_exist() {
    // One impossible-air rule. Whether a temperature, pressure and humidity
    // describe air that can exist (an uncapped water-vapour mole fraction below 1, and a
    // temperature below water's critical point) is decided by `StationAir::new` alone, and every
    // surface that refuses such air asks it. 0.49.0's `atmosphere.density_altitude` service
    // evaluated the vapour fraction itself, so its refusal and the solve's were two copies of
    // one rule that could drift apart.
    let predicates = [
        "water_vapor_mole_fraction_uncapped",
        "moist_air_vapor_mole_fraction",
        "vapor_mole_fraction",
        "enhanced_saturation_vapor_pressure",
    ];
    let mut offenders = Vec::new();
    for (path, text, production) in production_sources() {
        for needle in predicates {
            for at in calls(&production, needle) {
                offenders.push(format!("{path}:{}: calls {needle}", line_of(&text, at)));
            }
        }
    }
    assert!(
        offenders.is_empty(),
        "whether air can exist is StationAir::new's to decide; refuse impossible air through \
         it rather than evaluating the vapour fraction here:\n{}",
        offenders.join("\n")
    );
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

/// Offsets in `line` (production code) of a scaling by 100: `* 100`, `/ 100`, `100 *`, `100 /`,
/// with or without `.0` and spaces, the `100` not part of a longer number.
fn scalings_by_100(line: &str) -> Vec<usize> {
    let bytes = line.as_bytes();
    let is_num = |b: u8| b.is_ascii_digit() || b == b'.' || b == b'_';
    let mut out = Vec::new();
    for (at, _) in line.match_indices("100") {
        if at > 0 && is_num(bytes[at - 1]) {
            continue;
        }
        let mut end = at + 3;
        if line[end..].starts_with(".0") {
            end += 2;
            while bytes.get(end) == Some(&b'0') {
                end += 1;
            }
        }
        if bytes.get(end).is_some_and(|b| is_num(*b)) {
            continue;
        }
        let before = line[..at].trim_end();
        let after = line[end..].trim_start();
        let operator_before = before.ends_with('*') || before.ends_with('/');
        let operator_after = (after.starts_with('*') || after.starts_with('/'))
            && !after.starts_with("//")
            && !after.starts_with("/*");
        if operator_before || operator_after {
            out.push(at);
        }
    }
    out
}

#[test]
fn the_scaling_finder_finds_scalings_only() {
    for line in [
        "humidity / 100.0",
        "h*100.0",
        "x * 100",
        "100.0 * f",
        "100.00/ h",
        "(h / 100.0).clamp(0.0, 1.0)",
    ] {
        assert_eq!(scalings_by_100(line).len(), 1, "{line}");
    }
    for line in [
        "humidity: 100.0,",
        "1000.0 * h",
        "h / 1000.0",
        "h / 100.5",
        "0.100 * h",
        "x = 100.0;",
        "let pct = 100.0 - h;",
        "1_100.0 / x",
    ] {
        assert!(scalings_by_100(line).is_empty(), "{line}");
    }
}

#[test]
fn no_production_path_converts_a_humidity_by_hand() {
    // One unit per contract. Rust types and functions carry a relative humidity as a
    // PERCENT (`BallisticInputs.humidity_percent`, `AtmosphericConditions.humidity`), and the
    // fraction<->percent conversion is `RelativeHumidity::try_fraction` / `fraction()`, so no
    // production line outside src/atmosphere.rs scales a humidity by 100 itself, or clamps one
    // into a 0-1 fraction. 0.49.0 did both at every site that wrote the old fraction field.
    //
    // No site is allowed: the last one, the `atmosphere.density_altitude` service's SI-wire
    // fraction, converts through `RelativeHumidity::try_fraction` since 0.50.0.
    let allowed: BTreeMap<&str, usize> = BTreeMap::new();
    let mut found: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for (path, text, production) in production_sources() {
        for (index, line) in production.lines().enumerate() {
            if !line.to_ascii_lowercase().contains("humid") {
                continue;
            }
            let scales = !scalings_by_100(line).is_empty();
            let clamps_to_a_fraction = line.contains("clamp(0.0, 1.0)");
            if scales || clamps_to_a_fraction {
                let raw = text.lines().nth(index).unwrap_or("").trim();
                found
                    .entry(path.clone())
                    .or_default()
                    .push(format!("{path}:{}: {raw}", index + 1));
            }
        }
    }
    let mut offenders = Vec::new();
    for (path, sites) in &found {
        if allowed.get(path.as_str()) != Some(&sites.len()) {
            offenders.extend(sites.iter().cloned());
        }
    }
    for (path, count) in &allowed {
        if !found.contains_key(*path) {
            offenders.push(format!(
                "{path}: allowlisted for {count} but has none (update the list)"
            ));
        }
    }
    assert!(
        offenders.is_empty(),
        "a humidity is converted by hand; carry a percent, and convert a wire fraction with \
         RelativeHumidity::try_fraction:\n{}",
        offenders.join("\n")
    );
}

/// The byte range of the expression a literal at `at` belongs to: its field, argument list or
/// statement. Balanced `()`/`[]` groups are read as one unit; a literal inside a call's
/// arguments belongs to the whole call (`map_or(50.0, |p| p.humidity)`); the expression ends at
/// a `,`, `;`, `{` or `}` outside every group; and a block's tail expression
/// (`fn default_humidity() -> f64 { 50.0 }`) includes the item header before the block.
fn enclosing_expression(code: &str, at: usize) -> (usize, usize) {
    let b = code.as_bytes();
    // The outermost group enclosing `at` within its statement: the anchor.
    let (mut anchor_left, mut depth, mut i) = (at, 0i32, at);
    while i > 0 {
        i -= 1;
        match b[i] {
            b')' | b']' => depth += 1,
            b'(' | b'[' if depth > 0 => depth -= 1,
            b'(' | b'[' => anchor_left = i,
            b';' | b'{' | b'}' if depth == 0 => break,
            _ => {}
        }
    }
    let matching_close = |open: usize| {
        let mut depth = 0i32;
        for (j, &c) in b.iter().enumerate().skip(open) {
            match c {
                b'(' | b'[' => depth += 1,
                b')' | b']' => {
                    depth -= 1;
                    if depth == 0 {
                        return j + 1;
                    }
                }
                _ => {}
            }
        }
        b.len()
    };
    let anchor_right = if anchor_left < at {
        matching_close(anchor_left)
    } else {
        at
    };
    // Outward from the anchor, groups opaque, to the delimiters.
    let (mut start, mut left_stop, mut depth, mut i) = (0, b'\0', 0i32, anchor_left);
    while i > 0 {
        i -= 1;
        match b[i] {
            b')' | b']' => depth += 1,
            b'(' | b'[' => depth -= 1,
            c @ (b',' | b';' | b'{' | b'}') if depth == 0 => {
                start = i + 1;
                left_stop = c;
                break;
            }
            _ => {}
        }
    }
    let (mut end, mut right_stop, mut depth) = (b.len(), b'\0', 0i32);
    for (j, &c) in b.iter().enumerate().skip(anchor_right) {
        match c {
            b'(' | b'[' => depth += 1,
            b')' | b']' => depth -= 1,
            b',' | b';' | b'{' | b'}' if depth <= 0 => {
                end = j;
                right_stop = c;
                break;
            }
            _ => {}
        }
    }
    if left_stop == b'{' && right_stop == b'}' {
        start = code[..start - 1]
            .rfind([';', '{', '}'])
            .map_or(0, |n| n + 1);
    }
    (start, end)
}

/// The literals a humidity default is written as: 50 % as a percent (`50.0`), or as the SI
/// wire's fraction (`0.5`, `0.50`).
const HUMIDITY_DEFAULT_SPELLINGS: [&str; 3] = ["50.0", "0.5", "0.50"];

/// Offsets of a humidity-default literal ([`HUMIDITY_DEFAULT_SPELLINGS`], no further digits) in
/// `code` whose expression (see [`enclosing_expression`]) names a humidity, in offset order.
fn humidity_default_literals(code: &str) -> Vec<usize> {
    let mut found: Vec<usize> = HUMIDITY_DEFAULT_SPELLINGS
        .iter()
        .flat_map(|spelling| literal(code, spelling))
        .filter(|&at| {
            let (start, end) = enclosing_expression(code, at);
            code[start..end].to_ascii_lowercase().contains("humid")
        })
        .collect();
    found.sort_unstable();
    found
}

#[test]
fn the_default_finder_finds_humidity_defaults_only() {
    for code in [
        "let mut humidity = 50.0;",
        "let h = humidity.or_else(|| p.humidity)\n    .unwrap_or(50.0);",
        "saved.map_or(50.0, |p| p.humidity)",
        "fn default_humidity() -> f64 {\n    50.0\n}",
        "Inputs { humidity_percent: 50.0, x: 1.0 }",
        "csv(&data, &[\"H\"], saved.map_or(50.0, |p| p.humidity))",
        "let mut humidity: f64 = 50.0;",
        // A fraction default (solve v1's SI wire) is the same default in the other unit.
        "const DEFAULT_RELATIVE_HUMIDITY: f64 = 0.5;",
        "let rh = request.relative_humidity.unwrap_or(0.5);",
        "AtmosphereV1 { relative_humidity: 0.50, x: 1.0 }",
        "fn default_humidity_fraction() -> f64 {\n    0.5\n}",
    ] {
        assert_eq!(humidity_default_literals(code).len(), 1, "{code}");
    }
    for code in [
        "let wind = 50.0; let humidity = h;",
        "let humidity = 150.0;",
        "let x = 50.05; let humidity = 1.0;",
        "fn f() { 50.0 } fn humidity() {}",
        "Inputs { tipoff_decay_distance: 50.0, humidity_percent: h }",
        "call(a, 50.0); let humidity = h;",
        "let half = 0.5; let humidity = h;",
        "let humidity = 0.55;",
        "let humidity = 10.5;",
        "let humidity = 100.5;",
        "let humidity = 0.505;",
        "Atmo { lateral: 0.5, relative_humidity: h }",
        "const DEFAULT_RELATIVE_HUMIDITY: f64 = RelativeHumidity::DEFAULT.fraction();",
    ] {
        assert!(humidity_default_literals(code).is_empty(), "{code}");
    }
}

#[test]
fn every_production_humidity_default_is_the_one_constant() {
    // One default. An omitted humidity means `DEFAULT_HUMIDITY_PERCENT` (and
    // `RelativeHumidity::DEFAULT`) on every surface, so no production line outside
    // src/atmosphere.rs writes its own 50.0 for one, nor its own 0.5 for the SI wire's fraction
    // (solve v1's default is `RelativeHumidity::DEFAULT.fraction()`). (clap's
    // `default_value = "50.0"` strings are string literals, blanked here; `src/main.rs`'s
    // `humidity_flag_tests` pins them.)
    let mut offenders = Vec::new();
    for (path, text, production) in production_sources() {
        for at in humidity_default_literals(&production) {
            let line = line_of(&text, at);
            let raw = text.lines().nth(line - 1).unwrap_or("").trim();
            offenders.push(format!("{path}:{line}: {raw}"));
        }
    }
    assert!(
        offenders.is_empty(),
        "a humidity default is a bare 50.0 or 0.5; use DEFAULT_HUMIDITY_PERCENT, or \
         RelativeHumidity::DEFAULT.fraction() for a fraction:\n{}",
        offenders.join("\n")
    );
}

#[test]
fn no_production_path_clamps_a_humidity_into_range() {
    // Out-of-range input is REJECTED with a field-named error, never clamped. The
    // only clamps left are inside src/atmosphere.rs: the sensor band's snap of a reading to the
    // end it reads as (`RelativeHumidity`'s one rule), and the total-function guards behind
    // validated input and the legacy wrappers.
    let mut offenders = Vec::new();
    for (path, text, production) in production_sources() {
        for (index, line) in production.lines().enumerate() {
            let lower = line.to_ascii_lowercase();
            if !lower.contains("humid") {
                continue;
            }
            let compact: String = line.chars().filter(|c| !c.is_whitespace()).collect();
            if compact.contains("clamp(0.0,100.0)")
                || compact.contains(".min(100.0)")
                || compact.contains(".max(0.0).min(")
            {
                let raw = text.lines().nth(index).unwrap_or("").trim();
                offenders.push(format!("{path}:{}: {raw}", index + 1));
            }
        }
    }
    assert!(
        offenders.is_empty(),
        "a humidity is clamped; refuse it through RelativeHumidity instead:\n{}",
        offenders.join("\n")
    );
}

/// The values a hand-written humidity bound is written with, by magnitude: dry and saturated
/// air in either unit (0, 1 and 100) and the sensor band's edges in either unit (5, 105, 0.05
/// and 1.05).
const HUMIDITY_BOUND_VALUES: [f64; 7] = [0.0, 1.0, 100.0, 5.0, 105.0, 0.05, 1.05];

/// The FLOAT literal that starts at `at` in `code` (`0.0`, `1e2`, `100_f64`, `1.05`; not the
/// tail of an identifier, a field access or a longer number): its end offset and its value.
/// Every humidity is an `f64`, so a bound on one is a float literal, and an integer (`0..n`,
/// `len() > 0`) is never one.
fn float_literal_at(code: &str, at: usize) -> Option<(usize, f64)> {
    let b = code.as_bytes();
    if !b.get(at)?.is_ascii_digit() {
        return None;
    }
    // `x.0` is a field; `..1.0` (a range's end) is a literal.
    let after_a_field_dot = at > 0 && b[at - 1] == b'.' && !(at > 1 && b[at - 2] == b'.');
    if after_a_field_dot || (at > 0 && (b[at - 1].is_ascii_alphanumeric() || b[at - 1] == b'_')) {
        return None;
    }
    let digits = |mut j: usize| {
        while b.get(j).is_some_and(|c| c.is_ascii_digit() || *c == b'_') {
            j += 1;
        }
        j
    };
    let mut end = digits(at);
    let mut float = false;
    if b.get(end) == Some(&b'.') && b.get(end + 1).is_some_and(u8::is_ascii_digit) {
        end = digits(end + 1);
        float = true;
    }
    if matches!(b.get(end), Some(b'e' | b'E')) {
        let sign = usize::from(matches!(b.get(end + 1), Some(b'+' | b'-')));
        if b.get(end + 1 + sign).is_some_and(u8::is_ascii_digit) {
            end = digits(end + 1 + sign);
            float = true;
        }
    }
    let number_end = end;
    while b
        .get(end)
        .is_some_and(|c| c.is_ascii_alphanumeric() || *c == b'_')
    {
        end += 1;
    }
    match code[number_end..end].trim_start_matches('_') {
        "" => {}
        "f64" | "f32" => float = true,
        _ => return None,
    }
    if !float {
        return None;
    }
    let value = code[at..number_end].replace('_', "").parse().ok()?;
    Some((end, value))
}

/// The float literal that ENDS at `end` in `code`, whitespace before `end` skipped: its start
/// and its value.
fn float_literal_ending_at(code: &str, end: usize) -> Option<(usize, f64)> {
    let end = code[..end].trim_end().len();
    let b = code.as_bytes();
    let mut start = end;
    while start > 0
        && (b[start - 1].is_ascii_alphanumeric() || b[start - 1] == b'_' || b[start - 1] == b'.')
    {
        start -= 1;
    }
    let (literal_end, value) = float_literal_at(code, start)?;
    (literal_end == end).then_some((start, value))
}

/// The float literal that BEGINS the operand at `from` in `code`, whitespace and one unary
/// minus skipped: its start and its value.
fn float_literal_after(code: &str, from: usize) -> Option<(usize, f64)> {
    let skip_spaces = |at: usize| at + (code[at..].len() - code[at..].trim_start().len());
    let mut at = skip_spaces(from);
    if code.as_bytes().get(at) == Some(&b'-') {
        at = skip_spaces(at + 1);
    }
    float_literal_at(code, at).map(|(_, value)| (at, value))
}

/// Offsets in `code` (production code, `code_only`) of the float literals that bound a value by
/// hand with one of the [`HUMIDITY_BOUND_VALUES`]: an end of a range
/// (`(0.0..=1.0).contains(&v)`, `(0.0..).contains(&v)`), an operand of an ordering comparison
/// (`v > 100.0`, `1.0 < v`, `v >= -0.05`), an argument of `min` or `max` (`v.min(100.0)`,
/// `f64::max(v, 0.0)`), or a pair of such arguments of any call, next to each other, low then
/// high (`v.clamp(0.0, 1.0)`, `require_range(path, v, 0.0, 1.0)`, `in_band(v, -5.0, 105.0)`;
/// a run of zeros, `f(.., 0.0, 0.0, ..)`, is not a range). Sorted, without repeats.
fn hand_written_bounds(code: &str) -> Vec<usize> {
    let b = code.as_bytes();
    let bound = |found: Option<(usize, f64)>| {
        found
            .filter(|(_, value)| HUMIDITY_BOUND_VALUES.contains(&value.abs()))
            .map(|(at, _)| at)
    };
    let mut out = Vec::new();
    // Ranges: `a..b` and `a..=b`, never `...`.
    for (at, _) in code.match_indices("..") {
        if (at > 0 && b[at - 1] == b'.') || b.get(at + 2) == Some(&b'.') {
            continue;
        }
        let operator_end = if b.get(at + 2) == Some(&b'=') {
            at + 3
        } else {
            at + 2
        };
        out.extend(bound(float_literal_ending_at(code, at)));
        out.extend(bound(float_literal_after(code, operator_end)));
    }
    // Ordering comparisons: `<`, `<=`, `>`, `>=`, never a shift, `->` or `=>`.
    for (i, &c) in b.iter().enumerate() {
        if c != b'<' && c != b'>' {
            continue;
        }
        let previous = i.checked_sub(1).map(|j| b[j]);
        let next = b.get(i + 1).copied();
        if next == Some(c)
            || previous == Some(c)
            || (c == b'>' && matches!(previous, Some(b'-' | b'=')))
        {
            continue;
        }
        let operator_end = if next == Some(b'=') { i + 2 } else { i + 1 };
        out.extend(bound(float_literal_ending_at(code, i)));
        out.extend(bound(float_literal_after(code, operator_end)));
    }
    // Calls: `min`/`max` with one bound literal argument, any call with a low-high pair.
    for (open, _) in code.match_indices('(') {
        let name_end = code[..open].trim_end().len();
        let name_start = code[..name_end]
            .rfind(|c: char| !(c.is_ascii_alphanumeric() || c == '_'))
            .map_or(0, |n| n + 1);
        let name = &code[name_start..name_end];
        if name.is_empty() || name.as_bytes()[0].is_ascii_digit() {
            continue;
        }
        let (mut depth, mut argument_start, mut arguments) = (0i32, open + 1, Vec::new());
        for (j, &c) in b.iter().enumerate().skip(open) {
            match c {
                b'(' | b'[' | b'{' => depth += 1,
                b')' | b']' | b'}' => {
                    depth -= 1;
                    if depth == 0 {
                        arguments.push((argument_start, j));
                        break;
                    }
                }
                b',' if depth == 1 => {
                    arguments.push((argument_start, j));
                    argument_start = j + 1;
                }
                _ => {}
            }
        }
        // Each argument that is a bound literal and nothing else: its offset and signed value.
        let literals: Vec<Option<(usize, f64)>> = arguments
            .into_iter()
            .map(|(start, end)| {
                let (at, value) = float_literal_after(code, start)?;
                let literal_end = float_literal_at(code, at)?.0;
                let negative = code[start..at].trim() == "-";
                (code[literal_end..end].trim().is_empty() && HUMIDITY_BOUND_VALUES.contains(&value))
                    .then_some((at, if negative { -value } else { value }))
            })
            .collect();
        if matches!(name, "min" | "max") {
            out.extend(literals.iter().flatten().map(|&(at, _)| at));
        }
        for pair in literals.windows(2) {
            if let [Some((low_at, low)), Some((high_at, high))] = *pair {
                if low < high {
                    out.extend([low_at, high_at]);
                }
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

/// Offsets of the [`hand_written_bounds`] in the PRODUCTION code of `text` (a whole source
/// file: comments, string contents and `#[cfg(test)]` items blanked) whose line, or whose
/// expression ([`enclosing_expression`]), names a humidity. The name may be inside a string: a
/// wire path (`"$.atmosphere.relative_humidity"`) names one.
fn hand_written_humidity_bounds(text: &str) -> Vec<usize> {
    let code = code_only(text);
    let production = without_test_items(&code);
    let names = without_test_items_of(&without_comments(text), &code).to_ascii_lowercase();
    hand_written_bounds(&production)
        .into_iter()
        .filter(|&at| {
            let line_start = names[..at].rfind('\n').map_or(0, |n| n + 1);
            let line_end = names[at..].find('\n').map_or(names.len(), |n| at + n);
            let (start, end) = enclosing_expression(&production, at);
            names[line_start..line_end].contains("humid") || names[start..end].contains("humid")
        })
        .collect()
}

#[test]
fn the_bound_finder_finds_hand_written_humidity_bounds_only() {
    for (code, count) in [
        // The perturbation layer's sensor-band domain check, as first written, on one line and on
        // two.
        (
            "if axis == InputAxis::RelativeHumidity && !(0.0..=1.0).contains(&value) {}",
            2,
        ),
        (
            "if axis == InputAxis::RelativeHumidity\n    && !(0.0..=1.0).contains(&value)\n{}",
            2,
        ),
        ("let ok = (0.0..=100.0).contains(&humidity);", 2),
        ("let ok = (-5.0..=105.0).contains(&humidity_percent);", 2),
        ("let ok = (0.0..).contains(&humidity);", 1),
        ("if humidity_percent > 100.0 {}", 1),
        ("if 100.0 < humidity {}", 1),
        ("let ok = humidity_fraction >= 1.05;", 1),
        ("let ok = relative_humidity < -0.05;", 1),
        ("let ok = humidity <= 1e2;", 1),
        ("let h = humidity.min(100.0);", 1),
        ("let h = humidity.max(0.0_f64);", 1),
        ("let h = f64::max(humidity, 0.0);", 1),
        ("let h = humidity.clamp(-5.0, 105.0);", 2),
        // A wire path names the humidity from inside a string: solve v1's 0.49 decoder check.
        (
            "require_range(\"$.atmosphere.relative_humidity\", value, 0.0, 1.0)?;",
            2,
        ),
        // A block's tail expression includes its item header.
        (
            "fn saturated(humidity_percent: f64) -> bool {\n    p > 100.0\n}",
            1,
        ),
        (
            "let ok = (..=1.0).contains(&humidity) && (..100.0).contains(&humidity_percent);",
            2,
        ),
    ] {
        assert_eq!(hand_written_humidity_bounds(code).len(), count, "{code}");
    }
    for code in [
        "humidity: Option<f64>,",
        "fn f(humidity: f64) -> Result<f64, String> { g(humidity) }",
        "let humidities: Vec<Vec<f64>> = Vec::new();",
        "for i in 0..humidities.len() {}",
        "if humidities.len() > 0 {}",
        "if humidity > threshold {}",
        "if humidity > 50.0 {}",
        "if humidity_percent > 1000.0 {}",
        "if humidity < 0.5 {}",
        "let x = humidity_percent - 100.0;",
        "let band = humidity_band.0..=humidity_band.1;",
        "let h = match humidity { Some(h) => h, None => 0.0 };",
        "let h = RelativeHumidity::try_percent(100.0);",
        "let h = RelativeHumidity::DEFAULT.fraction();",
        // The humidity is named only in a comment, or on another line and expression.
        "let x = distance.max(0.0); // humidity",
        "if speed > 0.0 {\n    humidity = h;\n}",
        "let wind = Vector3::new(0.0, 0.0, 1.0);\nlet humidity = h;",
        // Positional zeros in a call that also takes the humidity are not a range.
        "fly(velocity, humidity, 0.0, 0.0, 1.5)",
        "fly(velocity, humidity, 1.0, 0.0)",
        "#[cfg(test)]\nfn t() { assert!(humidity > 100.0); }",
    ] {
        assert!(hand_written_humidity_bounds(code).is_empty(), "{code}");
    }
}

#[test]
fn no_production_path_bounds_a_humidity_by_hand() {
    // The sensor band (owner decision 2026-11-01): a relative humidity has ONE range,
    // RelativeHumidity's in src/atmosphere.rs: the sensor band a reading may lie in
    // (HUMIDITY_SENSOR_TOLERANCE_PERCENT, `RelativeHumidity::FRACTION_BAND`) and the 0..=1 physical
    // humidity it snaps into (`RelativeHumidity::is_physical_fraction`). So no production line or
    // expression outside src/atmosphere.rs that names a humidity checks one against 0, 1, 100 or a
    // band edge itself (a range, an ordering comparison, `min`/`max`, or a call with two such
    // bounds). A second range would drift from the band; one did while the band was being
    // written (the perturbation layer's `(0.0..=1.0).contains`), and solve v1's 0.49
    // `require_range(.., 0.0, 1.0)` was another.
    let mut offenders = Vec::new();
    for (path, text, _) in production_sources() {
        let mut lines: Vec<usize> = hand_written_humidity_bounds(&text)
            .into_iter()
            .map(|at| line_of(&text, at))
            .collect();
        lines.dedup();
        for line in lines {
            let raw = text.lines().nth(line - 1).unwrap_or("").trim();
            offenders.push(format!("{path}:{line}: {raw}"));
        }
    }
    assert!(
        offenders.is_empty(),
        "a humidity is bounded by hand; take its range from RelativeHumidity \
         (try_percent/try_fraction, FRACTION_BAND, is_physical_fraction):\n{}",
        offenders.join("\n")
    );
}

#[test]
fn every_wasm_terminal_humidity_flag_goes_through_the_one_validator() {
    // 0.49.0's browser terminal range-checked `--humidity` only in true-velocity
    // and true-wind; every other command flew 150 % clamped. Every `--humidity`-family arm in
    // src/wasm.rs must parse its value with `humidity_value`, the terminal's wrapper over
    // `atmosphere::parse_humidity_percent_arg` (the one validator, natively tested).
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let text = fs::read_to_string(root.join("src/wasm.rs")).expect("read src/wasm.rs");
    let lines: Vec<&str> = text.lines().collect();
    let mut arms = 0;
    let mut offenders = Vec::new();
    for (index, line) in lines.iter().enumerate() {
        let trimmed = line.trim_start();
        let is_arm = (trimmed.starts_with("\"--humidity\"")
            || trimmed.starts_with("\"--zero-humidity\""))
            && trimmed.contains("=>");
        if !is_arm {
            continue;
        }
        arms += 1;
        let body = lines[index..(index + 10).min(lines.len())].join("\n");
        let arm_end = body[1..]
            .find("\n                \"-")
            .map_or(body.len(), |n| n + 1);
        if !body[..arm_end].contains("humidity_value(") {
            offenders.push(format!("src/wasm.rs:{}: {}", index + 1, trimmed));
        }
    }
    assert!(
        offenders.is_empty(),
        "a WASM humidity flag is parsed without the one validator:\n{}",
        offenders.join("\n")
    );
    // trajectory --humidity and --zero-humidity, zero, lead, monte-carlo, true-velocity,
    // true-wind and estimate-bc.
    assert!(arms >= 8, "sanity: found {arms} humidity arms");
    let shared = calls(&text, "parse_humidity_percent_arg").len();
    assert!(
        shared >= 1,
        "humidity_value must call parse_humidity_percent_arg"
    );
}
