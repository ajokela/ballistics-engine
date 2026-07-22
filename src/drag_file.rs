//! Cleanroom parser for the `.drg` drag-curve text file format (MBA-1409).
//!
//! `.drg` is a small vendor text format used to distribute Doppler-radar-measured
//! drag-coefficient-vs-Mach curves for individual bullets (e.g. Lapua's free downloads
//! for QuickTARGET Unlimited). This module was written from **publicly observable file
//! structure only** — line layout, field separators, and column semantics inferred by
//! inspecting one real vendor sample locally for grammar-pinning purposes (see the task
//! report for structural findings). **No vendor drag-coefficient data is embedded here**:
//! every fixture in this module's tests is synthetic, invented content.
//!
//! Shared by the CLI (`main.rs`) and the WASM terminal (`wasm.rs`); lives in the library
//! crate — not gated behind any feature — so it compiles for `wasm32-unknown-unknown`.
//! This module is fs-free: it parses an in-memory `&str`, never touches `std::fs`. File
//! I/O (reading the path the user gave `--drag-table`, or the WASM text-ingestion glue)
//! stays in the caller.
//!
//! ## Grammar (tolerant superset covering the vendor layout and plain two-column CSV)
//!
//! - Zero or more leading non-data lines (title/header/blank/comment) are skipped until
//!   the first line that parses as exactly two finite numbers. The first non-empty
//!   skipped line (trimmed) is captured as [`ParsedDragCurve::name`]; if the very first
//!   line in the file is itself a data row (headerless CSV), `name` is `None`.
//! - From that point on, every non-blank line must parse as exactly two finite numbers —
//!   garbage after data has started is a hard error, not silently skipped.
//! - Fields are separated by a run of whitespace (space or tab — real vendor files use
//!   tab) or by a single `,`/`;` delimiter.
//! - Columns may appear as `(mach, cd)` *or* `(cd, mach)`. Column detection: if exactly
//!   one column is strictly ascending across every row, that column is mach (this is not
//!   merely a defensive fallback — some real vendor files use `(cd, mach)` column order,
//!   so this detection is load-bearing). If *both* columns ascend — e.g. a `(cd, mach)`
//!   file whose cd happens to be monotonic across a subsonic-only deck — the column with
//!   the larger maximum value is mach: mach spans past `1.0` in real decks, while cd
//!   stays well under `1.5`. If the two maxima are within 20% of each other, that's too
//!   close to call, and parsing fails with a dedicated "ambiguous columns" error rather
//!   than silently guessing. If neither column ascends, parsing fails with a dedicated
//!   error naming the problem.
//! - The decimal separator is `.` only. A line is only flagged as the decimal-comma
//!   error if it splits into exactly two tokens and *both* look like decimal-comma
//!   numbers (e.g. `"0,5 0,3"`); this produces a dedicated error naming the problem,
//!   rather than a generic parse failure or a comma-CSV misread. A line where only one
//!   token looks like a decimal-comma number — e.g. a header line mixing prose and a
//!   comma-formatted number such as `"Density 1,225"` — is not treated as this error; it
//!   is simply skipped like ordinary header text.
//! - Row count must be in `2..=4096` (`4096` matches `ffi::MAX_FFI_DRAG_TABLE_LEN`, the
//!   cap already enforced on the array-based FFI drag-table entry point).
//! - `mach` must be finite and `>= 0` (the real vendor sample's first row is exactly
//!   `mach = 0`); `cd` must be finite and `> 0`.

/// The maximum number of `(mach, cd)` rows a `.drg` file may contain.
/// Matches `ffi::MAX_FFI_DRAG_TABLE_LEN`, the cap already enforced on the array-based FFI
/// custom-drag-table entry point, so a `.drg` load can never produce a table too large for
/// that path to accept.
const MAX_DRG_ROWS: usize = 4096;
/// The minimum number of rows needed to define a drag curve.
const MIN_DRG_ROWS: usize = 2;

/// A drag curve parsed from `.drg` text: an optional vendor-supplied name/description and
/// the `(mach, cd)` points in ascending-mach order.
#[derive(Debug, Clone, PartialEq)]
pub struct ParsedDragCurve {
    pub name: Option<String>,
    pub points: Vec<(f64, f64)>,
}

/// How a single line of `.drg` text classifies during parsing.
enum LineKind {
    /// Blank (whitespace-only) line: always skippable, before or after data starts.
    Blank,
    /// A data row: two finite-parseable tokens, in file column order (not yet assigned
    /// to mach/cd — that happens once every row has been collected).
    Row(f64, f64),
    /// Looks like an attempted data row using a decimal comma (e.g. `"0,5 0,3"`) rather
    /// than a decimal point: exactly two tokens, *both* comma-decimal-shaped. Always a
    /// hard error, even during header-skip — unlike a merely-skippable header line that
    /// happens to contain one comma-formatted number amid prose (only one token would be
    /// comma-decimal-shaped there, which doesn't qualify).
    DecimalComma,
    /// Anything else: header/title/comment text, or (once data has started) garbage.
    Other,
}

/// Tries to interpret `tokens` as exactly two numeric fields. Returns `None` if `tokens`
/// doesn't have exactly two entries (the caller should try a different separator), or if
/// there aren't clearly two intended numbers at all.
fn try_two_tokens(tokens: &[&str]) -> Option<LineKind> {
    if tokens.len() != 2 {
        return None;
    }
    let t0 = tokens[0].trim();
    let t1 = tokens[1].trim();
    match (t0.parse::<f64>(), t1.parse::<f64>()) {
        (Ok(a), Ok(b)) => Some(LineKind::Row(a, b)),
        _ => {
            // Both tokens must look like decimal-comma numbers before we call this an
            // attempted (malformed) data row. If only one token looks comma-decimal-ish
            // (e.g. a header line mixing prose with a comma-formatted number, such as
            // "Density 1,225"), that's not a data row at all — leave it as `None` so the
            // caller falls through to ordinary header/garbage handling instead of a hard
            // decimal-comma error.
            if looks_like_comma_decimal(t0) && looks_like_comma_decimal(t1) {
                Some(LineKind::DecimalComma)
            } else {
                None
            }
        }
    }
}

/// True if `tok` looks like a number written with a decimal comma instead of a decimal
/// point (e.g. `"0,523"`): no `.` present, exactly one `,`, and swapping that `,` for `.`
/// makes it parse as `f64`.
fn looks_like_comma_decimal(tok: &str) -> bool {
    if tok.is_empty() || tok.contains('.') {
        return false;
    }
    if tok.matches(',').count() != 1 {
        return false;
    }
    tok.replace(',', ".").parse::<f64>().is_ok()
}

/// Classifies one line of `.drg` text. Tries whitespace-separated fields first (the real
/// vendor format uses tab), then comma, then semicolon.
fn classify_line(line: &str) -> LineKind {
    let trimmed = line.trim();
    if trimmed.is_empty() {
        return LineKind::Blank;
    }
    let ws_tokens: Vec<&str> = trimmed.split_whitespace().collect();
    if let Some(kind) = try_two_tokens(&ws_tokens) {
        return kind;
    }
    let comma_tokens: Vec<&str> = trimmed.split(',').map(str::trim).collect();
    if let Some(kind) = try_two_tokens(&comma_tokens) {
        return kind;
    }
    let semi_tokens: Vec<&str> = trimmed.split(';').map(str::trim).collect();
    if let Some(kind) = try_two_tokens(&semi_tokens) {
        return kind;
    }
    LineKind::Other
}

/// True if every consecutive pair in `vals` is strictly increasing.
fn is_strictly_ascending(vals: &[f64]) -> bool {
    vals.windows(2).all(|w| w[0] < w[1])
}

/// Parses `.drg` drag-curve text into a [`ParsedDragCurve`]. See the module docs for the
/// full tolerance contract. On any violation, returns `Err(String)` naming the offending
/// line number (1-indexed) and the problem.
pub fn parse_drg(text: &str) -> Result<ParsedDragCurve, String> {
    let mut name: Option<String> = None;
    let mut header_done = false;
    let mut raw: Vec<(usize, f64, f64)> = Vec::new();

    for (idx, line) in text.lines().enumerate() {
        let lineno = idx + 1;
        match classify_line(line) {
            LineKind::DecimalComma => {
                return Err(format!(
                    "line {lineno}: numbers appear to use a decimal comma (e.g. \"0,5\"); \
                     .drg files must use a decimal point"
                ));
            }
            LineKind::Row(a, b) => {
                header_done = true;
                if raw.len() >= MAX_DRG_ROWS {
                    return Err(format!(
                        "line {lineno}: too many data rows (more than {MAX_DRG_ROWS})"
                    ));
                }
                raw.push((lineno, a, b));
            }
            LineKind::Blank => {}
            LineKind::Other => {
                if header_done {
                    return Err(format!(
                        "line {lineno}: expected two numbers, found {:?}",
                        line.trim()
                    ));
                }
                let trimmed = line.trim();
                if name.is_none() && !trimmed.is_empty() {
                    name = Some(trimmed.to_string());
                }
            }
        }
    }

    if raw.len() < MIN_DRG_ROWS {
        return Err(format!(
            "found {} data row(s); need at least {MIN_DRG_ROWS}",
            raw.len()
        ));
    }

    let col0: Vec<f64> = raw.iter().map(|&(_, a, _)| a).collect();
    let col1: Vec<f64> = raw.iter().map(|&(_, _, b)| b).collect();

    // Column detection. Real vendor files store (cd, mach) — column 0 is *not* mach there
    // — so this detection must run for every file, not just as a defensive fallback.
    let col0_ascends = is_strictly_ascending(&col0);
    let col1_ascends = is_strictly_ascending(&col1);
    let mach_is_col0 = match (col0_ascends, col1_ascends) {
        (true, false) => true,
        (false, true) => false,
        (false, false) => {
            return Err(
                "neither column is strictly ascending; expected a mach column".to_string(),
            );
        }
        (true, true) => {
            // Both columns ascend (e.g. a (cd, mach) file whose cd happens to be
            // monotonic across a subsonic-only deck). Break the tie by magnitude: mach
            // spans past 1.0 in real decks, cd stays well under ~1.5, so the column with
            // the larger maximum is mach — unless the two maxima are too close to call.
            let max0 = col0.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let max1 = col1.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let (larger, smaller) = if max0 >= max1 { (max0, max1) } else { (max1, max0) };
            if larger > 0.0 && smaller / larger >= 0.8 {
                return Err(
                    "ambiguous columns: both ascend with similar ranges; cannot determine \
                     which is mach"
                        .to_string(),
                );
            }
            max0 > max1
        }
    };

    let mut points = Vec::with_capacity(raw.len());
    for &(lineno, a, b) in &raw {
        let (mach, cd) = if mach_is_col0 { (a, b) } else { (b, a) };
        if !mach.is_finite() || mach < 0.0 {
            return Err(format!(
                "line {lineno}: mach must be finite and >= 0, got {mach}"
            ));
        }
        if !cd.is_finite() || cd <= 0.0 {
            return Err(format!("line {lineno}: cd must be finite and > 0, got {cd}"));
        }
        points.push((mach, cd));
    }

    Ok(ParsedDragCurve { name, points })
}

/// Cheap sniff: does `text` look like `.drg` (or headerless two-column CSV of the same
/// shape)? Finds the first two data rows and checks that one of the two columns is
/// ascending between them — a bounded scan that does **not** allocate the full points
/// vector (unlike calling [`parse_drg`] and checking `is_ok()`), so it's cheap to use as
/// a pre-check before committing to a full parse.
pub fn looks_like_drg(text: &str) -> bool {
    let mut first: Option<(f64, f64)> = None;
    for line in text.lines() {
        match classify_line(line) {
            LineKind::Row(a, b) => match first {
                None => first = Some((a, b)),
                Some((pa, pb)) => return pa < a || pb < b,
            },
            LineKind::DecimalComma => return false,
            LineKind::Blank | LineKind::Other => {}
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    const SYNTH: &str = "SYN123, Synthetic Test Bullet 7.62mm 10.00g\n\
                         Doppler radar test data (synthetic)\n\
                         0.50 0.230\n\
                         0.80 0.210\n\
                         0.95 0.280\n\
                         1.05 0.450\n\
                         1.50 0.420\n\
                         2.50 0.300\n";

    #[test]
    fn parses_synthetic_drg_with_header() {
        let c = parse_drg(SYNTH).unwrap();
        assert_eq!(c.points.len(), 6);
        assert_eq!(c.points[0], (0.50, 0.230));
        assert_eq!(c.points[5], (2.50, 0.300));
        assert_eq!(c.name.as_deref(), Some("SYN123, Synthetic Test Bullet 7.62mm 10.00g"));
    }

    #[test]
    fn rejects_bad_decks() {
        // descending mach
        assert!(parse_drg("t\n1.0 0.3\n0.9 0.3\n").is_err());
        // single row
        assert!(parse_drg("t\n1.0 0.3\n").is_err());
        // non-positive cd
        assert!(parse_drg("t\n0.5 0.0\n1.0 0.3\n").is_err());
        // garbage after data starts
        assert!(parse_drg("t\n0.5 0.3\nnot numbers\n1.0 0.3\n").is_err());
        // decimal commas -> clear error mentioning decimal
        let e = parse_drg("t\n0,5 0,3\n1,0 0,31\n").unwrap_err();
        assert!(e.to_lowercase().contains("decimal") || e.to_lowercase().contains("comma"), "{e}");
        // empty / header-only
        assert!(parse_drg("").is_err());
        assert!(parse_drg("just a title\n").is_err());
    }

    #[test]
    fn sniff_distinguishes_drg_from_csv_and_junk() {
        assert!(looks_like_drg(SYNTH));
        assert!(!looks_like_drg("mach,cd\n")); // header-only csv
        assert!(!looks_like_drg("hello world\nthis is prose\n"));
        // a plain two-column csv WITHOUT header is also a valid drg shape — sniffing
        // may accept it; that is fine (same numbers either way). Just pin the behavior:
        assert!(looks_like_drg("0.5,0.3\n1.0,0.31\n"));
    }

    #[test]
    fn row_cap_enforced() {
        let mut s = String::from("t\n");
        for i in 0..4097 {
            s.push_str(&format!("{} 0.3\n", 0.1 + i as f64 * 0.001));
        }
        assert!(parse_drg(&s).is_err());
    }

    // --- Additional coverage beyond the brief's baseline, informed by the real vendor
    // sample inspected for this task (structure only; no vendor data reproduced here) ---

    #[test]
    fn accepts_cd_mach_column_order_like_the_real_vendor_format() {
        // Exercises column-order detection: some files store (cd, mach) rather than
        // (mach, cd). Everything here — header text, field values, separators — is
        // invented for this test; it shares no structure with any vendor file.
        let synth_cd_mach =
            "synthetic reversed-column test deck, invented values\r\n\
             0.230\t0.000\r\n\
             0.210\t0.400\r\n\
             0.280\t0.900\r\n\
             0.450\t1.050\r\n\
             0.300\t2.500\r\n";
        let c = parse_drg(synth_cd_mach).unwrap();
        assert_eq!(c.points.len(), 5);
        assert_eq!(c.points[0], (0.000, 0.230));
        assert_eq!(c.points[4], (2.500, 0.300));
    }

    #[test]
    fn semicolon_separated_rows_are_accepted() {
        let c = parse_drg("t\n0.5;0.3\n1.0;0.31\n").unwrap();
        assert_eq!(c.points, vec![(0.5, 0.3), (1.0, 0.31)]);
    }

    #[test]
    fn headerless_csv_has_no_name() {
        let c = parse_drg("0.5,0.3\n1.0,0.31\n").unwrap();
        assert_eq!(c.name, None);
    }

    #[test]
    fn mach_may_start_at_zero() {
        // The real vendor sample's first row is mach = 0 exactly; confirm the >= 0
        // (not > 0) tolerance is intentional and doesn't reject this.
        let c = parse_drg("t\n0.0 0.3\n1.0 0.31\n").unwrap();
        assert_eq!(c.points[0], (0.0, 0.3));
    }

    #[test]
    fn looks_like_drg_is_false_for_empty_and_single_row() {
        assert!(!looks_like_drg(""));
        assert!(!looks_like_drg("t\n1.0 0.3\n"));
    }

    // --- Column tie-break when both columns strictly ascend (MBA-1409 review finding 2) ---

    #[test]
    fn both_ascend_larger_max_col1_rescues_subsonic_only_cd_mach_deck() {
        // A (cd, mach) deck limited to the subsonic range can have a strictly ascending
        // cd column too (col0 here: 0.20 -> 0.30). Column 1 still wins because its
        // maximum (2.0) is far larger than column 0's (0.30) — mach spans past 1.0,
        // cd doesn't.
        let c = parse_drg("t\n0.20 0.5\n0.25 1.0\n0.30 2.0\n").unwrap();
        assert_eq!(c.points, vec![(0.5, 0.20), (1.0, 0.25), (2.0, 0.30)]);
    }

    #[test]
    fn both_ascend_larger_max_col0_is_mach() {
        // Ordinary (mach, cd) order where cd also happens to ascend (0.20 -> 0.30).
        // Column 0 wins because its maximum (2.0) is far larger than column 1's (0.30).
        let c = parse_drg("t\n0.5 0.20\n1.0 0.25\n2.0 0.30\n").unwrap();
        assert_eq!(c.points, vec![(0.5, 0.20), (1.0, 0.25), (2.0, 0.30)]);
    }

    #[test]
    fn both_ascend_similar_maxima_is_ambiguous_error() {
        // Both columns ascend and their maxima (1.0 vs 1.1) are within 20% of each
        // other -- too close to call which one is mach.
        let e = parse_drg("t\n0.5 0.6\n0.8 0.9\n1.0 1.1\n").unwrap_err();
        assert!(e.contains("ambiguous"), "{e}");
    }

    // --- Header-phase decimal-comma false positive (MBA-1409 review finding 3) ---

    #[test]
    fn header_line_with_prose_and_comma_number_is_skipped_not_flagged_as_decimal_comma() {
        // A header/metadata line mixing prose with a comma-formatted number (only one of
        // its two tokens looks comma-decimal-shaped) must be treated as ordinary
        // skippable header text, not misdiagnosed as an attempted decimal-comma data row.
        let c = parse_drg("Density 1,225\n0.5 0.3\n1.0 0.31\n").unwrap();
        assert_eq!(c.points, vec![(0.5, 0.3), (1.0, 0.31)]);
        assert_eq!(c.name.as_deref(), Some("Density 1,225"));
    }
}
