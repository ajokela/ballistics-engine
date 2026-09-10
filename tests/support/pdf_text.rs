//! Reading text back out of a rendered dope card, shared by every test that asserts on what
//! a PDF actually says (`card_pdf_bridge`, `card_row_range_invariance`).
//!
//! Not a test target of its own — it lives under `tests/support/` (cargo builds only the
//! top-level `tests/*.rs` as test binaries) and is pulled in with
//! `#[path = "support/pdf_text.rs"] mod pdf_text;`.
//!
//! `printpdf` 0.12 writes uncompressed content streams whose text is hex-encoded GLYPH IDS
//! of the embedded Liberation Sans subset, not ASCII — so a naive `grep` for "4.8" finds
//! nothing (an older comment in `tests/dope_card_units.rs` blames compression; the real
//! reason is the glyph encoding). Two extractors are used:
//!
//! * `pdftotext` (poppler) when it is on PATH — the real thing, via the font's ToUnicode map;
//! * otherwise a glyph scan: pull every `<hex> Tj` operand and map glyph id -> character.
//!   Liberation Sans lays its ASCII glyphs out contiguously, so ONE constant offset decodes
//!   the whole card, and the offset is calibrated from the card's own "Range" column header
//!   rather than hard-coded.
//!
//! `card_pdf_bridge::both_extractors_agree_on_the_same_card` cross-checks the two whenever
//! pdftotext is available, so the fallback is not an untested path on machines that have the
//! tool.
//!
//! This module was a private block inside `card_pdf_bridge.rs` until MBA-1477 needed the same
//! reading in the truncation tests. A printed card's contents are asserted from ONE reader,
//! not from a second copy of it that could drift into agreeing with a different document.

// Each including test file uses the part of this it needs.
#![allow(dead_code)]

/// Every `<hex> Tj` operand in the document, as glyph-id runs, in draw order.
pub fn glyph_runs(pdf: &[u8]) -> Vec<Vec<u32>> {
    fn hex_value(c: u8) -> Option<u32> {
        match c {
            b'0'..=b'9' => Some(u32::from(c - b'0')),
            b'a'..=b'f' => Some(u32::from(c - b'a') + 10),
            b'A'..=b'F' => Some(u32::from(c - b'A') + 10),
            _ => None,
        }
    }
    let mut runs = Vec::new();
    let mut i = 0;
    while i < pdf.len() {
        if pdf[i] != b'<' {
            i += 1;
            continue;
        }
        let mut digits = Vec::new();
        let mut j = i + 1;
        while j < pdf.len() {
            match hex_value(pdf[j]) {
                Some(v) => {
                    digits.push(v);
                    j += 1;
                }
                None => break,
            }
        }
        // Only a well-formed `<....> Tj` show-text operand counts; PDF dictionaries
        // (`<<`) and hex strings used for anything else are skipped.
        let closed = j < pdf.len() && pdf[j] == b'>' && !digits.is_empty() && digits.len() % 4 == 0;
        let mut k = j + 1;
        while k < pdf.len() && pdf[k].is_ascii_whitespace() {
            k += 1;
        }
        if closed && pdf[k..].starts_with(b"Tj") {
            runs.push(digits.chunks(4).map(|q| q.iter().fold(0, |acc, d| acc * 16 + d)).collect());
            i = k + 2;
        } else {
            i += 1;
        }
    }
    runs
}

pub fn decode_runs(runs: &[Vec<u32>], offset: u32) -> Vec<String> {
    runs.iter()
        .map(|run| run.iter().filter_map(|&g| char::from_u32(g + offset)).collect())
        .collect()
}

/// Fallback extraction: decode the glyph ids with the one constant offset that makes the
/// card's own "Range" column header appear. Self-calibrating, so a different font subset
/// ordering fails loudly here instead of silently producing garbage that no assertion
/// happens to notice.
pub fn glyph_scan(pdf: &[u8], label: &str) -> String {
    let runs = glyph_runs(pdf);
    assert!(!runs.is_empty(), "{label}: no show-text operands found in the PDF");
    let offset = (0u32..=0x2000)
        .find(|&offset| decode_runs(&runs, offset).iter().any(|s| s == "Range"))
        .unwrap_or_else(|| {
            panic!("{label}: could not calibrate the glyph offset against the \"Range\" header")
        });
    decode_runs(&runs, offset).join("\n")
}

/// Real extraction via poppler, or `None` when `pdftotext` is not installed. Feeds the
/// document on stdin (`pdftotext - -`) so no temporary file is involved.
pub fn pdftotext(pdf: &[u8]) -> Option<String> {
    use std::io::Write;
    use std::process::{Command, Stdio};

    let mut child = Command::new("pdftotext")
        .args(["-layout", "-", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .ok()?;
    child.stdin.take()?.write_all(pdf).ok()?;
    let out = child.wait_with_output().ok()?;
    if !out.status.success() {
        return None;
    }
    Some(String::from_utf8_lossy(&out.stdout).into_owned())
}

/// The card's drawn text as whitespace-separated tokens. Both the two-column table's
/// visual line order (`pdftotext -layout`) and the generator's draw order (glyph scan)
/// emit a row as `range drop wind lead`, so a row's cells are a contiguous token run
/// under either extractor.
pub fn tokens(pdf: &[u8], label: &str) -> Vec<String> {
    let text = pdftotext(pdf).unwrap_or_else(|| glyph_scan(pdf, label));
    text.split_whitespace().map(str::to_string).collect()
}


/// True when `phrase`'s whitespace-separated words appear consecutively in the card's text.
///
/// Asserted on TOKENS rather than the raw string because the two extractors disagree about
/// runs of spaces — `pdftotext -layout` reconstructs inter-word gaps from glyph positions,
/// while the glyph scan emits exactly what was drawn — and a footer line is one drawn string
/// under both.
pub fn contains_phrase(pdf: &[u8], label: &str, phrase: &str) -> bool {
    let want: Vec<String> = phrase.split_whitespace().map(str::to_string).collect();
    let have = tokens(pdf, label);
    !want.is_empty() && have.windows(want.len()).any(|w| w == want.as_slice())
}
