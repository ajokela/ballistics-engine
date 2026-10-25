//! Read text from generated PDFs for presentation assertions.
//!
//! Prefer poppler when installed; otherwise use printpdf's parser, which decodes
//! compressed content streams and embedded font mappings. The integration suite
//! compares both extractors on the same document when poppler is available.

#![allow(dead_code)]

/// Decode the document without requiring an external PDF utility.
pub fn native_text(pdf: &[u8], label: &str) -> String {
    let document =
        printpdf::PdfDocument::parse(pdf, &printpdf::PdfParseOptions::default(), &mut Vec::new())
            .unwrap_or_else(|error| panic!("{label}: could not parse PDF: {error}"));
    let text = document
        .extract_text()
        .into_iter()
        .flatten()
        .collect::<Vec<_>>()
        .join("\n");
    assert!(
        !text.trim().is_empty(),
        "{label}: PDF contains no readable text"
    );
    text
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
/// visual line order (`pdftotext -layout`) and the generator's draw order (native parser)
/// emit a row as `range drop wind lead`, so a row's cells are a contiguous token run
/// under either extractor.
pub fn tokens(pdf: &[u8], label: &str) -> Vec<String> {
    let text = pdftotext(pdf).unwrap_or_else(|| native_text(pdf, label));
    text.split_whitespace().map(str::to_string).collect()
}

/// True when `phrase`'s whitespace-separated words appear consecutively in the card's text.
///
/// Asserted on TOKENS rather than the raw string because the two extractors disagree about
/// runs of spaces — `pdftotext -layout` reconstructs inter-word gaps from glyph positions,
/// while the native parser emits exactly what was drawn — and a footer line is one drawn string
/// under both.
pub fn contains_phrase(pdf: &[u8], label: &str, phrase: &str) -> bool {
    let want: Vec<String> = phrase.split_whitespace().map(str::to_string).collect();
    let have = tokens(pdf, label);
    !want.is_empty() && have.windows(want.len()).any(|w| w == want.as_slice())
}
