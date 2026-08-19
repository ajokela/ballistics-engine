//! `card.pdf`: the printable PDF dope card through the JSON bridge.
//!
//! The load-bearing property is not "a PDF came back" — it is that **the printed numbers
//! are the numbers the shooter already read on screen**. An app stores one
//! `CardRequestV1` per saved card and replays it against `card.range_table` for the display
//! and `card.pdf` for the printout; if those two ever disagree, the card is lying to
//! somebody. So the tests here generate both from the SAME request and read the figures
//! back out of the PDF, and then cross-check the whole thing against the CLI's own dope
//! card for the same load.
//!
//! ## Reading text back out of a dope card
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
//! `both_extractors_agree_on_the_same_card` cross-checks the two whenever pdftotext is
//! available, so the fallback is not an untested path on machines that have the tool.

#![cfg(feature = "bridge")]

use serde_json::{json, Value};

/// One bridge exchange, returning the whole envelope (callers assert `ok` themselves).
fn call(command: &str, request: Value) -> Value {
    let envelope = json!({"api_version": 1, "command": command, "request": request});
    let raw = ballistics_engine::bridge::bridge_call(&envelope.to_string());
    serde_json::from_str(&raw).expect("bridge output must be valid JSON")
}

/// The shared fixture load: a 175gr .308 at 2600 fps, G7 0.243, 100 yd zero, full-value
/// 10 mph wind from the right, 100-600 yd in 100s, MIL on both axes, 10 mph crossing
/// target. `fixture_cli_args` mirrors it flag for flag.
fn fixture_request() -> Value {
    json!({
        "units": "imperial",
        "muzzle_velocity": 2600.0,
        "ballistic_coefficient": 0.243,
        "mass": 175.0,
        "diameter": 0.308,
        "drag_model": "g7",
        "sight_height": 1.5,
        "zero_distance": 100.0,
        "wind_speed": 10.0,
        "wind_direction_deg": 90.0,
        "start": 100.0,
        "end": 600.0,
        "step": 100.0,
        "adjustment_unit": "mil",
        "pdf": {
            "title": "Bridge Card",
            "location": "Test Range",
            "target_speed": 10.0
        }
    })
}

#[cfg(not(feature = "pdf"))]
mod pdf_absent {
    use super::*;

    /// A build without the `pdf` feature must not pretend: `card.pdf` is an unknown
    /// command, and the capability list says so before an app ever tries.
    #[test]
    fn card_pdf_is_unknown_and_unlisted() {
        let out = call("card.pdf", fixture_request());
        assert_eq!(out["ok"], false, "{out}");
        assert_eq!(out["error"]["code"], "unknown_command", "{out}");

        let caps = call("meta.capabilities", Value::Null);
        let commands: Vec<String> =
            serde_json::from_value(caps["result"]["commands"].clone()).unwrap();
        assert!(!commands.contains(&"card.pdf".to_string()), "{caps}");
    }

    /// ...but the request that carries a `pdf` presentation block still parses, so the one
    /// stored request an app keeps per card is not rejected by a pdf-less engine.
    #[test]
    fn the_stored_request_still_drives_the_on_screen_card() {
        let out = call("card.range_table", fixture_request());
        assert_eq!(out["ok"], true, "{out}");
        assert!(!out["result"]["rows"].as_array().unwrap().is_empty(), "{out}");
    }
}

#[cfg(feature = "pdf")]
mod pdf_present {
    use super::*;

    // -- transport helpers ---------------------------------------------------------------

    /// Decode the response's `pdf_base64`, checking it against the reported `byte_length`.
    fn pdf_bytes(result: &Value) -> Vec<u8> {
        const ALPHABET: &[u8; 64] =
            b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        let text = result["pdf_base64"].as_str().expect("pdf_base64 must be a string");
        let mut acc: u32 = 0;
        let mut bits: u32 = 0;
        let mut out = Vec::with_capacity(text.len() / 4 * 3);
        for c in text.bytes() {
            if c == b'=' {
                break;
            }
            let sextet = ALPHABET
                .iter()
                .position(|&a| a == c)
                .unwrap_or_else(|| panic!("pdf_base64 has a non-alphabet byte {c:?}"))
                as u32;
            acc = (acc << 6) | sextet;
            bits += 6;
            if bits >= 8 {
                bits -= 8;
                out.push((acc >> bits) as u8);
            }
        }
        assert_eq!(
            result["byte_length"].as_u64().expect("byte_length must be a number"),
            out.len() as u64,
            "byte_length must describe the DECODED document, not the base64 text"
        );
        out
    }

    /// Generate a card.pdf and return `(decoded PDF, reported page_count)`.
    fn generate(request: Value) -> (Vec<u8>, usize) {
        let out = call("card.pdf", request);
        assert_eq!(out["ok"], true, "card.pdf failed: {out}");
        let page_count = out["result"]["page_count"].as_u64().expect("page_count") as usize;
        (pdf_bytes(&out["result"]), page_count)
    }

    fn assert_is_a_pdf(bytes: &[u8], label: &str) {
        assert!(bytes.len() > 10_000, "{label}: PDF suspiciously small ({} bytes)", bytes.len());
        assert_eq!(&bytes[..5], b"%PDF-", "{label}: does not start with a PDF header");
    }

    /// Count non-overlapping-enough occurrences of `needle` (windows is fine here: the
    /// needles are PDF dictionary keys, which cannot overlap themselves).
    fn count_bytes(haystack: &[u8], needle: &[u8]) -> usize {
        haystack.windows(needle.len()).filter(|w| *w == needle).count()
    }

    /// Page objects in the document. `/Type/Pages` (the page TREE, one per document) shares
    /// the `/Type/Page` prefix, so it is subtracted out.
    fn page_object_count(bytes: &[u8]) -> usize {
        count_bytes(bytes, b"/Type/Page") - count_bytes(bytes, b"/Type/Pages")
    }

    // -- text extraction ----------------------------------------------------------------

    /// Every `<hex> Tj` operand in the document, as glyph-id runs, in draw order.
    fn glyph_runs(pdf: &[u8]) -> Vec<Vec<u32>> {
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

    fn decode_runs(runs: &[Vec<u32>], offset: u32) -> Vec<String> {
        runs.iter()
            .map(|run| run.iter().filter_map(|&g| char::from_u32(g + offset)).collect())
            .collect()
    }

    /// Fallback extraction: decode the glyph ids with the one constant offset that makes the
    /// card's own "Range" column header appear. Self-calibrating, so a different font subset
    /// ordering fails loudly here instead of silently producing garbage that no assertion
    /// happens to notice.
    fn glyph_scan(pdf: &[u8], label: &str) -> String {
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
    fn pdftotext(pdf: &[u8]) -> Option<String> {
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
    fn tokens(pdf: &[u8], label: &str) -> Vec<String> {
        let text = pdftotext(pdf).unwrap_or_else(|| glyph_scan(pdf, label));
        text.split_whitespace().map(str::to_string).collect()
    }

    fn assert_run_present(tokens: &[String], run: &[String], label: &str) {
        assert!(
            tokens.windows(run.len()).any(|w| w == run),
            "{label}: {run:?} does not appear as a row in the card's text {tokens:?}"
        );
    }

    // -- expected cell formatting -------------------------------------------------------
    //
    // Mirrors `pdf_dope_card::format_range` / `format_adjustment` for the MIL fixture (both
    // are private to that module). The clicks variant — integers, no decimal point — is
    // covered by that module's own unit tests, not here.

    fn cell_range(range: f64) -> String {
        format!("{range:.0}")
    }
    fn cell_adjustment(value: f64) -> String {
        format!("{value:.1}")
    }

    /// `[range, drop, wind]` per on-screen row, in the card's own display formatting.
    fn expected_runs(range_table: &Value) -> Vec<Vec<String>> {
        range_table["rows"]
            .as_array()
            .expect("range_table rows")
            .iter()
            .map(|row| {
                vec![
                    cell_range(row["range"].as_f64().expect("range")),
                    cell_adjustment(row["drop_adj"].as_f64().expect("drop_adj")),
                    cell_adjustment(row["wind_adj"].as_f64().expect("wind_adj")),
                ]
            })
            .collect()
    }

    // -- tests ---------------------------------------------------------------------------

    /// (a) The response really carries a PDF, and its self-description is consistent:
    /// `byte_length` is the decoded length (checked inside `pdf_bytes`) and `page_count`
    /// matches the page objects the document actually contains.
    #[test]
    fn a_generated_card_is_a_non_trivial_pdf_that_describes_itself_correctly() {
        let (bytes, page_count) = generate(fixture_request());
        assert_is_a_pdf(&bytes, "fixture card");
        assert_eq!(page_count, 1, "a six-row card is one page");
        assert_eq!(
            page_object_count(&bytes),
            page_count,
            "the reported page_count must match the document's page objects"
        );
    }

    /// (b) THE load-bearing test: one stored request, replayed against both surfaces, must
    /// put the same figures on the screen and on the paper. Every row's
    /// `range / drop / wind` triple from `card.range_table` has to appear in the PDF's text.
    #[test]
    fn the_printed_figures_are_the_on_screen_figures_for_the_same_request() {
        let request = fixture_request();

        let table = call("card.range_table", request.clone());
        assert_eq!(table["ok"], true, "card.range_table failed: {table}");
        let runs = expected_runs(&table["result"]);
        assert!(
            runs.len() >= 3,
            "fixture must produce at least three rows to be worth comparing, got {}",
            runs.len()
        );

        let (bytes, _pages) = generate(request);
        let printed = tokens(&bytes, "bridge card");
        for run in &runs {
            assert_run_present(&printed, run, "bridge card");
        }
    }

    /// (c) Golden cross-check against the CLI's own dope-card path (`trajectory -o pdf`) for
    /// the same load. Structural equivalence only — same page count, same printed cells —
    /// never byte equality, because three things differ by construction and none of them is
    /// a number a shooter dials:
    ///
    /// * the footer carries a generation timestamp, so two runs never match byte-for-byte;
    /// * the header/footer LABELS differ (the CLI defaults to "Rifle"/"Field"; the bridge
    ///   card is labelled from its own request), and those strings feed the document title,
    ///   hence object lengths and the xref offsets after them;
    /// * `printpdf` embeds font subsets and object offsets that shift with all of the above.
    ///
    /// `--ignore-ground-impact` is required on the CLI side and only there: the CLI's
    /// trajectory truncates at ground impact from the default 5 ft bore height (this load
    /// lands around 537 yd), while the card services sample a pure LOS-referenced flight out
    /// to the requested domain. Without it the CLI card simply stops at 400 yd.
    #[cfg(feature = "cli")]
    #[test]
    fn the_card_matches_the_cli_dope_card_for_the_same_load() {
        use std::process::Command;

        let out_path = std::env::temp_dir().join(format!("bx_card_pdf_golden_{}.pdf", std::process::id()));
        let status = Command::new(env!("CARGO_BIN_EXE_ballistics"))
            .args([
                "trajectory",
                "-v", "2600", "-b", "0.243", "-m", "175", "-d", "0.308",
                "--drag-model", "g7",
                "--sight-height", "1.5",
                "--auto-zero", "100",
                "--max-range", "600",
                "--ignore-ground-impact",
                "--wind-speed", "10", "--wind-direction", "90",
                "--adjustment-unit", "mil",
                "--target-speed", "10",
                // Metres, always — `--sample-interval` does not follow `--units`. 91.44 m is
                // exactly 100 yd, which puts the CLI's samples on the bridge card's 100 yd
                // grid so the two cards describe the same ranges.
                "--sample-trajectory", "--sample-interval", "91.44",
                "-o", "pdf",
                "--output-file", out_path.to_str().unwrap(),
            ])
            .output()
            .expect("run the ballistics CLI");
        assert!(
            status.status.success(),
            "CLI dope card failed: {}",
            String::from_utf8_lossy(&status.stderr)
        );
        let cli_bytes = std::fs::read(&out_path).expect("CLI wrote a PDF");
        let _ = std::fs::remove_file(&out_path);
        assert_is_a_pdf(&cli_bytes, "CLI card");

        let (bridge_bytes, bridge_pages) = generate(fixture_request());
        assert_eq!(
            page_object_count(&cli_bytes),
            bridge_pages,
            "the two cards must paginate the same"
        );

        let cli_tokens = tokens(&cli_bytes, "CLI card");
        let bridge_tokens = tokens(&bridge_bytes, "bridge card");
        // The full four-cell row, Lead column included: `--target-speed 10` on the CLI and
        // `pdf.target_speed: 10.0` on the bridge must land on the same hold.
        let mut compared = 0;
        for range in [100.0_f64, 200.0, 300.0, 400.0, 500.0, 600.0] {
            let anchor = cell_range(range);
            let start = bridge_tokens
                .iter()
                .position(|t| *t == anchor)
                .unwrap_or_else(|| panic!("bridge card has no {anchor} yd row: {bridge_tokens:?}"));
            assert!(
                start + 4 <= bridge_tokens.len(),
                "the {anchor} yd row is truncated at the end of the card's text: {bridge_tokens:?}"
            );
            let row: Vec<String> = bridge_tokens[start..start + 4].to_vec();
            assert_run_present(&cli_tokens, &row, "CLI card");
            compared += 1;
        }
        assert_eq!(compared, 6, "every fixture range must have been compared");
    }

    /// The fallback extractor must agree with poppler where poppler is available, so the
    /// glyph scan is trustworthy on machines (CI included) that lack `pdftotext`.
    ///
    /// Compared as a token MULTISET plus the data rows, not as an ordered sequence: the two
    /// disagree on the column-header block alone, and legitimately. `pdftotext -layout`
    /// reads it visually (eight headings, then eight unit sub-headings), while the glyph scan
    /// reads draw order (heading then sub-heading, per column). Every data row, and every
    /// header/footer string, is identical under both.
    #[test]
    fn both_extractors_agree_on_the_same_card() {
        let (bytes, _pages) = generate(fixture_request());
        let Some(real) = pdftotext(&bytes) else {
            eprintln!("pdftotext not installed; the glyph-scan fallback is what ran everywhere else");
            return;
        };
        let mut real: Vec<String> = real.split_whitespace().map(str::to_string).collect();
        let mut scanned: Vec<String> = glyph_scan(&bytes, "bridge card")
            .split_whitespace()
            .map(str::to_string)
            .collect();

        let table = call("card.range_table", fixture_request());
        for run in expected_runs(&table["result"]) {
            assert_run_present(&real, &run, "pdftotext");
            assert_run_present(&scanned, &run, "glyph scan");
        }

        real.sort();
        scanned.sort();
        assert_eq!(real, scanned, "the two extractors must read the same card");
    }

    /// The card's Range column follows the request's own distance unit — unlike
    /// `trajectory -o pdf`, whose dope card is always yards — and the imperial header block
    /// is converted for display, not re-solved. Both are properties of a stored metric
    /// request replayed for print.
    #[test]
    fn a_metric_request_prints_metres_and_an_imperial_header() {
        // The metric twin of the fixture: 792.48 m/s = 2600 fps, 11.34 g = 175 gr,
        // 7.82 mm = 0.3079 in, 4.4704 m/s = 10 mph.
        let (bytes, _pages) = generate(json!({
            "units": "metric",
            "muzzle_velocity": 792.48,
            "ballistic_coefficient": 0.243,
            "mass": 11.34,
            "diameter": 7.82,
            "drag_model": "g7",
            "sight_height": 38.0,
            "zero_distance": 100.0,
            "wind_speed": 4.4704,
            "wind_direction_deg": 90.0,
            "start": 100.0,
            "end": 500.0,
            "step": 100.0,
            "adjustment_unit": "mil",
            "pdf": {"title": "Metric Card", "target_speed": 4.4704}
        }));
        let printed = tokens(&bytes, "metric card");
        assert!(printed.iter().any(|t| t == "M"), "Range sub-header must be M: {printed:?}");
        assert!(!printed.iter().any(|t| t == "Yd"), "no yards on a metric card: {printed:?}");
        // Header/footer stay imperial on both CLI PDF call sites; the bridge matches.
        assert!(
            printed.iter().any(|t| t == "Weight:175gr"),
            "grams must print as grains in the footer: {printed:?}"
        );
        assert!(
            printed.iter().any(|t| t == "Vel:2600fps"),
            "m/s must print as fps in the footer: {printed:?}"
        );
    }

    /// `font_preset` is a presentation option that has to actually reach the layout: 100 rows
    /// fit on one page at `small` (104 rows/page) and need two at `medium` (98).
    #[test]
    fn the_font_preset_changes_the_pagination() {
        let request = |preset: &str| {
            let mut request = fixture_request();
            request["end"] = json!(1090.0);
            request["step"] = json!(10.0);
            request["pdf"] = json!({"font_preset": preset, "target_speed": 10.0});
            request
        };
        let (_small_bytes, small_pages) = generate(request("small"));
        let (_medium_bytes, medium_pages) = generate(request("medium"));
        assert_eq!(small_pages, 1, "100 rows fit one page at small");
        assert_eq!(medium_pages, 2, "100 rows need two pages at medium");
    }

    /// Absent `target_speed` means "no lead data", which the generator renders as em-dashes
    /// rather than a column of confident-looking zeroes. Only the successful render is
    /// asserted here: the em-dash is outside the contiguous ASCII glyph run the fallback
    /// extractor decodes, and the exact string is pinned by `pdf_dope_card`'s own
    /// `format_adjustment_cell` unit test.
    #[test]
    fn a_card_without_a_target_speed_still_renders() {
        let mut request = fixture_request();
        request["pdf"] = json!({});
        let (bytes, pages) = generate(request.clone());
        assert_is_a_pdf(&bytes, "no-lead card");
        assert_eq!(pages, 1);

        // No `pdf` block at all is the same story, and is what an app that never set a
        // presentation option stores.
        request.as_object_mut().unwrap().remove("pdf");
        let (bytes, _pages) = generate(request);
        assert_is_a_pdf(&bytes, "no-pdf-block card");
    }

    /// (d) The output cap, end to end. The `pdf` block's labels are drawn verbatim, so a
    /// several-hundred-KiB label is the cheapest way to push a document past
    /// `MAX_PDF_BYTES` — and it must come back as a typed `resource_limit` refusal rather
    /// than a multi-megabyte payload. (The boundary itself, at exactly the cap and one byte
    /// over, is unit-tested in `bridge::tests::pdf_output_cap_refuses_only_over_the_limit`.)
    #[test]
    fn an_oversize_card_is_refused_rather_than_returned() {
        let mut request = fixture_request();
        // Each character costs 4 bytes of hex-encoded glyph ids in the content stream, on
        // top of the ~815 KiB of embedded fonts every card carries.
        request["pdf"] = json!({"powder": "P".repeat(900_000)});
        let out = call("card.pdf", request);
        assert_eq!(out["ok"], false, "an over-cap card must not be returned: {}", out["ok"]);
        assert_eq!(out["error"]["code"], "resource_limit", "{}", out["error"]);
        assert!(
            out["error"]["message"].as_str().unwrap().contains("the limit is"),
            "{}",
            out["error"]
        );
    }

    /// Presentation options are validated rather than silently coerced: a stored request
    /// that asks for an impossible font size gets told so, instead of quietly rendering at
    /// some other size.
    #[test]
    fn contradictory_or_out_of_band_presentation_options_are_refused() {
        for (label, pdf_block) in [
            ("both scale and preset", json!({"font_scale": 1.2, "font_preset": "large"})),
            ("scale above the band", json!({"font_scale": 9.0})),
            ("scale below the band", json!({"font_scale": 0.1})),
            ("non-finite scale", json!({"font_scale": f64::MAX})),
            ("unknown preset", json!({"font_preset": "gigantic"})),
        ] {
            let mut request = fixture_request();
            request["pdf"] = pdf_block;
            let out = call("card.pdf", request);
            assert_eq!(out["ok"], false, "{label} must be refused: {out}");
            assert_eq!(out["error"]["code"], "command_failed", "{label}: {out}");
        }
    }

    /// An unknown key inside the `pdf` block is rejected, like every other request payload
    /// on this bridge — a misspelled presentation option should not silently do nothing.
    #[test]
    fn an_unknown_presentation_option_is_an_invalid_request() {
        let mut request = fixture_request();
        request["pdf"] = json!({"titel": "typo"});
        let out = call("card.pdf", request);
        assert_eq!(out["ok"], false, "{out}");
        assert_eq!(out["error"]["code"], "invalid_request", "{out}");
    }
}
