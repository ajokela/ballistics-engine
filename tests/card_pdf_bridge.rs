//! `card.pdf`: the printable PDF dope card through the JSON bridge.
//!
//! The load-bearing property is not "a PDF came back" — it is that **the printed numbers
//! are the numbers the shooter already read on screen**. An app stores one `CardRequestV1`
//! per saved card, plus the `card.range_table` response it displayed; at export time it
//! hands both back, and `card.pdf` prints THOSE rows (`stored_card`) rather than solving
//! again. Re-solving is what let a reprint drift: the rows are a function of the engine build
//! and of the correction-table file at the stored path, and both move under a saved card.
//!
//! So the tests here work from a stored response and read the figures back out of the PDF —
//! including with the correction table the request names deliberately absent, which a
//! re-solve cannot survive and a reprint must not notice. The solve path is still exercised
//! (it is what a CLI-shaped caller gets) and cross-checked against the CLI's own dope card
//! for the same load.
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

    // -- the stored card an app replays for print (E1) -----------------------------------

    /// The `card.range_table` response an app has stored for a saved card, exactly as the
    /// engine emitted it — this is the document `stored_card.card` takes, pasted verbatim.
    ///
    /// The adjustments are hand-picked rather than solved: each one pins a case of the
    /// one-decimal-place rule both surfaces must agree on, ties included (see
    /// `the_printed_adjustments_are_the_stored_rows_at_the_apps_precision`). They are
    /// physically ordered so the card still reads like a card.
    fn stored_range_table() -> Value {
        let row = |range: f64, drop_linear: f64, drop_adj: f64, wind_linear: f64, wind_adj: f64,
                   velocity: f64, energy: f64, time: f64| {
            json!({
                "range": range, "drop_linear": drop_linear, "drop_adj": drop_adj,
                "wind_linear": wind_linear, "wind_adj": wind_adj,
                "velocity": velocity, "energy": energy, "time": time
            })
        };
        json!({
            "schema_version": 1,
            "kind": "range_table",
            "zero_distance": 100.0,
            "bc_for_solve": 0.2381,
            "units": {
                "distance": "yd", "velocity": "fps", "energy": "ft-lb", "drop": "in",
                "wind_speed": "mph",
                "elevation_adjustment": "MIL", "windage_adjustment": "MIL"
            },
            "rows": [
                row(100.0,    0.0, 0.0,     0.0,  0.0,    2447.0, 2327.0, 0.117),
                row(200.0,   -4.4, 0.65,   -1.1, -0.105,  2302.0, 2059.0, 0.244),
                row(300.0,  -15.7, 1.4551, -2.4, -0.21,   2162.0, 1816.0, 0.372),
                row(400.0,  -35.2, 2.4478, -4.6, -0.315,  2026.0, 1595.0, 0.52),
                row(500.0,  -64.3, 3.575,  -7.4, -0.42,   1894.0, 1394.0, 0.68),
                row(600.0, -104.6, 4.8412, -11.0, -0.55,  1766.0, 1212.0, 0.86),
                row(700.0, -158.0, 6.25,  -15.0, -0.65,   1642.0, 1047.0, 1.06),
                row(800.0, -225.0, 7.35,  -20.0, -0.75,   1522.0,  900.0, 1.28)
            ]
        })
    }

    /// `[range, drop, wind]` for every stored row, formatted the way the apps format a
    /// stored row for the screen: whole yards, and **one** decimal place on a MIL/MOA/SMOA/
    /// IPHY adjustment (turret resolution is 0.1 — the apps' old two decimals are what put
    /// `2.45` on screen against `2.4` on paper).
    ///
    /// These are hard-coded strings, not `format!("{:.1}")` of the fixture: a test that
    /// derives the expectation from the same rule it is checking cannot fail. The two
    /// pinned near-ties are the cross-platform hazard — `6.25` is an exact binary tie
    /// (renders `6.2`, half-to-even) and `7.35` is just BELOW the tie (renders `7.3`), so a
    /// client that rounds the shortest decimal string half-up prints `6.3`/`7.4` and
    /// disagrees with the paper by a click.
    fn expected_stored_runs() -> Vec<Vec<String>> {
        [
            ["100", "0.0", "0.0"],
            ["200", "0.7", "-0.1"],
            ["300", "1.5", "-0.2"],
            ["400", "2.4", "-0.3"],
            ["500", "3.6", "-0.4"],
            ["600", "4.8", "-0.6"],
            ["700", "6.2", "-0.7"],
            ["800", "7.3", "-0.8"],
        ]
        .iter()
        .map(|run| run.iter().map(|s| s.to_string()).collect())
        .collect()
    }

    /// The full `stored_card` block: the stored response plus the provenance the footer
    /// prints, so a reprint can be reconciled with the screen afterwards.
    fn stored_card_block() -> Value {
        json!({
            "card": stored_range_table(),
            "engine_version": "0.34.1",
            "bc5d_table_version": "2.5.0"
        })
    }

    /// The fixture request with the stored card attached — an app's export call.
    fn stored_request() -> Value {
        let mut request = fixture_request();
        request["end"] = json!(800.0);
        request["stored_card"] = stored_card_block();
        request
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

    /// (b) The solve path's row-for-row correspondence: `card.range_table` and `card.pdf`,
    /// given the same request in the same build, must put the same figures on the screen and
    /// on the paper.
    ///
    /// NOTE what this test is not. It derives its expectation with `cell_adjustment`, i.e.
    /// with the PDF's own rounding rule, so it cannot detect a PRECISION disagreement — that
    /// is what `the_printed_adjustments_are_the_stored_rows_at_the_apps_precision` is for,
    /// with hard-coded strings. And it holds only within one build for one request, which is
    /// why a saved card is reprinted from `stored_card` instead.
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
        let message = out["error"]["message"].as_str().unwrap();
        assert!(message.contains("the limit is"), "{message}");
        // The refusal describes the document — six rows on one page, made huge by a label —
        // and names no control a saved card does not have.
        assert!(message.contains("6 rows") && message.contains("1 pages"), "{message}");
        for absent in ["coarsen", "shorten"] {
            assert!(!message.contains(absent), "{message}");
        }
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

    // -- printing the stored rows (E1) ----------------------------------------------------

    /// THE falsifiable test for "a reprint is a reprint": the stored request points
    /// `bc5d_table_path` at a table that is not on this machine.
    ///
    /// * Without `stored_card`, `card.pdf` re-solves — so it opens that table and fails.
    /// * With `stored_card`, nothing is solved and no table is opened, so the SAME request
    ///   succeeds and prints the stored rows.
    ///
    /// A shooter who deletes a correction table (or whose app refreshes the table set,
    /// overwriting `bc5d_<caliber>.bin` in place) can still reprint the card in his pocket.
    #[test]
    fn stored_rows_print_without_a_solve_or_a_correction_table() {
        let missing_table = std::env::temp_dir().join("bx_no_such_bc5d_308.bin");
        assert!(!missing_table.exists(), "fixture path must not exist");

        let mut solving = fixture_request();
        solving["end"] = json!(800.0);
        solving["bc5d_table_path"] = json!(missing_table.to_str().unwrap());
        let out = call("card.pdf", solving.clone());
        assert_eq!(
            out["ok"], false,
            "a re-solving export must open the table (and fail when it is gone): {out}"
        );

        let mut reprinting = solving;
        reprinting["stored_card"] = stored_card_block();
        let out = call("card.pdf", reprinting);
        assert_eq!(out["ok"], true, "printing stored rows must not open the table: {out}");
        assert_eq!(out["result"]["source"], "stored_rows", "{out}");
        assert_eq!(out["result"]["row_count"], 8, "{out}");

        let bytes = pdf_bytes(&out["result"]);
        assert_is_a_pdf(&bytes, "reprinted card");
        let printed = tokens(&bytes, "reprinted card");
        for run in expected_stored_runs() {
            assert_run_present(&printed, &run, "reprinted card");
        }
    }

    /// (b) One precision per axis unit, on both surfaces. The stored rows are formatted the
    /// way the apps format them for the screen — one decimal place on an angular
    /// adjustment — and those exact strings must be the ones on the paper.
    ///
    /// This is the test the old parity check could not be: it formatted the on-screen rows
    /// *to the PDF's own precision* before searching for them, so a screen that printed two
    /// decimals still passed. `expected_stored_runs` is a table of literals instead.
    #[test]
    fn the_printed_adjustments_are_the_stored_rows_at_the_apps_precision() {
        let (bytes, pages) = generate(stored_request());
        assert_eq!(pages, 1, "eight rows are one page");
        let printed = tokens(&bytes, "stored card");
        for run in expected_stored_runs() {
            assert_run_present(&printed, &run, "stored card");
        }
        // ...and nothing on the paper carries a second decimal place for an adjustment: the
        // two-decimal spellings of the same rows must be absent.
        for absent in ["2.45", "1.46", "3.57", "-0.31", "6.25", "7.35"] {
            assert!(
                !printed.iter().any(|t| t == absent),
                "{absent} must not appear: the card prints one decimal place, {printed:?}"
            );
        }
    }

    /// The response says which card it printed and where the rows came from, so a client can
    /// verify it got a reprint rather than a re-solve (and that every stored row was drawn).
    #[test]
    fn the_response_states_which_card_it_printed_and_where_the_rows_came_from() {
        let out = call("card.pdf", fixture_request());
        assert_eq!(out["ok"], true, "{out}");
        assert_eq!(out["result"]["kind"], "range_table", "{out}");
        assert_eq!(out["result"]["source"], "solve", "{out}");
        assert_eq!(out["result"]["row_count"], 6, "{out}");

        let out = call("card.pdf", stored_request());
        assert_eq!(out["result"]["kind"], "range_table", "{out}");
        assert_eq!(out["result"]["source"], "stored_rows", "{out}");
        assert_eq!(out["result"]["row_count"], 8, "{out}");
    }

    /// (E2) A stored `card.wind` request must not come back as a range-table PDF with a
    /// Wind column of zeroes. `wind_speeds`/`wind_angles_deg` are the wind card's defining
    /// fields; a range-table PDF ignores them, so the request is refused rather than
    /// silently answered.
    #[test]
    fn a_wind_card_request_is_refused_instead_of_printed_as_a_range_table() {
        for field in ["wind_speeds", "wind_angles_deg"] {
            let mut request = fixture_request();
            request["wind_speed"] = json!(0.0);
            request[field] = json!([5.0, 10.0, 15.0, 20.0]);
            let out = call("card.pdf", request);
            assert_eq!(out["ok"], false, "{field} must be refused, not ignored: {out}");
            // The service's own typed rejections ride as `command_failed`, exactly like the
            // sibling card commands; `invalid_request` is the bridge's own decode failure.
            assert_eq!(out["error"]["code"], "command_failed", "{field}: {out}");
            let message = out["error"]["message"].as_str().unwrap();
            assert!(message.contains(field), "the refusal must name {field}: {message}");
            assert!(
                message.contains("range_table"),
                "the refusal must say which card card.pdf prints: {message}"
            );
        }
    }

    /// (E2) ...and the same for a stored card of a kind this surface cannot print: the
    /// stored response's own `kind` is consulted, so a come-ups card is refused instead of
    /// being reprinted as a range table with a Wind column the shooter never saw.
    #[test]
    fn a_stored_card_of_a_kind_this_surface_cannot_print_is_refused() {
        for kind in ["come_ups", "wind_card", "truing"] {
            let mut request = stored_request();
            request["stored_card"]["card"]["kind"] = json!(kind);
            let out = call("card.pdf", request);
            assert_eq!(out["ok"], false, "a {kind} card must be refused: {out}");
            assert_eq!(out["error"]["code"], "command_failed", "{kind}: {out}");
            assert!(
                out["error"]["message"].as_str().unwrap().contains(kind),
                "the refusal must name the stored kind: {out}"
            );
        }
    }

    /// (E1 fix 3) Provenance on the paper: the footer states the engine version and the
    /// correction-table version the rows came from, so a printed card and a screen can be
    /// reconciled later. A solve prints THIS build's version; a reprint prints the stored
    /// card's. The footer BC is the stored card's own `bc_for_solve`, not the request's
    /// published BC — a BC5D-corrected card's footer states the BC its numbers used.
    #[test]
    fn the_footer_states_the_engine_and_table_versions() {
        let (bytes, _pages) = generate(stored_request());
        let printed = tokens(&bytes, "stored card");
        assert!(
            printed.iter().any(|t| t == "Engine:0.34.1"),
            "the stored engine version must be printed: {printed:?}"
        );
        assert!(
            printed.iter().any(|t| t == "Table:2.5.0"),
            "the stored table version must be printed: {printed:?}"
        );
        assert!(
            printed.iter().any(|t| t == "BC:0.238"),
            "the footer BC must be the stored card's bc_for_solve: {printed:?}"
        );

        let (bytes, _pages) = generate(fixture_request());
        let printed = tokens(&bytes, "solved card");
        let this_build = format!("Engine:{}", env!("CARGO_PKG_VERSION"));
        assert!(
            printed.contains(&this_build),
            "a freshly solved card must print this build's version ({this_build}): {printed:?}"
        );
        assert!(
            !printed.iter().any(|t| t.starts_with("Table:")),
            "no table version is known for a solve, so none is claimed: {printed:?}"
        );
    }

    /// The on-screen card reports the BC its rows were computed with, which is what an app
    /// stores and what a reprint's footer then states. Without it there is no way for a
    /// saved card to print the corrected BC its numbers came from.
    #[test]
    fn the_on_screen_card_reports_the_bc_its_rows_used() {
        let out = call("card.range_table", fixture_request());
        assert_eq!(out["ok"], true, "{out}");
        assert_eq!(
            out["result"]["bc_for_solve"].as_f64(),
            Some(0.243),
            "an uncorrected card reports its published BC: {out}"
        );
    }

    /// A stored response and a request that disagree about what a row MEANS cannot be the
    /// same card: refuse, rather than print MIL numbers under a MOA heading.
    #[test]
    fn stored_rows_and_a_request_that_disagree_on_units_are_refused() {
        let mut request = stored_request();
        request["adjustment_unit"] = json!("moa");
        let out = call("card.pdf", request);
        assert_eq!(out["ok"], false, "{out}");
        assert_eq!(out["error"]["code"], "command_failed", "{out}");
        let message = out["error"]["message"].as_str().unwrap();
        assert!(message.contains("MOA") && message.contains("MIL"), "{message}");

        let mut request = stored_request();
        request["units"] = json!("metric");
        let out = call("card.pdf", request);
        assert_eq!(out["ok"], false, "a metric request cannot own yard rows: {out}");
        assert_eq!(out["error"]["code"], "command_failed", "{out}");
    }

    /// (E4/E5) Too many rows is refused from the row count itself — before a document is
    /// built — and the refusal states what is true of a stored card: how many rows and how
    /// many pages. It names no "coarsen the step" control, because a saved card's domain is
    /// immutable.
    #[test]
    fn too_many_rows_is_refused_with_the_row_and_page_counts() {
        let rows: Vec<Value> = (0..6_000)
            .map(|i| json!({"range": 100.0 + i as f64, "drop_adj": 1.0, "wind_adj": 0.1, "time": 0.5}))
            .collect();
        let mut request = stored_request();
        request["end"] = json!(6_100.0);
        request["stored_card"]["card"]["rows"] = json!(rows);

        let out = call("card.pdf", request);
        assert_eq!(out["ok"], false, "{}", out["result"]["page_count"]);
        assert_eq!(out["error"]["code"], "resource_limit", "{}", out["error"]);
        let message = out["error"]["message"].as_str().unwrap();
        assert!(message.contains("6000 rows"), "the refusal must state the rows: {message}");
        assert!(message.contains("pages"), "the refusal must state the pages: {message}");
        for absent in ["coarsen", "shorten"] {
            assert!(
                !message.contains(absent),
                "the refusal must not name a control a saved card does not have: {message}"
            );
        }
    }

    /// The same refusal on the solve path, which is the one an app can reach by accident: a
    /// free-text `step` of 0.1 over a long domain. The rows are counted the moment they exist
    /// and the document is never built.
    #[test]
    fn a_solved_card_with_too_many_rows_is_refused_before_the_document_is_built() {
        let mut request = fixture_request();
        request["end"] = json!(800.0);
        request["step"] = json!(0.1);
        let out = call("card.pdf", request);
        assert_eq!(out["ok"], false, "{}", out["result"]["page_count"]);
        assert_eq!(out["error"]["code"], "resource_limit", "{}", out["error"]);
        let message = out["error"]["message"].as_str().unwrap();
        assert!(message.contains("7001 rows"), "{message}");
        assert!(message.contains("72 pages"), "{message}");
    }

    /// `stored_card` is a `card.pdf` key, not part of a saved card: it is the card's stored
    /// RESPONSE, attached at export time. An app that stored it inside the request would find
    /// the request no longer replayable, so the on-screen command says so loudly instead of
    /// ignoring it.
    #[test]
    fn stored_card_is_not_a_field_of_the_request_the_screen_replays() {
        let out = call("card.range_table", stored_request());
        assert_eq!(out["ok"], false, "{out}");
        assert_eq!(out["error"]["code"], "invalid_request", "{out}");
        assert!(
            out["error"]["message"].as_str().unwrap().contains("stored_card"),
            "{out}"
        );
    }

    /// The Lead column of a reprint is derived from the STORED time of flight — arithmetic
    /// on the numbers the shooter already has, not a new trajectory. 10 mph across a 400 yd
    /// row whose stored ToF is 0.52 s is 6.4 MIL of hold.
    #[test]
    fn the_lead_column_of_a_reprint_comes_from_the_stored_time_of_flight() {
        let (bytes, _pages) = generate(stored_request());
        let printed = tokens(&bytes, "stored card");
        let anchor = printed
            .iter()
            .position(|t| t == "400")
            .unwrap_or_else(|| panic!("no 400 yd row: {printed:?}"));
        assert_eq!(
            printed[anchor..anchor + 4],
            ["400", "2.4", "-0.3", "6.4"].map(str::to_string),
            "the 400 yd row, Lead included: {printed:?}"
        );
    }

    /// A stored card written by a NEWER engine still reprints its own cells.
    ///
    /// The two apps carry independent engine pins and ship through separate stores, so a
    /// staggered release is the normal case: a card saved on the platform that already
    /// bumped is synced to the one that has not. The document being read back here is the
    /// engine's own output, never user input, so a field this build has not heard of is
    /// something to IGNORE, not to refuse — refusing makes the card unprintable on the
    /// older platform with a serde message listing internal struct fields, which is the one
    /// outcome the reprint path exists to prevent.
    #[test]
    fn a_stored_card_from_a_newer_engine_still_prints_its_stored_rows() {
        // A field added to the response, one added to the units block, and one added to a
        // row: every level of the stored document a future `CardResponseV1` can grow.
        let mut request = stored_request();
        request["stored_card"]["card"]["muzzle_regime"] = json!("supersonic");
        request["stored_card"]["card"]["units"]["lead"] = json!("MIL");
        request["stored_card"]["card"]["rows"][3]["bc_segments_used"] = json!(2);

        let out = call("card.pdf", request);
        assert_eq!(out["ok"], true, "a newer engine's stored card must still print: {out}");
        assert_eq!(out["result"]["source"], "stored_rows", "{out}");
        assert_eq!(out["result"]["row_count"], 8, "{out}");

        // ...and it prints the STORED cells, not a re-solve of the request.
        let bytes = pdf_bytes(&out["result"]);
        let printed = tokens(&bytes, "stored card from a newer engine");
        for run in expected_stored_runs() {
            assert_run_present(&printed, &run, "stored card from a newer engine");
        }
    }

    /// The stored-side "this is not a range table" guard must refuse the same fields the
    /// request-side one does. `wind_angles_deg` on the REQUEST is refused; on the stored
    /// card it used to be accepted, so the two halves of one guard disagreed — and the
    /// stored half is the one that has to hold when a future kind starts emitting it.
    #[test]
    fn a_stored_card_carrying_a_wind_matrix_is_refused_field_by_field() {
        for (field, value) in [
            ("wind_speeds", json!([5.0, 10.0])),
            ("wind_angles_deg", json!([90.0])),
            ("extra_angle_rows", json!([[]])),
        ] {
            let mut request = stored_request();
            request["stored_card"]["card"][field] = value;
            let out = call("card.pdf", request);
            assert_eq!(out["ok"], false, "stored {field} must be refused: {out}");
            assert_eq!(out["error"]["code"], "command_failed", "{field}: {out}");
            let message = out["error"]["message"].as_str().unwrap();
            assert!(message.contains(field), "the refusal must name {field}: {message}");
        }
    }

    /// A card name in a script Liberation Sans does not cover is REPORTED, not silently
    /// dropped.
    ///
    /// Both apps pass the card's name through as `pdf.title` and both accept any non-empty
    /// name, so a Japanese, Arabic, Hebrew, Thai or emoji-bearing name reached the header
    /// and vanished from it — an unidentifiable dope card in the shooter's pocket, with
    /// `ok: true` and no warning anywhere in the result. The characters now come back in
    /// `unprintable_title_chars`, and the header prints a visible substitute where each one
    /// stood rather than closing up over the gap.
    #[test]
    fn a_title_the_font_cannot_print_is_reported_and_substituted() {
        // (title, the characters that must be reported, the header token they become)
        for (title, reported, substituted) in [
            // Japanese: the Latin part used to survive and the rest to vanish.
            ("射撃カード 308", "射撃カード", "?????"),
            // Arabic: EVERY character was unprintable, so the header came out blank.
            ("بطاقة الرمي", "بطاقةلرمي", "?????"),
            // An emoji among Latin text — one dropped glyph, no error, no warning.
            ("Match 🎯 .308", "🎯", "?"),
        ] {
            let mut request = stored_request();
            request["pdf"]["title"] = json!(title);

            let out = call("card.pdf", request);
            assert_eq!(out["ok"], true, "{title} must still print: {out}");
            assert_eq!(
                out["result"]["unprintable_title_chars"],
                json!(reported),
                "{title}: the result must name the characters that could not be printed: {out}"
            );

            let bytes = pdf_bytes(&out["result"]);
            let printed = tokens(&bytes, "unprintable title");
            assert!(
                printed.iter().any(|t| t == substituted),
                "{title}: the header must show a substitute where the dropped glyphs stood: \
                 {printed:?}"
            );
        }
    }

    /// A name the font DOES cover reports nothing — the field is a warning, so it must stay
    /// empty for the Latin, Cyrillic and punctuation names that print correctly today.
    #[test]
    fn a_printable_title_reports_no_unprintable_characters() {
        for title in ["Bridge Card", "Карточка .308", "Ø7.82 – 175gr ™ ±½"] {
            let mut request = stored_request();
            request["pdf"]["title"] = json!(title);
            let out = call("card.pdf", request);
            assert_eq!(out["ok"], true, "{title}: {out}");
            assert_eq!(
                out["result"]["unprintable_title_chars"],
                json!(""),
                "{title} prints correctly, so nothing may be reported: {out}"
            );
        }
    }
}
