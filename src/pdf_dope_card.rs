//! PDF Dope Card Generation Module
//!
//! Generates printable dope cards in Glenn's proven field-ready format.
//! Format: Two-column layout with Range (yd or m) and Drop/Wind/Lead in MIL or MOA
//! (selected via DopeCardConfig::unit_label; values pre-converted by the caller)
//! Color coding: Black=Range, Red=Drop, Green=Wind, Blue=Lead
//! Row striping for improved readability.
//!
//! 0.33.0 decision-support Plan B Task 10: moved from a `ballistics`-binary-private module
//! (`mod pdf_dope_card;` in `main.rs`) into this library behind the `pdf` feature, and
//! rewritten to consume `&[crate::card::CardRow]` instead of this module's own
//! `DopeCardRow { range_yd: u32, .. }` -- the u32-yards field couldn't express the
//! non-integer ranges Task 11's adaptive card engine (`crate::card`) produces. The caller
//! now states the row range's unit explicitly via `RangeUnit`; see its doc comment.

use crate::card::CardRow;
use printpdf::*;

// Embed Liberation Sans fonts directly into the binary (SIL Open Font License)
static FONT_REGULAR: &[u8] = include_bytes!("../fonts/LiberationSans-Regular.ttf");
static FONT_BOLD: &[u8] = include_bytes!("../fonts/LiberationSans-Bold.ttf");

/// Configuration for the dope card PDF
#[derive(Debug, Clone)]
pub struct DopeCardConfig {
    pub rifle_name: String,
    pub location: String,
    pub density_altitude_ft: f64,
    pub pressure_inhg: f64,
    pub pressure_hpa: f64,
    pub temperature_f: f64,
    pub altitude_ft: f64,
    pub wind_speed_mph: f64,
    pub target_speed_mph: f64,
    pub solver_mode: String,
    pub powder: String,
    pub bullet: String,
    pub weight_gr: f64,
    pub bc: f64,
    pub drag_model: String,
    pub velocity_fps: f64,
    pub font_scale: f32,
    pub bold_data: bool,
    /// Angular unit label shown in the Drop column sub-header ("MIL", "MOA", "SMOA",
    /// "IPHY", or "CLICKS"). `CardRow::drop_adj` is already expressed in this unit
    /// (MBA-1410: independent from `windage_unit_label` below).
    pub elevation_unit_label: String,
    /// Angular unit label shown in the Wind/Lead column sub-headers. `CardRow::
    /// wind_adj`/`lead_adj` are already expressed in this unit -- may differ from
    /// `elevation_unit_label` (MBA-1410 independent elevation/windage unit selection).
    pub windage_unit_label: String,
    /// Engine version that produced these rows, printed in the footer as `Engine:<v>`.
    ///
    /// A card in a shooter's pocket is otherwise impossible to reconcile with a screen: the
    /// rows are a function of the engine build and of the correction table, and both move.
    /// An EMPTY string prints nothing at all -- the same rule the apps' provenance line
    /// follows, because a placeholder ("unknown") on a printed card is worse than silence.
    pub engine_version: String,
    /// Correction-table version these rows were solved against, printed as `Table:<v>`.
    /// Empty prints nothing, which is the honest rendering of "no correction table".
    pub table_version: String,
}

/// Preset font size profiles for dope cards
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FontSizePreset {
    Small,
    Medium,
    Large,
}

impl FontSizePreset {
    pub fn scale(&self) -> f32 {
        match self {
            Self::Small => 0.8,
            Self::Medium => 1.0,
            Self::Large => 1.4,
        }
    }

    // Promoting this module to `pub` in the library (Task 10) makes this method part of
    // the crate's exported API for the first time, which is why clippy only starts
    // flagging it here: it returns `Option<Self>`, not `Result<Self, Self::Err>`, so it
    // doesn't fit std::str::FromStr's contract, and FontSizePreset's shape is preserved
    // verbatim rather than reworked to satisfy the lint.
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "small" | "s" => Some(Self::Small),
            "medium" | "m" => Some(Self::Medium),
            "large" | "l" => Some(Self::Large),
            _ => None,
        }
    }
}

/// Which unit `CardRow::range` is expressed in for this card. Selects only the Range
/// column's sub-header text ("Yd" / "M") -- row values are never converted here, same
/// "the caller already converted it" convention `CardRow` itself documents.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RangeUnit {
    Yards,
    Meters,
}

impl RangeUnit {
    fn label(&self) -> &'static str {
        match self {
            RangeUnit::Yards => "Yd",
            RangeUnit::Meters => "M",
        }
    }
}

/// The `font_scale` band [`generate_dope_card_pdf`] honours; anything outside it is
/// clamped into it (see [`dope_card_rows_per_page`], which clamps identically so a page
/// count and the document it describes cannot disagree).
pub const FONT_SCALE_RANGE: std::ops::RangeInclusive<f32> = 0.5..=3.0;

// Page dimensions (Letter size in mm)
const PAGE_WIDTH: f32 = 215.9;
const PAGE_HEIGHT: f32 = 279.4;
const MARGIN: f32 = 10.0;

// Font sizes
const HEADER_FONT_SIZE: f32 = 9.0;
const TABLE_FONT_SIZE: f32 = 8.0;
const FOOTER_FONT_SIZE: f32 = 8.0;

// Table layout
const ROW_HEIGHT: f32 = 4.5;
const COL_WIDTH: f32 = 24.0; // Width per column (8 columns total)

// Colors (RGB 0.0-1.0)
const COLOR_BLACK: (f32, f32, f32) = (0.0, 0.0, 0.0);
const COLOR_RED: (f32, f32, f32) = (0.78, 0.0, 0.0);
const COLOR_GREEN: (f32, f32, f32) = (0.0, 0.5, 0.0);
const COLOR_BLUE: (f32, f32, f32) = (0.0, 0.0, 0.78);
const COLOR_STRIPE: (f32, f32, f32) = (0.94, 0.94, 0.94); // Light gray for alternating rows
const INHG_TO_HPA: f64 = 33.863_886_666_667;

// The angular conversion (drop_yd/range_yd -> MIL or MOA) and moving-target lead now
// live in main.rs::drop_to_adjustment, so both card units share one code path; the
// caller fills CardRow's drop_adj/wind_adj/lead_adj (as Some(..)) already in the chosen
// unit. Fix-round I-1: a None on any of those three renders as an em-dash (see
// `format_adjustment_cell`), never a fake 0.0 -- Task 11's adaptive engine, the stated
// reason this module takes CardRow at all, emits `lead_adj: None` on every row it
// produces, and a plausible-looking dialed zero would be dangerously wrong there.

/// Calculate density altitude from environmental conditions
///
/// MBA-643: Fixed to interpret pressure as STATION PRESSURE (actual local pressure),
/// not altimeter setting (sea-level corrected). This matches how weather stations
/// and most ballistic tools report pressure.
///
/// Pressure altitude follows the published NWS station-pressure equation:
/// PA = 145366.45 * (1 - (P_hPa/1013.25)^0.190284)
/// <https://www.weather.gov/media/epz/wxcalc/pressureAltitude.pdf>
///
/// ```text
/// DA = PA + 66.7 * (OAT_F - ISA_temp_F)
/// ```
pub fn calculate_density_altitude(_altitude_ft: f64, pressure_inhg: f64, temp_f: f64) -> f64 {
    // The NWS equation is defined in hPa (equivalently millibars), so convert before
    // applying its matched coefficient, reference pressure, and exponent.
    let pressure_hpa = pressure_inhg * INHG_TO_HPA;
    let pressure_alt = 145_366.45 * (1.0 - (pressure_hpa / 1013.25).powf(0.190_284));

    // ISA temperature at pressure altitude (lapse rate: 3.57°F per 1000 ft)
    let isa_temp_f = 59.0 - (pressure_alt / 1000.0) * 3.57;

    // Density altitude = pressure altitude + temperature correction.
    // The common 120 ft/degree rule is per degree Celsius; these values are Fahrenheit.
    pressure_alt + (120.0 * 5.0 / 9.0) * (temp_f - isa_temp_f)
}

/// Find font file - tries external locations first (for user overrides),
/// then falls back to embedded fonts compiled into the binary.
fn find_font_file(font_name: &str) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
    let ttf = format!("{}.ttf", font_name);

    // Try exe directory
    if let Ok(exe_path) = std::env::current_exe() {
        if let Some(exe_dir) = exe_path.parent() {
            let font_path = exe_dir.join("fonts").join(&ttf);
            if font_path.exists() {
                return Ok(std::fs::read(font_path)?);
            }
        }
    }

    // Try home directory
    if let Some(home) = dirs::home_dir() {
        let font_path = home.join(".ballistics").join("fonts").join(&ttf);
        if font_path.exists() {
            return Ok(std::fs::read(font_path)?);
        }
    }

    // Try working directory
    for prefix in &["./fonts", "../fonts"] {
        let font_path = std::path::Path::new(prefix).join(&ttf);
        if font_path.exists() {
            return Ok(std::fs::read(font_path)?);
        }
    }

    // Try system font directories
    #[cfg(target_os = "linux")]
    {
        for dir in &["/usr/share/fonts", "/usr/local/share/fonts"] {
            if let Some(path) = find_in_directory(dir, &ttf) {
                return Ok(std::fs::read(path)?);
            }
        }
    }

    #[cfg(target_os = "macos")]
    {
        for dir in &["/Library/Fonts", "/System/Library/Fonts"] {
            let font_path = std::path::Path::new(dir).join(&ttf);
            if font_path.exists() {
                return Ok(std::fs::read(font_path)?);
            }
        }
    }

    #[cfg(target_os = "windows")]
    {
        if let Ok(windir) = std::env::var("WINDIR") {
            let font_path = std::path::Path::new(&windir).join("Fonts").join(&ttf);
            if font_path.exists() {
                return Ok(std::fs::read(font_path)?);
            }
        }
    }

    // Fall back to embedded fonts
    match font_name {
        "LiberationSans-Regular" => Ok(FONT_REGULAR.to_vec()),
        "LiberationSans-Bold" => Ok(FONT_BOLD.to_vec()),
        _ => Err(format!("Font {} not found", font_name).into()),
    }
}

/// Recursively search a directory for a font file by name
#[cfg(target_os = "linux")]
fn find_in_directory(dir: &str, filename: &str) -> Option<std::path::PathBuf> {
    let dir_path = std::path::Path::new(dir);
    if !dir_path.is_dir() {
        return None;
    }
    for entry in std::fs::read_dir(dir_path).ok()?.flatten() {
        let path = entry.path();
        if path.is_file() && path.file_name().is_some_and(|n| n == filename) {
            return Some(path);
        }
        if path.is_dir() {
            if let Some(found) = find_in_directory(path.to_str()?, filename) {
                return Some(found);
            }
        }
    }
    None
}

/// Truncate a string for header display, appending "..." if too long
fn truncate_for_header(s: &str, max_chars: usize) -> String {
    // Count/truncate by CHARACTERS, not bytes. The header concatenates user-controlled
    // rifle/location names; byte-slicing a multi-byte UTF-8 string at an offset that isn't a
    // char boundary panics. Identical output for ASCII (byte len == char count).
    if s.chars().count() <= max_chars {
        s.to_string()
    } else if max_chars <= 3 {
        s.chars().take(max_chars).collect()
    } else {
        let head: String = s.chars().take(max_chars - 3).collect();
        format!("{head}...")
    }
}

/// Draw a light gray separator line across the page width
fn draw_separator_line(ops: &mut Vec<Op>, y: f32) {
    ops.push(Op::SetOutlineColor {
        col: Color::Rgb(Rgb::new(0.7, 0.7, 0.7, None)),
    });
    ops.push(Op::SetOutlineThickness { pt: Pt(0.3) });
    ops.push(Op::DrawLine {
        line: Line {
            points: vec![
                LinePoint {
                    p: Point::new(Mm(MARGIN), Mm(y)),
                    bezier: false,
                },
                LinePoint {
                    p: Point::new(Mm(PAGE_WIDTH - MARGIN), Mm(y)),
                    bezier: false,
                },
            ],
            is_closed: false,
        },
    });
}

/// Data rows the two-column table fits on one page at `font_scale`.
///
/// Split out of [`generate_dope_card_pdf`] (which now calls it, so there is exactly one
/// copy of this arithmetic) for callers that must report a page count without holding the
/// document: the bridge's `card.pdf` returns `page_count` in its response, and a second,
/// independent copy of the layout maths there could silently drift from the pagination the
/// generator actually performed.
///
/// `font_scale` is clamped to [`FONT_SCALE_RANGE`] exactly as the generator clamps it.
pub fn dope_card_rows_per_page(font_scale: f32) -> usize {
    let row_height = ROW_HEIGHT * font_scale.clamp(*FONT_SCALE_RANGE.start(), *FONT_SCALE_RANGE.end());
    // Leave space for header/footer + separators
    let usable_height = PAGE_HEIGHT - (2.0 * MARGIN) - 36.0;
    // The clamp above bounds row_height to 2.25..=13.5 mm against a ~223 mm usable
    // height, so this is 16..=52 in practice; `.max(1)` only guarantees the caller's
    // `div_ceil` below can never divide by zero if those page constants are ever edited.
    let visual_rows_per_page = ((usable_height / row_height) as usize).clamp(1, 52);
    // Each visual row shows 2 data points (left + right columns)
    visual_rows_per_page * 2
}

/// Pages [`generate_dope_card_pdf`] emits for `row_count` rows at `font_scale`.
///
/// `0` rows is `0` pages — that call errors rather than producing an empty document, so a
/// zero here is a caller's row set to reject, not a document to describe.
pub fn dope_card_page_count(row_count: usize, font_scale: f32) -> usize {
    row_count.div_ceil(dope_card_rows_per_page(font_scale))
}

/// Generate a dope card PDF matching Glenn's format with row striping.
///
/// `rows` is display-ready per `CardRow`'s convention (already converted to the card's
/// chosen angular unit, and to `range_unit`); `range_unit` only selects the Range
/// column's sub-header text ("Yd" or "M") -- it does not convert `row.range`.
pub fn generate_dope_card_pdf(
    config: &DopeCardConfig,
    rows: &[CardRow],
    range_unit: RangeUnit,
) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
    if rows.is_empty() {
        return Err("generate_dope_card_pdf: rows must not be empty".into());
    }
    let mut doc = PdfDocument::new(&format!("{} Dope Card", config.rifle_name));

    // Load and register fonts
    let font_data = find_font_file("LiberationSans-Regular")?;
    let mut font_warnings = Vec::new();
    let parsed_font = ParsedFont::from_bytes(&font_data, 0, &mut font_warnings)
        .ok_or("Failed to parse LiberationSans-Regular font")?;
    let font = doc.add_font(&parsed_font);

    let font_bold_data = find_font_file("LiberationSans-Bold")?;
    let mut font_bold_warnings = Vec::new();
    let parsed_font_bold = ParsedFont::from_bytes(&font_bold_data, 0, &mut font_bold_warnings)
        .ok_or("Failed to parse LiberationSans-Bold font")?;
    let font_bold = doc.add_font(&parsed_font_bold);

    // Only scale the data table — header/footer stay at base size
    // so they don't overflow or consume disproportionate page space
    let font_scale = config
        .font_scale
        .clamp(*FONT_SCALE_RANGE.start(), *FONT_SCALE_RANGE.end());
    let header_size = HEADER_FONT_SIZE; // UNSCALED
    let table_size = TABLE_FONT_SIZE * font_scale; // SCALED
    let footer_size = FOOTER_FONT_SIZE; // UNSCALED
    let row_height = ROW_HEIGHT * font_scale; // SCALED

    // Pagination lives in `dope_card_rows_per_page`/`dope_card_page_count` so a caller
    // that must report the page count (the bridge's `card.pdf`) reads the same numbers
    // this loop paginates by, instead of reimplementing them.
    let data_rows_per_page = dope_card_rows_per_page(config.font_scale);
    let total_pages = dope_card_page_count(rows.len(), config.font_scale);

    let mut pages = Vec::with_capacity(total_pages);

    for page_num in 0..total_pages {
        let start_idx = page_num * data_rows_per_page;
        let end_idx = std::cmp::min(start_idx + data_rows_per_page, rows.len());
        let page_rows = &rows[start_idx..end_idx];

        let mut ops = Vec::new();

        render_page(
            &mut ops,
            &font,
            &font_bold,
            config,
            page_rows,
            range_unit,
            page_num + 1,
            total_pages,
            header_size,
            table_size,
            footer_size,
            row_height,
            font_scale,
            config.bold_data,
        );

        pages.push(PdfPage::new(Mm(PAGE_WIDTH), Mm(PAGE_HEIGHT), ops));
    }

    let mut save_warnings = Vec::new();
    let bytes = doc
        .with_pages(pages)
        .save(&PdfSaveOptions::default(), &mut save_warnings);
    Ok(bytes)
}

#[allow(clippy::too_many_arguments)] // Fixed-layout renderer keeps its page metrics explicit.
fn render_page(
    ops: &mut Vec<Op>,
    font: &FontId,
    font_bold: &FontId,
    config: &DopeCardConfig,
    rows: &[CardRow],
    range_unit: RangeUnit,
    page: usize,
    _total_pages: usize,
    header_size: f32,
    table_size: f32,
    footer_size: f32,
    row_height: f32,
    font_scale: f32,
    bold_data: bool,
) {
    let mut y = PAGE_HEIGHT - MARGIN;

    // Header line 1 (auto-truncate long text)
    let header1 = truncate_for_header(
        &format!(
            "{} Loc: {} DA:{:.0} ft Pressure:{:.2}/{:.0} Temp:{:.0} Alt:{:.0} Wind:{:.0} Mph",
            config.rifle_name,
            config.location,
            config.density_altitude_ft,
            config.pressure_inhg,
            config.pressure_hpa,
            config.temperature_f,
            config.altitude_ft,
            config.wind_speed_mph
        ),
        77,
    );
    draw_centered_text(ops, font, header_size, y, &header1, COLOR_BLACK);
    y -= 4.0;

    // Header line 2
    let header2 = format!(
        "TargetSpeed:{:.0} Solver: {} - Pg {}",
        config.target_speed_mph, config.solver_mode, page
    );
    draw_centered_text(ops, font, header_size, y, &header2, COLOR_BLACK);
    y -= 1.0;

    // Separator line after header
    draw_separator_line(ops, y);
    y -= 5.0;

    // Table start position
    let table_x = (PAGE_WIDTH - (8.0 * COL_WIDTH)) / 2.0;

    // Draw table header
    draw_table_header(
        ops,
        font_bold,
        table_x,
        y,
        table_size,
        font_scale,
        range_unit,
        &config.elevation_unit_label,
        &config.windage_unit_label,
    );
    y -= row_height;

    // Split rows into left and right columns
    let mid = rows.len().div_ceil(2);
    let (left_rows, right_rows) = rows.split_at(mid);

    // Select font for data rows (bold or regular)
    let data_font = if bold_data { font_bold } else { font };

    // Draw data rows with striping
    for (i, left) in left_rows.iter().enumerate() {
        let right = right_rows.get(i);

        // Draw stripe background for alternating rows
        if i % 2 == 1 {
            draw_row_stripe(ops, table_x, y, 8.0 * COL_WIDTH, row_height);
        }

        // Draw left side data
        draw_data_row(
            ops, data_font, table_x, y, left, true, table_size, font_scale,
            &config.elevation_unit_label, &config.windage_unit_label,
        );

        // Draw right side data
        if let Some(r) = right {
            draw_data_row(
                ops,
                data_font,
                table_x + 4.0 * COL_WIDTH,
                y,
                r,
                false,
                table_size,
                font_scale,
                &config.elevation_unit_label,
                &config.windage_unit_label,
            );
        }

        y -= row_height;
    }

    // Separator line before footer
    draw_separator_line(ops, y - 1.0);
    y -= 5.0;

    // Footer line 1: load data
    let footer1 = format!(
        "Powder:{} Bullet:{} Weight:{:.0}gr BC:{:.3} ({}) Vel:{:.0}fps",
        config.powder,
        config.bullet,
        config.weight_gr,
        config.bc,
        config.drag_model.to_lowercase(),
        config.velocity_fps,
    );
    draw_centered_text(ops, font, footer_size, y, &footer1, COLOR_BLACK);
    y -= 4.0;

    // Footer line 2: timestamp, plus the provenance of the numbers above it. Truncated
    // because both strings are caller-supplied and drawn on every page (the same reason the
    // header truncates), and omitted entirely when empty rather than printing a placeholder.
    let mut footer2 = get_timestamp();
    if !config.engine_version.is_empty() {
        footer2.push_str(&format!(
            " Engine:{}",
            truncate_for_header(&config.engine_version, 24)
        ));
    }
    if !config.table_version.is_empty() {
        footer2.push_str(&format!(
            " Table:{}",
            truncate_for_header(&config.table_version, 24)
        ));
    }
    draw_centered_text(ops, font, footer_size, y, &footer2, COLOR_BLACK);
}

fn draw_row_stripe(ops: &mut Vec<Op>, x: f32, y: f32, width: f32, height: f32) {
    let points = vec![
        LinePoint {
            p: Point::new(Mm(x), Mm(y)),
            bezier: false,
        },
        LinePoint {
            p: Point::new(Mm(x + width), Mm(y)),
            bezier: false,
        },
        LinePoint {
            p: Point::new(Mm(x + width), Mm(y - height)),
            bezier: false,
        },
        LinePoint {
            p: Point::new(Mm(x), Mm(y - height)),
            bezier: false,
        },
    ];

    ops.push(Op::SetFillColor {
        col: Color::Rgb(Rgb::new(
            COLOR_STRIPE.0,
            COLOR_STRIPE.1,
            COLOR_STRIPE.2,
            None,
        )),
    });
    ops.push(Op::DrawPolygon {
        polygon: Polygon {
            rings: vec![PolygonRing { points }],
            mode: PaintMode::Fill,
            winding_order: WindingOrder::NonZero,
        },
    });
}

#[allow(clippy::too_many_arguments)] // Drawing primitive mirrors the PDF text/layout parameters.
fn draw_table_header(
    ops: &mut Vec<Op>,
    font: &FontId,
    x: f32,
    y: f32,
    table_size: f32,
    font_scale: f32,
    range_unit: RangeUnit,
    elevation_unit: &str,
    windage_unit: &str,
) {
    let headers = [
        ("Range", COLOR_BLACK),
        ("Drop", COLOR_RED),
        ("Wind", COLOR_GREEN),
        ("Lead", COLOR_BLUE),
        ("Range", COLOR_BLACK),
        ("Drop", COLOR_RED),
        ("Wind", COLOR_GREEN),
        ("Lead", COLOR_BLUE),
    ];
    // Range columns use the card's range unit (Yd/M); Drop is the elevation unit,
    // Wind/Lead the (possibly different, MBA-1410) windage unit.
    let range_label = range_unit.label();
    let sub_headers = [
        range_label,
        elevation_unit,
        windage_unit,
        windage_unit,
        range_label,
        elevation_unit,
        windage_unit,
        windage_unit,
    ];

    for (i, ((header, color), sub)) in headers.iter().zip(sub_headers.iter()).enumerate() {
        let col_x = x + (i as f32 * COL_WIDTH) + (COL_WIDTH / 2.0);
        draw_text(ops, font, table_size, col_x, y, header, *color, true);
        draw_text(
            ops,
            font,
            table_size - 1.0,
            col_x,
            y - 3.0 * font_scale,
            sub,
            *color,
            true,
        );
    }
}

/// Formats one Drop/Wind/Lead cell value for its column's angular unit (MBA-1410 fold-in
/// of an MBA-1355 backlog minor): whole turret clicks are integers, so `unit_label ==
/// "CLICKS"` prints with no decimal point -- every other unit (MIL/MOA/SMOA/IPHY) keeps
/// the pre-existing one-decimal-place format. Before this fix, a clicks dope card printed
/// e.g. "5.0" instead of "5" for every cell.
///
/// ONE decimal place is the contract for an angular adjustment on EVERY surface that shows
/// these rows, screen included: a turret's resolution is 0.1, and a second decimal is a
/// precision the shooter cannot dial. It is also how a printed card and a screen came to
/// disagree -- 2.4478 MIL read `2.45` on screen against `2.4` on paper, half a click apart
/// on a 0.1-mil turret. Linear columns are unaffected; they keep their own precision.
///
/// The rounding is Rust's: correct rounding of the value's exact binary expansion, with
/// ties-to-even. That differs from rounding the SHORTEST DECIMAL spelling of the same
/// double, which is what most platform number formatters do -- 7.35 is just below the tie
/// and prints `7.3` here, while a decimal half-up formatter prints `7.4`. Any client that
/// must agree with the paper has to match this, which is why
/// `format_adjustment_pins_one_decimal_place_including_the_near_ties` pins the vectors.
fn format_adjustment(value: f64, unit_label: &str) -> String {
    if unit_label.eq_ignore_ascii_case("clicks") {
        format!("{:.0}", value)
    } else {
        format!("{:.1}", value)
    }
}

/// Fix-round I-1: renders an `Option<f64>` dope-card cell honestly. `Some(v)` delegates to
/// `format_adjustment` exactly as before; `None` renders as an em-dash rather than a
/// plausible-looking (but fake) `0.0`. `DopeCardRow` made drop/wind/lead mandatory, so
/// this distinction didn't exist before Task 10 promoted the card onto `CardRow`, whose
/// equivalent fields are `Option`. It matters because Task 11's adaptive card engine --
/// the reason this module accepts `CardRow` at all -- always emits `lead_adj: None`; an
/// `unwrap_or(0.0)` would print a full Lead column of confident-looking zeroes for a card
/// that carries no lead data whatsoever.
fn format_adjustment_cell(value: Option<f64>, unit_label: &str) -> String {
    match value {
        Some(v) => format_adjustment(v, unit_label),
        None => "—".to_string(),
    }
}

/// Formats `CardRow::range` for the dope card's Range column. The pre-Task-10
/// `DopeCardRow` stored range as `u32` yards -- always a bare integer. `CardRow::range`
/// is `f64`: the `trajectory -o pdf` call site now explicitly rounds to whole yards
/// before building each row (`main.rs::dope_card_row_from_sample`, fix-round C-1 --
/// sampled ranges are exact multiples of a metre-denominated sample interval, so their
/// yard conversion is genuinely fractional and would NOT fall into the noise band below
/// on its own), so its rows always land in the integer branch here; a future caller with
/// genuinely fractional ranges (Task 11/12's adaptive cards) renders with one decimal
/// instead, e.g. `417.3` -> "417.3". The `< 0.05` band exists only to absorb ordinary
/// floating-point rounding noise around an already-whole number, not to "round" a
/// meaningfully fractional value down to the nearest integer.
fn format_range(range: f64) -> String {
    let rounded = range.round();
    if (range - rounded).abs() < 0.05 {
        format!("{:.0}", rounded)
    } else {
        format!("{:.1}", range)
    }
}

#[allow(clippy::too_many_arguments)] // Drawing primitive mirrors the PDF text/layout parameters.
fn draw_data_row(
    ops: &mut Vec<Op>,
    font: &FontId,
    x: f32,
    y: f32,
    row: &CardRow,
    _is_left: bool,
    table_size: f32,
    font_scale: f32,
    elevation_unit: &str,
    windage_unit: &str,
) {
    let values = [
        (format_range(row.range), COLOR_BLACK),
        (format_adjustment_cell(row.drop_adj, elevation_unit), COLOR_RED),
        (format_adjustment_cell(row.wind_adj, windage_unit), COLOR_GREEN),
        (format_adjustment_cell(row.lead_adj, windage_unit), COLOR_BLUE),
    ];

    for (i, (value, color)) in values.iter().enumerate() {
        let col_x = x + (i as f32 * COL_WIDTH) + (COL_WIDTH / 2.0);
        draw_text(
            ops,
            font,
            table_size,
            col_x,
            y - 2.5 * font_scale,
            value,
            *color,
            true,
        );
    }
}

/// Push the PDF ops for a run of text, wrapped in its own text section.
///
/// `x`/`y` are the baseline origin in mm (matches printpdf 0.7's `use_text` convention).
fn draw_text_ops(
    ops: &mut Vec<Op>,
    font: &FontId,
    size: f32,
    x: f32,
    y: f32,
    text: &str,
    color: (f32, f32, f32),
) {
    ops.push(Op::StartTextSection);
    ops.push(Op::SetFillColor {
        col: Color::Rgb(Rgb::new(color.0, color.1, color.2, None)),
    });
    ops.push(Op::SetFont {
        font: PdfFontHandle::External(font.clone()),
        size: Pt(size),
    });
    ops.push(Op::SetLineHeight { lh: Pt(size) });
    ops.push(Op::SetTextCursor {
        pos: Point::new(Mm(x), Mm(y)),
    });
    ops.push(Op::ShowText {
        items: vec![TextItem::Text(text.to_string())],
    });
    ops.push(Op::EndTextSection);
}

#[allow(clippy::too_many_arguments)] // Drawing primitive mirrors the PDF text/layout parameters.
fn draw_text(
    ops: &mut Vec<Op>,
    font: &FontId,
    size: f32,
    x: f32,
    y: f32,
    text: &str,
    color: (f32, f32, f32),
    center: bool,
) {
    // Approximate centering by estimating text width (printpdf 0.10 has no simple
    // string-width lookup for an arbitrary font+size, so keep the same heuristic
    // the 0.7-era code used rather than reaching for glyph-level metrics).
    let text_width = if center {
        text.len() as f32 * size * 0.3 // Rough approximation
    } else {
        0.0
    };

    draw_text_ops(ops, font, size, x - text_width / 2.0, y, text, color);
}

fn draw_centered_text(
    ops: &mut Vec<Op>,
    font: &FontId,
    size: f32,
    y: f32,
    text: &str,
    color: (f32, f32, f32),
) {
    // Center on page (same width approximation as draw_text)
    let text_width = text.len() as f32 * size * 0.28;
    let x = (PAGE_WIDTH - text_width) / 2.0;

    draw_text_ops(ops, font, size, x, y, text, color);
}

fn get_timestamp() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};

    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();

    let secs_per_day = 86400u64;
    let secs_per_hour = 3600u64;
    let secs_per_min = 60u64;

    let days_since_epoch = now / secs_per_day;
    let time_of_day = now % secs_per_day;

    let hours = time_of_day / secs_per_hour;
    let minutes = (time_of_day % secs_per_hour) / secs_per_min;
    let seconds = time_of_day % secs_per_min;

    let mut year = 1970;
    let mut remaining_days = days_since_epoch as i64;

    loop {
        let days_in_year = if is_leap_year(year) { 366 } else { 365 };
        if remaining_days < days_in_year {
            break;
        }
        remaining_days -= days_in_year;
        year += 1;
    }

    let days_in_months: [i64; 12] = if is_leap_year(year) {
        [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    } else {
        [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    };

    let month_names = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ];
    let day_names = ["Thu", "Fri", "Sat", "Sun", "Mon", "Tue", "Wed"];

    let mut month = 0;
    for (i, &days) in days_in_months.iter().enumerate() {
        if remaining_days < days {
            month = i;
            break;
        }
        remaining_days -= days;
    }

    let day = remaining_days + 1;
    // The Unix epoch (1970-01-01) was a Thursday and day_names is Thursday-first, so no offset.
    let day_of_week = (days_since_epoch % 7) as usize;

    let (hour_12, am_pm) = if hours == 0 {
        (12, "AM")
    } else if hours < 12 {
        (hours, "AM")
    } else if hours == 12 {
        (12, "PM")
    } else {
        (hours - 12, "PM")
    };

    format!(
        "{} {} {:02} {:02}:{:02}:{:02} {} UTC {}",
        day_names[day_of_week], month_names[month], day, hour_12, minutes, seconds, am_pm, year
    )
}

fn is_leap_year(year: i64) -> bool {
    (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_density_altitude() {
        // MBA-643: Test with Glenn's reference conditions
        // altitude=2500ft, pressure=27.32inHg (station), temp=55°F
        // Glenn's tool: DA ≈ 2835 ft
        let da = calculate_density_altitude(2500.0, 27.32, 55.0);
        assert!(
            da > 2500.0 && da < 3500.0,
            "DA should be ~3000 ft for near-standard conditions, got {}",
            da
        );

        // Higher temp should increase DA
        let da_hot = calculate_density_altitude(2500.0, 27.32, 95.0);
        assert!(da_hot > da, "Higher temp should increase DA");

        // Lower pressure (thinner air) should increase DA
        let da_low_press = calculate_density_altitude(2500.0, 25.0, 55.0);
        assert!(da_low_press > da, "Lower pressure should increase DA");

        // Standard conditions at sea level: DA ≈ 0
        let da_standard = calculate_density_altitude(0.0, 29.92, 59.0);
        assert!(
            da_standard.abs() < 100.0,
            "Standard conditions should give DA near 0, got {}",
            da_standard
        );
    }

    /// MBA-1410 fold-in: a clicks dope-card cell must print as a bare integer ("5"), not
    /// "5.0" -- the trailing-`.0` bug tracked as an MBA-1355 backlog minor. Every other
    /// unit keeps its pre-existing one-decimal-place format.
    #[test]
    fn format_adjustment_drops_the_decimal_for_clicks_only() {
        assert_eq!(format_adjustment(5.0, "CLICKS"), "5");
        assert_eq!(format_adjustment(-3.0, "CLICKS"), "-3");
        assert_eq!(format_adjustment(5.0, "clicks"), "5", "case-insensitive");
        assert_eq!(format_adjustment(2.34, "MIL"), "2.3");
        assert_eq!(format_adjustment(2.34, "MOA"), "2.3");
        assert_eq!(format_adjustment(2.34, "SMOA"), "2.3");
        assert_eq!(format_adjustment(2.34, "IPHY"), "2.3");
    }

    /// The one-decimal contract, pinned as vectors any other surface can be held to.
    ///
    /// Every client that renders these same rows -- the on-screen card in both mobile apps
    /// -- must produce these strings from these doubles, or a shooter's screen and his paper
    /// disagree. The last four are the cases that catch a mismatch: `6.25` is an exact binary
    /// tie (ties-to-even, so `6.2`), `7.35` and `0.15` are just BELOW their ties (`7.3`,
    /// `0.1`) while `2.35` is just above (`2.4`). A formatter that rounds the shortest
    /// decimal spelling half-up prints `6.3`, `7.4` and `0.2` for the first three.
    #[test]
    fn format_adjustment_pins_one_decimal_place_including_the_near_ties() {
        for (value, expected) in [
            (0.0_f64, "0.0"),
            (0.65, "0.7"),
            (1.4551, "1.5"),
            (2.4478, "2.4"),
            (3.575, "3.6"),
            (4.8412, "4.8"),
            (-0.105, "-0.1"),
            (-0.21, "-0.2"),
            (-0.315, "-0.3"),
            (-0.42, "-0.4"),
            (-0.55, "-0.6"),
            (-0.65, "-0.7"),
            (-0.75, "-0.8"),
            (6.25, "6.2"),
            (7.35, "7.3"),
            (0.15, "0.1"),
            (2.35, "2.4"),
        ] {
            assert_eq!(format_adjustment(value, "MIL"), expected, "{value:?}");
            assert_eq!(format_adjustment(value, "MOA"), expected, "{value:?}");
        }
    }

    /// Fix-round I-1: a missing column must render as an honest em-dash, never a
    /// plausible-looking fake `0.0` -- pinned for both a decimal unit and a clicks unit,
    /// since a naive fix might special-case `None` only inside one branch of
    /// `format_adjustment`'s clicks/non-clicks split.
    #[test]
    fn format_adjustment_cell_renders_none_as_an_em_dash_not_a_fake_zero() {
        assert_eq!(format_adjustment_cell(Some(2.34), "MIL"), "2.3");
        assert_eq!(format_adjustment_cell(Some(5.0), "CLICKS"), "5");
        assert_eq!(format_adjustment_cell(None, "MIL"), "—");
        assert_eq!(format_adjustment_cell(None, "CLICKS"), "—");
    }

    #[test]
    fn density_altitude_uses_published_nws_pressure_altitude_equation() {
        // 20.670988150011322 inHg is exactly 700 hPa under the standard conversion.
        // Keep the public-unit fixture literal independent of the production conversion constant.
        let pressure_inhg = 20.670_988_150_011_322;
        let density_altitude = calculate_density_altitude(0.0, pressure_inhg, 55.0);

        assert!((density_altitude - 11_962.774_017_764_264).abs() < 1e-6);

        let standard_pressure_inhg = 29.921_255_347_141_39;
        let standard_density_altitude =
            calculate_density_altitude(0.0, standard_pressure_inhg, 59.0);
        assert!(standard_density_altitude.abs() < 1e-9);
    }

    // -----------------------------------------------------------------------------------
    // Task 10 (0.33.0 decision-support Plan B): CardRow-based range formatting + PDF
    // generation smoke tests.
    // -----------------------------------------------------------------------------------

    /// Pinned brief examples: a legacy-integer range renders with no decimal point, a
    /// genuinely fractional one (as Task 11's adaptive engine can now produce) renders
    /// with exactly one.
    #[test]
    fn format_range_pinned_examples() {
        assert_eq!(format_range(400.0), "400");
        assert_eq!(format_range(417.3), "417.3");
    }

    /// The pre-Task-10 call site rounded to `u32` before storing, so its yard conversion
    /// was never exactly integral (floating-point noise). The new `f64`-range path must
    /// still collapse that noise to the same bare integer -- this is the `< 0.05`
    /// tolerance's whole reason to exist.
    #[test]
    fn format_range_absorbs_floating_point_noise_around_a_whole_number() {
        assert_eq!(format_range(400.02), "400");
        assert_eq!(format_range(399.98), "400");
        assert_eq!(format_range(100.0 + 1e-9), "100");
        // The band rule now exists in two copies (this function and
        // `format_adaptive_card_range` in `main.rs`), so these fixtures are the sync
        // mechanism between them -- more coverage right at and near the +-0.05 edge.
        assert_eq!(format_range(400.04), "400");
        assert_eq!(format_range(399.96), "400");
        assert_eq!(format_range(0.0), "0");
    }

    #[test]
    fn format_range_renders_one_decimal_outside_the_noise_band() {
        assert_eq!(format_range(417.3), "417.3");
        assert_eq!(format_range(100.5), "100.5");
        assert_eq!(format_range(400.08), "400.1");
        // 400.05 is exactly the band's edge in decimal, but not in binary: the nearest f64 to
        // 400.05 is a hair above it, so `(range - rounded).abs() < 0.05` is false and this
        // renders with one decimal rather than collapsing to "400".
        assert_eq!(format_range(400.05), "400.1");
    }

    fn test_config() -> DopeCardConfig {
        DopeCardConfig {
            rifle_name: "Test Rifle".to_string(),
            location: "Test Range".to_string(),
            density_altitude_ft: 1500.0,
            pressure_inhg: 29.92,
            pressure_hpa: 1013.25,
            temperature_f: 59.0,
            altitude_ft: 1000.0,
            wind_speed_mph: 5.0,
            target_speed_mph: 0.0,
            solver_mode: "offline".to_string(),
            powder: "Test Powder".to_string(),
            bullet: "Test Bullet".to_string(),
            weight_gr: 175.0,
            bc: 0.5,
            drag_model: "g7".to_string(),
            velocity_fps: 2700.0,
            font_scale: 1.0,
            bold_data: false,
            elevation_unit_label: "MIL".to_string(),
            windage_unit_label: "MIL".to_string(),
            engine_version: "0.0.0-test".to_string(),
            table_version: String::new(),
        }
    }

    fn test_row(range: f64) -> CardRow {
        CardRow {
            range,
            drop_linear: None,
            drop_adj: Some(1.0),
            come_up: None,
            wind_linear: None,
            wind_adj: Some(0.5),
            velocity: None,
            energy: None,
            time: None,
            lead_adj: Some(0.2),
            wind_columns: Vec::new(),
        }
    }

    fn assert_valid_pdf_bytes(bytes: &[u8]) {
        assert!(!bytes.is_empty(), "PDF output must not be empty");
        assert_eq!(&bytes[..5], b"%PDF-", "output must start with a PDF header");
    }

    #[test]
    fn generate_dope_card_pdf_is_non_empty_for_a_single_row() {
        let config = test_config();
        let rows = vec![test_row(100.0)];
        let bytes = generate_dope_card_pdf(&config, &rows, RangeUnit::Yards)
            .expect("single-row dope card should generate");
        assert_valid_pdf_bytes(&bytes);
    }

    /// 120 rows forces the two-column, ~98-rows-per-page layout past a single page
    /// (`data_rows_per_page` is well under 120 at the default font scale), exercising the
    /// pagination path `generate_dope_card_pdf` walks via `total_pages`/`page_rows`.
    #[test]
    fn generate_dope_card_pdf_is_non_empty_for_120_rows_pagination() {
        let config = test_config();
        let rows: Vec<CardRow> = (0..120)
            .map(|i| test_row(100.0 + i as f64 * 25.0))
            .collect();
        let bytes = generate_dope_card_pdf(&config, &rows, RangeUnit::Meters)
            .expect("120-row dope card should generate and paginate");
        assert_valid_pdf_bytes(&bytes);
    }

    /// Fix-round I-1: Task 11's adaptive card engine -- the reason this module accepts
    /// `CardRow` at all -- always emits `lead_adj: None`. The full renderer must still
    /// produce a valid PDF for that row shape (via `format_adjustment_cell`'s em-dash,
    /// not `unwrap_or(0.0)`'s fake zero); this exercises `draw_data_row` end to end,
    /// while `format_adjustment_cell_renders_none_as_an_em_dash_not_a_fake_zero` above
    /// pins the exact string.
    #[test]
    fn generate_dope_card_pdf_succeeds_when_lead_adj_is_none() {
        let config = test_config();
        let mut row = test_row(100.0);
        row.lead_adj = None;
        let bytes = generate_dope_card_pdf(&config, &[row], RangeUnit::Yards)
            .expect("a row missing lead_adj must still render, not error");
        assert_valid_pdf_bytes(&bytes);
    }

    /// A zero-page PDF is not a useful answer from a public library API: both in-tree callers
    /// (`trajectory -o pdf`, `adaptive-card -o pdf`) already guard against empty rows
    /// themselves, but a binding calling this function directly should get a named error, not
    /// a silently-empty document.
    #[test]
    fn generate_dope_card_pdf_rejects_empty_rows() {
        let config = test_config();
        let err = generate_dope_card_pdf(&config, &[], RangeUnit::Yards)
            .expect_err("an empty row set must not silently produce a PDF");
        assert!(err.to_string().contains("rows"), "{err}");
    }

    /// The pagination `generate_dope_card_pdf` uses is now a public function (the bridge's
    /// `card.pdf` reports a page count from it), so pin its numbers: they are part of that
    /// response contract, and a silent change here would silently mislabel every PDF the
    /// bridge returns.
    #[test]
    fn dope_card_rows_per_page_pins_the_preset_scales() {
        // 279.4 - 20 - 36 = 223.4 mm usable / (4.5 * scale) mm per visual row, floored,
        // capped at 52 visual rows, doubled for the two-column layout.
        // The 52-visual-row cap bites below ~0.955 scale (223.4 / 4.5 / 52), so Small and
        // every scale under it are cap-limited rather than height-limited — which is why
        // Small and the 0.5 floor agree.
        assert_eq!(dope_card_rows_per_page(FontSizePreset::Small.scale()), 104);
        assert_eq!(dope_card_rows_per_page(FontSizePreset::Medium.scale()), 98);
        assert_eq!(dope_card_rows_per_page(FontSizePreset::Large.scale()), 70);
        assert_eq!(dope_card_rows_per_page(0.5), 104);
        // Out-of-band scales are clamped, exactly as the generator clamps them.
        assert_eq!(dope_card_rows_per_page(0.01), dope_card_rows_per_page(0.5));
        assert_eq!(dope_card_rows_per_page(99.0), dope_card_rows_per_page(3.0));
    }

    #[test]
    fn dope_card_page_count_rounds_up_and_reports_zero_for_no_rows() {
        let per_page = dope_card_rows_per_page(1.0);
        assert_eq!(dope_card_page_count(0, 1.0), 0);
        assert_eq!(dope_card_page_count(1, 1.0), 1);
        assert_eq!(dope_card_page_count(per_page, 1.0), 1);
        assert_eq!(dope_card_page_count(per_page + 1, 1.0), 2);
        assert_eq!(dope_card_page_count(per_page * 3, 1.0), 3);
    }

    /// The generator must paginate by the same function it reports: a PDF built from
    /// `per_page + 1` rows has to carry a second page's worth of ops, which shows up as a
    /// materially larger document than the single-page case.
    #[test]
    fn generated_pdf_grows_when_page_count_grows() {
        let config = test_config();
        let per_page = dope_card_rows_per_page(config.font_scale);
        let rows: Vec<CardRow> = (0..per_page).map(|i| test_row(100.0 + i as f64)).collect();
        let one_page = generate_dope_card_pdf(&config, &rows, RangeUnit::Yards).unwrap();
        assert_eq!(dope_card_page_count(rows.len(), config.font_scale), 1);

        let mut rows_plus_one = rows.clone();
        rows_plus_one.push(test_row(100.0 + per_page as f64));
        let two_pages = generate_dope_card_pdf(&config, &rows_plus_one, RangeUnit::Yards).unwrap();
        assert_eq!(
            dope_card_page_count(rows_plus_one.len(), config.font_scale),
            2
        );
        assert!(
            two_pages.len() > one_page.len(),
            "a second page must add bytes: {} vs {}",
            one_page.len(),
            two_pages.len()
        );
    }

    /// `RangeUnit::label()` untested was Task 10 review Minor #3: a swapped `Yd`/`M` would
    /// ship silently on `adaptive-card -o pdf`'s Range column sub-header.
    #[test]
    fn range_unit_label_matches_its_variant() {
        assert_eq!(RangeUnit::Yards.label(), "Yd");
        assert_eq!(RangeUnit::Meters.label(), "M");
    }
}
