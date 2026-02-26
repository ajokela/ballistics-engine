//! PDF Dope Card Generation Module
//!
//! Generates printable dope cards in Glenn's proven field-ready format.
//! Format: Two-column layout with Range, Drop MIL, Wind MIL, Lead MIL
//! Color coding: Black=Range, Red=Drop, Green=Wind, Blue=Lead
//! Row striping for improved readability.

use printpdf::*;
use std::io::BufWriter;

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

    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "small" | "s" => Some(Self::Small),
            "medium" | "m" => Some(Self::Medium),
            "large" | "l" => Some(Self::Large),
            _ => None,
        }
    }
}

/// A single row in the dope card table
#[derive(Debug, Clone)]
pub struct DopeCardRow {
    pub range_yd: u32,
    pub drop_mil: f64,
    pub wind_mil: f64,
    pub lead_mil: f64,
}

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

/// Convert drop in yards to MILs
pub fn yards_to_mil(drop_yd: f64, range_yd: f64) -> f64 {
    if range_yd < 1.0 {
        return 0.0;
    }
    (drop_yd / range_yd) * 1000.0
}

/// Calculate lead MIL for moving target
pub fn calculate_lead_mil(target_speed_mph: f64, time_of_flight_s: f64, range_yd: f64) -> f64 {
    if range_yd < 1.0 || target_speed_mph < 0.001 {
        return 0.0;
    }
    let target_speed_yps = target_speed_mph * 1760.0 / 3600.0;
    let target_movement_yd = target_speed_yps * time_of_flight_s;
    yards_to_mil(target_movement_yd, range_yd)
}

/// Calculate density altitude from environmental conditions
///
/// MBA-643: Fixed to interpret pressure as STATION PRESSURE (actual local pressure),
/// not altimeter setting (sea-level corrected). This matches how weather stations
/// and most ballistic tools report pressure.
///
/// Formula: PA = 145442 * (1 - (P/29.92)^0.190284)
///          DA = PA + 120 * (OAT - ISA_temp)
pub fn calculate_density_altitude(_altitude_ft: f64, pressure_inhg: f64, temp_f: f64) -> f64 {
    // Calculate pressure altitude from station pressure using barometric formula
    // This is the altitude in standard atmosphere that has the given pressure
    let pressure_alt = 145442.0 * (1.0 - (pressure_inhg / 29.92_f64).powf(0.190284));

    // ISA temperature at pressure altitude (lapse rate: 3.57°F per 1000 ft)
    let isa_temp_f = 59.0 - (pressure_alt / 1000.0) * 3.57;

    // Density altitude = pressure altitude + temperature correction
    pressure_alt + 120.0 * (temp_f - isa_temp_f)
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
        if path.is_file() && path.file_name().map_or(false, |n| n == filename) {
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
    if s.len() <= max_chars {
        s.to_string()
    } else if max_chars <= 3 {
        s[..max_chars].to_string()
    } else {
        format!("{}...", &s[..max_chars - 3])
    }
}

/// Draw a light gray separator line across the page width
fn draw_separator_line(layer: &PdfLayerReference, y: f32) {
    let line = Line {
        points: vec![
            (Point::new(Mm(MARGIN), Mm(y)), false),
            (Point::new(Mm(PAGE_WIDTH - MARGIN), Mm(y)), false),
        ],
        is_closed: false,
    };
    layer.set_outline_color(Color::Rgb(Rgb::new(0.7, 0.7, 0.7, None)));
    layer.set_outline_thickness(0.3);
    layer.add_line(line);
}

/// Generate a dope card PDF matching Glenn's format with row striping
pub fn generate_dope_card_pdf(
    config: &DopeCardConfig,
    rows: &[DopeCardRow],
) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
    let (doc, page1, layer1) =
        PdfDocument::new(&format!("{} Dope Card", config.rifle_name), Mm(PAGE_WIDTH), Mm(PAGE_HEIGHT), "Layer 1");

    // Load font
    let font_data = find_font_file("LiberationSans-Regular")?;
    let font = doc.add_external_font(&*font_data)?;

    let font_bold_data = find_font_file("LiberationSans-Bold")?;
    let font_bold = doc.add_external_font(&*font_bold_data)?;

    // Only scale the data table — header/footer stay at base size
    // so they don't overflow or consume disproportionate page space
    let font_scale = config.font_scale.clamp(0.5, 3.0);
    let header_size = HEADER_FONT_SIZE;           // UNSCALED
    let table_size = TABLE_FONT_SIZE * font_scale; // SCALED
    let footer_size = FOOTER_FONT_SIZE;           // UNSCALED
    let row_height = ROW_HEIGHT * font_scale;     // SCALED

    // Calculate visual rows per page (accounting for header and footer)
    // Each visual row shows 2 data points (left + right columns)
    let usable_height = PAGE_HEIGHT - (2.0 * MARGIN) - 36.0; // Leave space for header/footer + separators
    let visual_rows_per_page = ((usable_height / row_height) as usize).min(52);
    let data_rows_per_page = visual_rows_per_page * 2; // Two-column layout
    let total_pages = (rows.len() + data_rows_per_page - 1) / data_rows_per_page;

    for page_num in 0..total_pages {
        let start_idx = page_num * data_rows_per_page;
        let end_idx = std::cmp::min(start_idx + data_rows_per_page, rows.len());
        let page_rows = &rows[start_idx..end_idx];

        let (current_page, current_layer) = if page_num == 0 {
            (page1, layer1)
        } else {
            doc.add_page(Mm(PAGE_WIDTH), Mm(PAGE_HEIGHT), &format!("Page {}", page_num + 1))
        };

        let layer = doc.get_page(current_page).get_layer(current_layer);

        render_page(&layer, &font, &font_bold, config, page_rows, page_num + 1, total_pages,
                     header_size, table_size, footer_size, row_height, font_scale, config.bold_data)?;
    }

    let mut buffer = Vec::new();
    doc.save(&mut BufWriter::new(&mut buffer))?;
    Ok(buffer)
}

fn render_page(
    layer: &PdfLayerReference,
    font: &IndirectFontRef,
    font_bold: &IndirectFontRef,
    config: &DopeCardConfig,
    rows: &[DopeCardRow],
    page: usize,
    _total_pages: usize,
    header_size: f32,
    table_size: f32,
    footer_size: f32,
    row_height: f32,
    font_scale: f32,
    bold_data: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut y = PAGE_HEIGHT - MARGIN;

    // Header line 1 (auto-truncate long text)
    let header1 = truncate_for_header(&format!(
        "{} Loc: {} DA:{:.0} ft Pressure:{:.2}/{:.0} Temp:{:.0} Alt:{:.0} Wind:{:.0} Mph",
        config.rifle_name,
        config.location,
        config.density_altitude_ft,
        config.pressure_inhg,
        config.pressure_hpa,
        config.temperature_f,
        config.altitude_ft,
        config.wind_speed_mph
    ), 77);
    draw_centered_text(layer, font, header_size, y, &header1, COLOR_BLACK);
    y -= 4.0;

    // Header line 2
    let header2 = format!(
        "TargetSpeed:{:.0} Solver: {} - Pg {}",
        config.target_speed_mph, config.solver_mode, page
    );
    draw_centered_text(layer, font, header_size, y, &header2, COLOR_BLACK);
    y -= 1.0;

    // Separator line after header
    draw_separator_line(layer, y);
    y -= 5.0;

    // Table start position
    let table_x = (PAGE_WIDTH - (8.0 * COL_WIDTH)) / 2.0;

    // Draw table header
    draw_table_header(layer, font_bold, table_x, y, table_size, font_scale);
    y -= row_height;

    // Split rows into left and right columns
    let mid = (rows.len() + 1) / 2;
    let (left_rows, right_rows) = rows.split_at(mid);

    // Select font for data rows (bold or regular)
    let data_font = if bold_data { font_bold } else { font };

    // Draw data rows with striping
    for i in 0..left_rows.len() {
        let left = &left_rows[i];
        let right = right_rows.get(i);

        // Draw stripe background for alternating rows
        if i % 2 == 1 {
            draw_row_stripe(layer, table_x, y, 8.0 * COL_WIDTH, row_height);
        }

        // Draw left side data
        draw_data_row(layer, data_font, table_x, y, left, true, table_size, font_scale);

        // Draw right side data
        if let Some(r) = right {
            draw_data_row(layer, data_font, table_x + 4.0 * COL_WIDTH, y, r, false, table_size, font_scale);
        }

        y -= row_height;
    }

    // Separator line before footer
    draw_separator_line(layer, y - 1.0);
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
    draw_centered_text(layer, font, footer_size, y, &footer1, COLOR_BLACK);
    y -= 4.0;

    // Footer line 2: timestamp
    let timestamp = get_timestamp();
    draw_centered_text(layer, font, footer_size, y, &timestamp, COLOR_BLACK);

    Ok(())
}

fn draw_row_stripe(layer: &PdfLayerReference, x: f32, y: f32, width: f32, height: f32) {
    use printpdf::path::PaintMode;

    let points = vec![
        (Point::new(Mm(x), Mm(y)), false),
        (Point::new(Mm(x + width), Mm(y)), false),
        (Point::new(Mm(x + width), Mm(y - height)), false),
        (Point::new(Mm(x), Mm(y - height)), false),
    ];

    let rect = Polygon {
        rings: vec![points],
        mode: PaintMode::Fill,
        winding_order: printpdf::path::WindingOrder::NonZero,
    };

    layer.set_fill_color(Color::Rgb(Rgb::new(
        COLOR_STRIPE.0,
        COLOR_STRIPE.1,
        COLOR_STRIPE.2,
        None,
    )));
    layer.add_polygon(rect);
}

fn draw_table_header(layer: &PdfLayerReference, font: &IndirectFontRef, x: f32, y: f32, table_size: f32, font_scale: f32) {
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
    let sub_headers = ["Yd", "MIL", "MIL", "MIL", "Yd", "MIL", "MIL", "MIL"];

    for (i, ((header, color), sub)) in headers.iter().zip(sub_headers.iter()).enumerate() {
        let col_x = x + (i as f32 * COL_WIDTH) + (COL_WIDTH / 2.0);
        draw_text(layer, font, table_size, col_x, y, header, *color, true);
        draw_text(layer, font, table_size - 1.0, col_x, y - 3.0 * font_scale, sub, *color, true);
    }
}

fn draw_data_row(layer: &PdfLayerReference, font: &IndirectFontRef, x: f32, y: f32, row: &DopeCardRow, _is_left: bool, table_size: f32, font_scale: f32) {
    let values = [
        (row.range_yd.to_string(), COLOR_BLACK),
        (format!("{:.1}", row.drop_mil), COLOR_RED),
        (format!("{:.1}", row.wind_mil), COLOR_GREEN),
        (format!("{:.1}", row.lead_mil), COLOR_BLUE),
    ];

    for (i, (value, color)) in values.iter().enumerate() {
        let col_x = x + (i as f32 * COL_WIDTH) + (COL_WIDTH / 2.0);
        draw_text(layer, font, table_size, col_x, y - 2.5 * font_scale, value, *color, true);
    }
}

fn draw_text(
    layer: &PdfLayerReference,
    font: &IndirectFontRef,
    size: f32,
    x: f32,
    y: f32,
    text: &str,
    color: (f32, f32, f32),
    center: bool,
) {
    layer.set_fill_color(Color::Rgb(Rgb::new(color.0, color.1, color.2, None)));

    // Approximate centering by estimating text width
    let text_width = if center {
        text.len() as f32 * size * 0.3  // Rough approximation
    } else {
        0.0
    };

    layer.use_text(text, size, Mm(x - text_width / 2.0), Mm(y), font);
}

fn draw_centered_text(
    layer: &PdfLayerReference,
    font: &IndirectFontRef,
    size: f32,
    y: f32,
    text: &str,
    color: (f32, f32, f32),
) {
    layer.set_fill_color(Color::Rgb(Rgb::new(color.0, color.1, color.2, None)));

    // Center on page
    let text_width = text.len() as f32 * size * 0.28;
    let x = (PAGE_WIDTH - text_width) / 2.0;

    layer.use_text(text, size, Mm(x), Mm(y), font);
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
    let day_of_week = ((days_since_epoch + 4) % 7) as usize;

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
        day_names[day_of_week],
        month_names[month],
        day,
        hour_12,
        minutes,
        seconds,
        am_pm,
        year
    )
}

fn is_leap_year(year: i64) -> bool {
    (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_yards_to_mil() {
        let mil = yards_to_mil(0.1, 100.0);
        assert!((mil - 1.0).abs() < 0.01);

        let mil = yards_to_mil(1.78, 500.0);
        assert!((mil - 3.56).abs() < 0.1);
    }

    #[test]
    fn test_lead_mil() {
        let lead = calculate_lead_mil(4.0, 0.73, 500.0);
        assert!((lead - 2.86).abs() < 0.2);
    }

    #[test]
    fn test_density_altitude() {
        // MBA-643: Test with Glenn's reference conditions
        // altitude=2500ft, pressure=27.32inHg (station), temp=55°F
        // Glenn's tool: DA ≈ 2835 ft
        let da = calculate_density_altitude(2500.0, 27.32, 55.0);
        assert!(da > 2500.0 && da < 3500.0,
            "DA should be ~3000 ft for near-standard conditions, got {}", da);

        // Higher temp should increase DA
        let da_hot = calculate_density_altitude(2500.0, 27.32, 95.0);
        assert!(da_hot > da, "Higher temp should increase DA");

        // Lower pressure (thinner air) should increase DA
        let da_low_press = calculate_density_altitude(2500.0, 25.0, 55.0);
        assert!(da_low_press > da, "Lower pressure should increase DA");

        // Standard conditions at sea level: DA ≈ 0
        let da_standard = calculate_density_altitude(0.0, 29.92, 59.0);
        assert!(da_standard.abs() < 100.0,
            "Standard conditions should give DA near 0, got {}", da_standard);
    }
}
