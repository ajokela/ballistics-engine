//! Zero-dependency terminal chart rendering for `trajectory --plot` (MBA-1320).
//!
//! CLI-only presentation code: this module is a `main.rs`-local `mod` (like
//! [`crate::pdf_dope_card`] and [`crate::solve_json_command`]), not part of the published
//! `ballistics_engine` library API — bindings/consumers embedding the crate have no use for
//! a terminal chart renderer.
//!
//! Two canvas backends share one dot-addressed API (`new`, `width_dots`/`height_dots`,
//! `set`, `render`):
//!
//! - [`BrailleCanvas`] packs a 2 (wide) x 4 (tall) grid of dots into a single Unicode
//!   Braille Patterns character (`U+2800`..=`U+28FF`), giving roughly 2x the horizontal and
//!   4x the vertical resolution of one terminal character cell.
//! - [`AsciiCanvas`] uses the *same* 2x4 dot-per-cell addressing so callers (in particular
//!   [`render_chart`]) don't need to special-case it, but collapses each cell to a single
//!   `'*'` (any dot set) or `' '` (no dots set) — for terminals/fonts without braille glyph
//!   coverage.
//!
//! Deliberately monochrome: no ANSI color/SGR codes appear anywhere in this module. That
//! sidesteps `NO_COLOR` (<https://no-color.org/>) entirely — there is nothing to suppress —
//! and keeps output byte-identical whether the terminal honors color, redirects to a file,
//! or is a dumb pipe. Multiple series plotted on one canvas are therefore visually
//! indistinguishable where their dots overlap; [`render_chart`] lists every series label in
//! the frame instead of color-coding them. Callers that need series to stay visually
//! separable should render one series per chart (which is what `trajectory --plot` does:
//! one chart for drop, one for lateral drift).

/// Bit weight of each of the 8 sub-dot positions within one braille cell, indexed
/// `[row][col]` with row `0..4` (top-to-bottom) and col `0..2` (left-to-right). This is the
/// standard Unicode braille dot numbering:
///
/// ```text
/// 1 4      0x01 0x08
/// 2 5  ->  0x02 0x10
/// 3 6      0x04 0x20
/// 7 8      0x40 0x80
/// ```
const BRAILLE_BITS: [[u8; 2]; 4] = [[0x01, 0x08], [0x02, 0x10], [0x04, 0x20], [0x40, 0x80]];

/// First codepoint of the Unicode Braille Patterns block; a cell's mask (0..=255) added to
/// this base gives the glyph for that cell.
const BRAILLE_BASE: u32 = 0x2800;

/// A monochrome canvas addressed in DOT coordinates, shared by [`BrailleCanvas`] and
/// [`AsciiCanvas`]: `(0, 0)` is the top-left dot, `x` grows right, `y` grows down. Each
/// implementor packs a 2 (wide) x 4 (tall) block of dots into one printable character
/// cell — `width_dots() == width_cells * 2`, `height_dots() == height_cells * 4`.
trait DotCanvas {
    fn width_dots(&self) -> usize;
    fn height_dots(&self) -> usize;
    /// Turn on the dot at `(x, y)`. Coordinates outside the canvas are silently dropped
    /// (callers are expected to clamp before calling; this just keeps the API panic-free).
    fn set(&mut self, x: usize, y: usize);
    /// Render to a string of `height_cells` lines, each `width_cells` characters wide, no
    /// trailing newline.
    fn render(&self) -> String;
}

/// Braille-dot terminal canvas: each character cell packs a 2x4 grid of dots into one
/// Unicode Braille Patterns codepoint.
#[derive(Debug, Clone)]
pub struct BrailleCanvas {
    width_cells: usize,
    height_cells: usize,
    cells: Vec<u8>,
}

impl BrailleCanvas {
    pub fn new(width_cells: usize, height_cells: usize) -> Self {
        Self {
            width_cells,
            height_cells,
            cells: vec![0u8; width_cells * height_cells],
        }
    }

    pub fn width_dots(&self) -> usize {
        self.width_cells * 2
    }

    pub fn height_dots(&self) -> usize {
        self.height_cells * 4
    }

    /// Turn on the dot at dot-coordinate `(x, y)`; `(0, 0)` is top-left. Coordinates
    /// outside the canvas are silently dropped.
    pub fn set(&mut self, x: usize, y: usize) {
        if x >= self.width_dots() || y >= self.height_dots() {
            return;
        }
        let cell_x = x / 2;
        let cell_y = y / 4;
        let bit = BRAILLE_BITS[y % 4][x % 2];
        self.cells[cell_y * self.width_cells + cell_x] |= bit;
    }

    /// Render to a string of `height_cells` lines, each `width_cells` characters wide (no
    /// trailing newline). A cell with no dots set renders as `U+2800` (BRAILLE PATTERN
    /// BLANK), not an ASCII space — it is a real glyph, not empty output.
    pub fn render(&self) -> String {
        let mut out = String::with_capacity((self.width_cells + 1) * self.height_cells);
        for row in 0..self.height_cells {
            if row > 0 {
                out.push('\n');
            }
            for col in 0..self.width_cells {
                let mask = self.cells[row * self.width_cells + col];
                // Every value 0..=255 maps to a valid codepoint in the braille block, so
                // this can never actually fall back to '?'.
                let ch = char::from_u32(BRAILLE_BASE + mask as u32).unwrap_or('?');
                out.push(ch);
            }
        }
        out
    }
}

impl DotCanvas for BrailleCanvas {
    fn width_dots(&self) -> usize {
        BrailleCanvas::width_dots(self)
    }
    fn height_dots(&self) -> usize {
        BrailleCanvas::height_dots(self)
    }
    fn set(&mut self, x: usize, y: usize) {
        BrailleCanvas::set(self, x, y)
    }
    fn render(&self) -> String {
        BrailleCanvas::render(self)
    }
}

/// ASCII fallback canvas: the same dot-addressed API as [`BrailleCanvas`] (a 2x4 dot grid
/// per character cell), but rendered with a single `'*'` whenever ANY dot in the cell is
/// set, `' '` otherwise — for terminals/fonts without braille glyph coverage.
#[derive(Debug, Clone)]
pub struct AsciiCanvas {
    width_cells: usize,
    height_cells: usize,
    cells: Vec<bool>,
}

impl AsciiCanvas {
    pub fn new(width_cells: usize, height_cells: usize) -> Self {
        Self {
            width_cells,
            height_cells,
            cells: vec![false; width_cells * height_cells],
        }
    }

    pub fn width_dots(&self) -> usize {
        self.width_cells * 2
    }

    pub fn height_dots(&self) -> usize {
        self.height_cells * 4
    }

    pub fn set(&mut self, x: usize, y: usize) {
        if x >= self.width_dots() || y >= self.height_dots() {
            return;
        }
        let cell_x = x / 2;
        let cell_y = y / 4;
        self.cells[cell_y * self.width_cells + cell_x] = true;
    }

    pub fn render(&self) -> String {
        let mut out = String::with_capacity((self.width_cells + 1) * self.height_cells);
        for row in 0..self.height_cells {
            if row > 0 {
                out.push('\n');
            }
            for col in 0..self.width_cells {
                out.push(if self.cells[row * self.width_cells + col] {
                    '*'
                } else {
                    ' '
                });
            }
        }
        out
    }
}

impl DotCanvas for AsciiCanvas {
    fn width_dots(&self) -> usize {
        AsciiCanvas::width_dots(self)
    }
    fn height_dots(&self) -> usize {
        AsciiCanvas::height_dots(self)
    }
    fn set(&mut self, x: usize, y: usize) {
        AsciiCanvas::set(self, x, y)
    }
    fn render(&self) -> String {
        AsciiCanvas::render(self)
    }
}

/// Which canvas backend [`render_chart`] should use.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CanvasStyle {
    Braille,
    Ascii,
}

fn make_canvas(style: CanvasStyle, width_cells: usize, height_cells: usize) -> Box<dyn DotCanvas> {
    match style {
        CanvasStyle::Braille => Box::new(BrailleCanvas::new(width_cells, height_cells)),
        CanvasStyle::Ascii => Box::new(AsciiCanvas::new(width_cells, height_cells)),
    }
}

/// Build one border line exactly `width` characters between `left` and `right`: as much of
/// `text` as fits, then `fill` out to `width`. `text` longer than `width` is truncated (no
/// ellipsis — this is a fixed-width terminal frame, not prose).
fn frame_line(left: char, fill: char, right: char, text: &str, width: usize) -> String {
    let mut s = String::with_capacity(width + 2);
    s.push(left);
    let mut used = 0usize;
    for c in text.chars() {
        if used >= width {
            break;
        }
        s.push(c);
        used += 1;
    }
    for _ in used..width {
        s.push(fill);
    }
    s.push(right);
    s
}

/// Render one or more named `(label, points)` series into a single framed chart.
///
/// - Axis ranges are auto-scaled from the data: min/max `x` and min/max `y` across every
///   series. The `y` range and every series label are printed on the top border; the `x`
///   range is printed on the bottom border — "framed chart string with x/y axis ranges
///   printed on the border", per the ticket.
/// - Points are plotted as dots only; there is no line interpolation between consecutive
///   samples. Dense series (trajectory integration output, which this is built for, easily
///   has hundreds to thousands of points over a 72-cell-wide canvas) render as a
///   continuous-looking curve purely from sample density; sparse series will look like
///   literally scattered dots. Adding Bresenham-style line segments was considered and
///   dropped: it's not needed for the trajectory use case and would be extra surface to get
///   right for no visible benefit there.
/// - `width_cells`/`height_cells` size only the PLOT AREA. The returned string is
///   `width_cells + 2` characters wide (a left/right frame column) and `height_cells + 2`
///   lines tall (top/bottom frame rows), with no trailing newline.
/// - A degenerate axis (every point shares the same `x`, or the same `y`) is expanded by
///   ±0.5 around that single value so the scale never divides by zero.
/// - Non-finite (`NaN`/`inf`) coordinates are skipped entirely — they neither affect the
///   autoscale range nor get plotted — rather than corrupting the whole chart's scale.
pub fn render_chart(
    series: &[(&str, &[(f64, f64)])],
    width_cells: usize,
    height_cells: usize,
    style: CanvasStyle,
) -> String {
    let finite_points = || {
        series
            .iter()
            .flat_map(|(_, pts)| pts.iter())
            .filter(|(x, y)| x.is_finite() && y.is_finite())
    };

    let (mut xmin, mut xmax) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut ymin, mut ymax) = (f64::INFINITY, f64::NEG_INFINITY);
    for &(x, y) in finite_points() {
        xmin = xmin.min(x);
        xmax = xmax.max(x);
        ymin = ymin.min(y);
        ymax = ymax.max(y);
    }
    // No finite data at all: fall back to a fixed 0..1 range so the frame/axis text still
    // renders sensibly instead of showing +inf/-inf.
    if !xmin.is_finite() || !xmax.is_finite() {
        xmin = 0.0;
        xmax = 1.0;
    }
    if !ymin.is_finite() || !ymax.is_finite() {
        ymin = 0.0;
        ymax = 1.0;
    }
    // Degenerate (single-valued) axis: expand around the value so the scale below never
    // divides by zero.
    if xmax <= xmin {
        xmin -= 0.5;
        xmax += 0.5;
    }
    if ymax <= ymin {
        ymin -= 0.5;
        ymax += 0.5;
    }
    let xspan = xmax - xmin;
    let yspan = ymax - ymin;

    let mut canvas = make_canvas(style, width_cells, height_cells);
    let width_dots = canvas.width_dots();
    let height_dots = canvas.height_dots();
    if width_dots > 0 && height_dots > 0 {
        let max_dx = (width_dots - 1) as f64;
        let max_dy = (height_dots - 1) as f64;
        for &(x, y) in finite_points() {
            let dx = ((x - xmin) / xspan * max_dx).round().clamp(0.0, max_dx) as usize;
            // Row 0 is the TOP of the canvas, so larger y (up) must map to a SMALLER row.
            let dy = ((ymax - y) / yspan * max_dy).round().clamp(0.0, max_dy) as usize;
            canvas.set(dx, dy);
        }
    }

    let labels = series
        .iter()
        .map(|(label, _)| *label)
        .collect::<Vec<_>>()
        .join(", ");
    let top_text = format!(" {} \u{2014} y:[{:.2}, {:.2}] ", labels, ymin, ymax);
    let bottom_text = format!(" x:[{:.2}, {:.2}] ", xmin, xmax);

    let mut out = String::new();
    out.push_str(&frame_line('\u{250C}', '\u{2500}', '\u{2510}', &top_text, width_cells));
    out.push('\n');
    for line in canvas.render().lines() {
        out.push('\u{2502}');
        out.push_str(line);
        out.push('\u{2502}');
        out.push('\n');
    }
    out.push_str(&frame_line(
        '\u{2514}',
        '\u{2500}',
        '\u{2518}',
        &bottom_text,
        width_cells,
    ));
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- Canvas primitives -------------------------------------------------------------

    #[test]
    fn braille_blank_cell_is_u2800() {
        let c = BrailleCanvas::new(1, 1);
        assert_eq!(c.render(), "\u{2800}");
    }

    #[test]
    fn braille_single_dot_top_left() {
        let mut c = BrailleCanvas::new(1, 1);
        c.set(0, 0); // dot 1 -> bit 0x01
        assert_eq!(c.render(), "\u{2801}");
    }

    #[test]
    fn braille_single_dot_bottom_right() {
        let mut c = BrailleCanvas::new(1, 1);
        c.set(1, 3); // dot 8 -> bit 0x80
        assert_eq!(c.render(), "\u{2880}");
    }

    #[test]
    fn braille_full_cell_is_all_eight_dots() {
        let mut c = BrailleCanvas::new(1, 1);
        for y in 0..4 {
            for x in 0..2 {
                c.set(x, y);
            }
        }
        // 0xFF -> U+28FF, the last codepoint in the braille block.
        assert_eq!(c.render(), "\u{28FF}");
    }

    #[test]
    fn braille_two_cell_row_addresses_columns_independently() {
        let mut c = BrailleCanvas::new(2, 1);
        c.set(0, 0); // left cell, dot 1
        c.set(2, 0); // right cell (x=2 is dot-column 0 of cell 1), dot 1
        assert_eq!(c.render(), "\u{2801}\u{2801}");
    }

    #[test]
    fn braille_out_of_range_set_is_ignored_not_panicking() {
        let mut c = BrailleCanvas::new(1, 1);
        c.set(99, 99);
        c.set(2, 0);
        c.set(0, 4);
        assert_eq!(c.render(), "\u{2800}");
    }

    #[test]
    fn braille_multirow_multicol_layout() {
        let mut c = BrailleCanvas::new(2, 2);
        c.set(0, 0); // cell (0,0), dot 1
        c.set(3, 7); // cell (1,1), dot 8 (x=3 -> col 1 of cell 1; y=7 -> row 3 of cell 1)
        let rendered = c.render();
        let lines: Vec<&str> = rendered.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "\u{2801}\u{2800}");
        assert_eq!(lines[1], "\u{2800}\u{2880}");
    }

    #[test]
    fn ascii_blank_cell_is_space() {
        let c = AsciiCanvas::new(1, 1);
        assert_eq!(c.render(), " ");
    }

    #[test]
    fn ascii_any_dot_in_cell_renders_star() {
        let mut c = AsciiCanvas::new(1, 1);
        c.set(1, 2); // a single dot anywhere in the cell is enough
        assert_eq!(c.render(), "*");
    }

    #[test]
    fn ascii_two_cell_row() {
        let mut c = AsciiCanvas::new(3, 1);
        c.set(0, 0);
        c.set(5, 3); // cell 2 (x=5 -> col 5/2=2)
        assert_eq!(c.render(), "* *");
    }

    #[test]
    fn ascii_out_of_range_set_is_ignored() {
        let mut c = AsciiCanvas::new(1, 1);
        c.set(100, 100);
        assert_eq!(c.render(), " ");
    }

    // ---- render_chart: framing, axis text, autoscale -----------------------------------

    #[test]
    fn render_chart_frame_dimensions_and_borders() {
        let pts = [(0.0, 0.0), (1.0, 1.0)];
        let chart = render_chart(&[("s", &pts)], 6, 2, CanvasStyle::Ascii);
        let lines: Vec<&str> = chart.lines().collect();
        // top border + height_cells canvas rows + bottom border
        assert_eq!(lines.len(), 2 + 2);
        for line in &lines {
            assert_eq!(line.chars().count(), 6 + 2);
        }
        assert!(lines[0].starts_with('\u{250C}') && lines[0].ends_with('\u{2510}'));
        assert!(lines[3].starts_with('\u{2514}') && lines[3].ends_with('\u{2518}'));
        for row in &lines[1..3] {
            assert!(row.starts_with('\u{2502}') && row.ends_with('\u{2502}'));
        }
    }

    #[test]
    fn render_chart_exact_text_two_point_diagonal_ascii() {
        // Two points at the extreme corners of a canvas wide enough that the axis-range
        // border text isn't truncated (22 and 17 characters respectively — see the
        // hand-counted fill lengths below). width_cells=24 -> 48 dot columns,
        // height_cells=1 -> 4 dot rows.
        let pts: [(f64, f64); 2] = [(0.0, 0.0), (10.0, 100.0)];
        let chart = render_chart(&[("s", &pts)], 24, 1, CanvasStyle::Ascii);
        // (0.0, 0.0) is the min on both axes -> dot column 0, dot row 3 (bottom).
        // (10.0, 100.0) is the max on both axes -> dot column 47 (last), dot row 0 (top).
        // Both dots fall in cell_y=0 (the only row), at cell_x=0 and cell_x=23 (the first
        // and last of the 24 cells) respectively, so the row is '*' + 22 spaces + '*'.
        let expected = format!(
            "\u{250C} s \u{2014} y:[0.00, 100.00] {}\u{2510}\n\
             \u{2502}*{}*\u{2502}\n\
             \u{2514} x:[0.00, 10.00] {}\u{2518}",
            "\u{2500}".repeat(2),
            " ".repeat(22),
            "\u{2500}".repeat(7)
        );
        assert_eq!(chart, expected);
    }

    #[test]
    fn render_chart_exact_text_single_point_braille() {
        // One point, one cell, width_cells=1 (so the border text is truncated down to a
        // single character — exercised deliberately here; long_label truncation is covered
        // separately below). The point itself lands dead-center of the 2x4 dot cell: a
        // single point always makes both axes degenerate (min==max), which
        // render_chart expands to +/-0.5 around the value, so the point is always
        // EXACTLY at the midpoint of both the dot-column and dot-row ranges.
        // dx = (5.0 - 4.5) / 1.0 * (2-1) = 0.5 -> f64::round rounds half away from zero -> 1
        // dy = (5.5 - 5.0) / 1.0 * (4-1) = 1.5 -> rounds to 2
        // set(1, 2) -> cell (0,0), dot bit BRAILLE_BITS[2][1] = 0x20 -> U+2820.
        let pts: [(f64, f64); 1] = [(5.0, 5.0)];
        let chart = render_chart(&[("p", &pts)], 1, 1, CanvasStyle::Braille);
        let expected = "\u{250C} \u{2510}\n\
                         \u{2502}\u{2820}\u{2502}\n\
                         \u{2514} \u{2518}";
        assert_eq!(chart, expected);
    }

    #[test]
    fn render_chart_multiple_series_labels_joined_in_border() {
        let a = [(0.0, 0.0)];
        let b = [(1.0, 1.0)];
        let chart = render_chart(&[("alpha", &a), ("beta", &b)], 20, 2, CanvasStyle::Ascii);
        let top = chart.lines().next().unwrap();
        assert!(top.contains("alpha, beta"));
    }

    #[test]
    fn render_chart_empty_series_still_frames_with_placeholder_range() {
        // width_cells=24 is wide enough that the 23-character axis-range border text
        // below isn't truncated (truncation itself is covered separately by
        // render_chart_long_label_is_truncated_not_panicking).
        let empty: [(f64, f64); 0] = [];
        let chart = render_chart(&[("none", &empty)], 24, 1, CanvasStyle::Ascii);
        let lines: Vec<&str> = chart.lines().collect();
        assert_eq!(lines.len(), 3);
        assert!(lines[0].contains("y:[0.00, 1.00]"));
        assert!(lines[2].contains("x:[0.00, 1.00]"));
        // No dots set anywhere.
        assert_eq!(lines[1], format!("\u{2502}{}\u{2502}", " ".repeat(24)));
    }

    #[test]
    fn render_chart_ignores_non_finite_points() {
        // width_cells=24 for the same untruncated-border-text reason as above.
        let pts: [(f64, f64); 3] = [(0.0, 0.0), (f64::NAN, 1.0), (f64::INFINITY, 2.0)];
        let chart = render_chart(&[("s", &pts)], 24, 1, CanvasStyle::Ascii);
        // Only the finite (0.0, 0.0) point survives, so both axes collapse to the
        // degenerate +/-0.5 expansion around 0.0.
        assert!(chart.contains("y:[-0.50, 0.50]"));
        assert!(chart.contains("x:[-0.50, 0.50]"));
    }

    #[test]
    fn render_chart_long_label_is_truncated_not_panicking() {
        let pts = [(0.0, 0.0)];
        let label = "x".repeat(100);
        let chart = render_chart(&[(label.as_str(), &pts)], 5, 1, CanvasStyle::Ascii);
        let top = chart.lines().next().unwrap();
        assert_eq!(top.chars().count(), 5 + 2);
    }
}
