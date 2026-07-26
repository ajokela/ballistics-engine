use crate::transonic_drag::{get_projectile_shape, transonic_correction, ProjectileShape};
use crate::DragModel;
use ndarray::ArrayD;
use std::sync::LazyLock;
/// Drag coefficient calculations for ballistics using actual drag table data
use std::path::Path;

/// Drag table data structure
#[derive(Debug, Clone)]
pub struct DragTable {
    pub mach_values: Vec<f64>,
    pub cd_values: Vec<f64>,
}

impl DragTable {
    /// Create a new drag table from mach and cd arrays
    pub fn new(mach_values: Vec<f64>, cd_values: Vec<f64>) -> Self {
        Self {
            mach_values,
            cd_values,
        }
    }

    /// Validated constructor for user-supplied drag decks. Enforces: equal-length axes,
    /// at least 2 points, strictly-ascending finite non-negative Mach, and finite positive Cd.
    /// Returns a descriptive, 1-based-row error string on failure (never panics).
    pub fn try_new(mach_values: Vec<f64>, cd_values: Vec<f64>) -> Result<Self, String> {
        if mach_values.len() != cd_values.len() {
            return Err(format!(
                "drag table has {} Mach values but {} Cd values; the columns must be equal length",
                mach_values.len(),
                cd_values.len()
            ));
        }
        if mach_values.len() < 2 {
            return Err(format!(
                "drag table needs at least 2 points, got {}",
                mach_values.len()
            ));
        }
        for (i, &m) in mach_values.iter().enumerate() {
            if !m.is_finite() || m < 0.0 {
                return Err(format!(
                    "drag table Mach at row {} must be finite and >= 0, got {m}",
                    i + 1
                ));
            }
            if i > 0 && m <= mach_values[i - 1] {
                return Err(format!(
                    "drag table Mach must strictly ascend; row {} ({m}) <= row {} ({})",
                    i + 1,
                    i,
                    mach_values[i - 1]
                ));
            }
        }
        for (i, &cd) in cd_values.iter().enumerate() {
            if !cd.is_finite() || cd <= 0.0 {
                return Err(format!(
                    "drag table Cd at row {} must be finite and > 0, got {cd}",
                    i + 1
                ));
            }
        }
        Ok(Self { mach_values, cd_values })
    }

    /// Parse a user drag deck from CSV text: two columns `mach,cd` per line. Blank lines and
    /// lines starting with `#` are ignored; a single leading header row is skipped once, but only
    /// when its first column is not itself a valid number (e.g. `mach,cd`) — a first row whose
    /// first column does parse as a float (e.g. `0.5` or `0.5,O.2`) is data, not a header, so a
    /// missing/invalid second column there is a hard error, not a silent skip. Any unparseable
    /// row is a hard error citing its 1-based line number. Values are validated via `try_new`.
    pub fn from_csv_str(csv: &str) -> Result<Self, String> {
        let mut mach_values = Vec::new();
        let mut cd_values = Vec::new();
        let mut header_skipped = false;
        for (lineno, raw) in csv.lines().enumerate() {
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut cols = line.split(',');
            let m = cols.next().map(str::trim);
            let cd = cols.next().map(str::trim);
            let m_parsed = m.and_then(|s| s.parse::<f64>().ok());
            match (m_parsed, cd.and_then(|s| s.parse::<f64>().ok())) {
                (Some(m), Some(cd)) => {
                    mach_values.push(m);
                    cd_values.push(cd);
                }
                _ => {
                    if !header_skipped && mach_values.is_empty() && m_parsed.is_none() {
                        // Tolerate one leading header row (e.g. "mach,cd") — but only when its
                        // first column gives no numeric evidence of being a data row. A row whose
                        // first column *does* parse (e.g. "0.5" or "0.5,O.2") is malformed data,
                        // not a header, and must error rather than be silently discarded.
                        header_skipped = true;
                        continue;
                    }
                    return Err(format!(
                        "drag table CSV: could not parse two numbers from line {}: {:?}",
                        lineno + 1,
                        raw
                    ));
                }
            }
        }
        if mach_values.is_empty() {
            return Err("drag table CSV contained no data rows".to_string());
        }
        Self::try_new(mach_values, cd_values)
    }

    /// Load and validate a user drag deck from a CSV file path.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self, String> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .map_err(|e| format!("could not read drag table {}: {e}", path.display()))?;
        Self::from_csv_str(&text)
    }

    /// Interpolate drag coefficient for a Mach number, holding the nearest tabulated endpoint
    /// outside the table's measured domain.
    pub fn interpolate(&self, mach: f64) -> f64 {
        let n = self.mach_values.len();

        if n == 0 {
            return 0.5; // Fallback
        }

        if n == 1 {
            return self.cd_values.first().copied().unwrap_or(0.5);
        }

        // A table has no information beyond its measured Mach domain. Hold the nearest endpoint
        // rather than extending the local edge slope indefinitely (which can drive Cd to 0.01).
        if mach <= self.mach_values[0] {
            return self.cd_values.first().copied().unwrap_or(0.5);
        }

        if mach >= self.mach_values[n - 1] {
            // Guard against a caller-built mismatched table (`new` is infallible): index the Cd
            // axis defensively rather than trusting the Mach-derived length.
            return self.cd_values.get(n - 1).copied()
                .or_else(|| self.cd_values.last().copied())
                .unwrap_or(0.5);
        }

        // Find the segment containing the mach value. Binary search over the
        // strictly-ascending mach axis; bit-identical to the previous linear scan
        // (first segment [i, i+1] with m[i] <= mach <= m[i+1]) but O(log n).
        let idx = self
            .mach_values
            .partition_point(|&m| m < mach)
            .saturating_sub(1)
            .min(n - 2);

        // Use cubic interpolation if we have enough points, otherwise linear
        if idx > 0 && idx < n - 2 {
            // Cubic interpolation using 4 points
            self.cubic_interpolate(mach, idx)
        } else {
            // Linear interpolation for edge cases
            self.linear_interpolate(mach, idx)
        }
    }

    /// Linear interpolation between two points
    pub fn linear_interpolate(&self, mach: f64, idx: usize) -> f64 {
        // Bounds check
        if idx + 1 >= self.mach_values.len() || idx + 1 >= self.cd_values.len() {
            return self.cd_values.get(idx).copied().unwrap_or(0.5);
        }

        let x0 = self.mach_values[idx];
        let x1 = self.mach_values[idx + 1];
        let y0 = self.cd_values[idx];
        let y1 = self.cd_values[idx + 1];

        if (x1 - x0).abs() < crate::constants::MIN_DIVISION_THRESHOLD {
            return y0;
        }

        let t = (mach - x0) / (x1 - x0);
        y0 + t * (y1 - y0)
    }

    /// Cubic Hermite interpolation using four points and centered chord-slope tangents.
    pub fn cubic_interpolate(&self, mach: f64, idx: usize) -> f64 {
        // Ensure we have enough points for cubic interpolation
        if idx == 0 || idx + 1 >= self.mach_values.len() || idx + 1 >= self.cd_values.len() {
            // Fall back to linear interpolation if not enough points
            return self.linear_interpolate(mach, idx);
        }

        // Use points at idx-1, idx, idx+1, idx+2
        let x = [
            self.mach_values[idx - 1],
            self.mach_values[idx],
            self.mach_values[idx + 1],
            if idx + 2 < self.mach_values.len() {
                self.mach_values[idx + 2]
            } else {
                self.mach_values[idx + 1]
            },
        ];
        let y = [
            self.cd_values[idx - 1],
            self.cd_values[idx],
            self.cd_values[idx + 1],
            if idx + 2 < self.cd_values.len() {
                self.cd_values[idx + 2]
            } else {
                self.cd_values[idx + 1]
            },
        ];

        // Scale centered chord-slope tangents by this segment's actual width. This Hermite
        // construction remains C1 across non-uniform knots; using the fixed uniform Catmull-Rom
        // coefficient matrix here bends even affine data when adjacent Mach intervals differ.
        let segment_width = x[2] - x[1];
        let left_chord_width = x[2] - x[0];
        let right_chord_width = x[3] - x[1];
        if segment_width.abs() < crate::constants::MIN_DIVISION_THRESHOLD
            || left_chord_width.abs() < crate::constants::MIN_DIVISION_THRESHOLD
            || right_chord_width.abs() < crate::constants::MIN_DIVISION_THRESHOLD
        {
            return self.linear_interpolate(mach, idx);
        }
        let t = (mach - x[1]) / segment_width;
        let t2 = t * t;
        let t3 = t2 * t;

        let tangent1 = segment_width * (y[2] - y[0]) / left_chord_width;
        let tangent2 = segment_width * (y[3] - y[1]) / right_chord_width;
        let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
        let h10 = t3 - 2.0 * t2 + t;
        let h01 = -2.0 * t3 + 3.0 * t2;
        let h11 = t3 - t2;

        h00 * y[1] + h10 * tangent1 + h01 * y[2] + h11 * tangent2
    }
}

/// Load drag table from NumPy binary file or CSV fallback
pub fn load_drag_table(
    drag_tables_dir: &Path,
    filename: &str,
    fallback_data: &[(f64, f64)],
) -> DragTable {
    // Try to load NumPy binary file first
    let npy_path = drag_tables_dir.join(format!("{filename}.npy"));
    if let Ok(array) = ndarray_npy::read_npy::<_, ArrayD<f64>>(&npy_path) {
        if let Ok(array_2d) = array.into_dimensionality::<ndarray::Ix2>() {
            let mach_values: Vec<f64> = array_2d.column(0).to_vec();
            let cd_values: Vec<f64> = array_2d.column(1).to_vec();
            return DragTable::new(mach_values, cd_values);
        }
    }

    // Fallback to CSV file. Hand-parsed (MBA-1331: the `csv` crate is now a
    // cli-feature dependency, and this was its only lib call site): any line whose
    // first two comma-separated fields parse as f64 is a data row, everything else
    // (headers, blanks, comments) is skipped. NOTE the old csv::Reader consumed the
    // FIRST row as headers unconditionally, so a headerless file silently lost its
    // first data point — parse-based skipping keeps that point.
    let csv_path = drag_tables_dir.join(format!("{filename}.csv"));
    if let Ok(bytes) = std::fs::read(&csv_path) {
        // Lossy UTF-8 (a stray Latin-1 byte in a header comment must not reject the
        // whole file) and bare-CR tolerance — the old csv::Reader accepted both.
        let text = String::from_utf8_lossy(&bytes).replace('\r', "\n");
        let mut mach_values = Vec::new();
        let mut cd_values = Vec::new();

        for line in text.lines() {
            let mut fields = line.split(',');
            if let (Some(m_str), Some(cd_str)) = (fields.next(), fields.next()) {
                // trim_matches('"'): the old csv::Reader unquoted "1.05","0.42"-style
                // rows (Excel exports); keep accepting them.
                if let (Ok(mach), Ok(cd)) = (
                    m_str.trim().trim_matches('"').trim().parse::<f64>(),
                    cd_str.trim().trim_matches('"').trim().parse::<f64>(),
                ) {
                    mach_values.push(mach);
                    cd_values.push(cd);
                }
            }
        }

        if !mach_values.is_empty() {
            return DragTable::new(mach_values, cd_values);
        }
    }

    // Use fallback data if both file loading methods fail
    let mach_values: Vec<f64> = fallback_data.iter().map(|(m, _)| *m).collect();
    let cd_values: Vec<f64> = fallback_data.iter().map(|(_, cd)| *cd).collect();
    DragTable::new(mach_values, cd_values)
}

/// Find the drag tables directory relative to the current location
fn find_drag_tables_dir() -> Option<std::path::PathBuf> {
    // Try common relative paths from the Rust crate location
    let candidates = [
        "../drag_tables",
        "../../drag_tables",
        "../../../drag_tables",
        "drag_tables",
    ];

    for candidate in &candidates {
        let path = Path::new(candidate);
        if path.exists() && path.is_dir() {
            return Some(path.to_path_buf());
        }
    }

    None
}

/// Parse an embedded CSV drag table (`mach,cd` per line, header tolerated). Used to bake the
/// high-resolution G1/G7 tables (data/*.csv) into the binary so the engine never depends on a
/// runtime `drag_tables/` directory existing. Falls back to the supplied coarse table only if
/// parsing yields no points (the shipped data files always parse).
fn parse_embedded_drag_table(csv: &str, fallback: &[(f64, f64)]) -> DragTable {
    let mut mach_values = Vec::new();
    let mut cd_values = Vec::new();
    for line in csv.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let mut cols = line.split(',');
        if let (Some(m), Some(cd)) = (cols.next(), cols.next()) {
            if let (Ok(m), Ok(cd)) = (m.trim().parse::<f64>(), cd.trim().parse::<f64>()) {
                mach_values.push(m);
                cd_values.push(cd);
            }
        }
    }
    if mach_values.is_empty() {
        mach_values = fallback.iter().map(|(m, _)| *m).collect();
        cd_values = fallback.iter().map(|(_, cd)| *cd).collect();
    }
    DragTable::new(mach_values, cd_values)
}

/// G1 drag table — high-resolution data baked in from data/g1.csv at compile time (MBA-939).
/// The previous runtime loader searched for a `drag_tables/` directory that does not exist when
/// the binary runs (the tables ship under data/), so the engine silently used the coarse 21-point
/// fallback below, flattening the transonic drag rise. include_str! guarantees the full table.
static G1_DRAG_TABLE: LazyLock<DragTable> = LazyLock::new(|| {
    // Coarse 21-point fallback, retained only for the impossible parse-failure path.
    let fallback_data = [
        (0.0, 0.2629),
        (0.5, 0.2695),
        (0.6, 0.2752),
        (0.7, 0.2817),
        (0.8, 0.2902),
        (0.9, 0.3012),
        (1.0, 0.4805),
        (1.1, 0.5933),
        (1.2, 0.6318),
        (1.3, 0.6440),
        (1.4, 0.6444),
        (1.5, 0.6372),
        (1.6, 0.6252),
        (1.7, 0.6105),
        (1.8, 0.5956),
        (1.9, 0.5815),
        (2.0, 0.5934),
        (2.5, 0.5598),
        (3.0, 0.5133),
        (4.0, 0.4811),
        (5.0, 0.4988),
    ];

    parse_embedded_drag_table(include_str!("../data/g1.csv"), &fallback_data)
});

/// G7 drag table — high-resolution data baked in from data/g7.csv at compile time (MBA-939).
/// Same root cause as G1: the runtime `drag_tables/` loader never resolved, so the coarse
/// 21-point fallback was used, missing the Mach 0.9->1.0 transonic knee (the embedded 0.9 point
/// was even wrong: 0.1294 vs the true 0.1464). include_str! bakes in the full 84-point table.
static G7_DRAG_TABLE: LazyLock<DragTable> = LazyLock::new(|| {
    // Coarse 21-point fallback, retained only for the impossible parse-failure path.
    let fallback_data = [
        (0.0, 0.1198),
        (0.5, 0.1197),
        (0.6, 0.1202),
        (0.7, 0.1213),
        (0.8, 0.1240),
        (0.9, 0.1294),
        (1.0, 0.3803),
        (1.1, 0.4015),
        (1.2, 0.4043),
        (1.3, 0.3956),
        (1.4, 0.3814),
        (1.5, 0.3663),
        (1.6, 0.3520),
        (1.7, 0.3398),
        (1.8, 0.3297),
        (1.9, 0.3221),
        (2.0, 0.2980),
        (2.5, 0.2731),
        (3.0, 0.2424),
        (4.0, 0.2196),
        (5.0, 0.1618),
    ];

    parse_embedded_drag_table(include_str!("../data/g7.csv"), &fallback_data)
});

/// G6 drag table - flat-base with 6 caliber secant ogive (military FMJ bullets)
/// MBA-156: Added for completeness with ballistics_rust
static G6_DRAG_TABLE: LazyLock<DragTable> = LazyLock::new(|| {
    let fallback_data = [
        (0.0, 0.2617),
        (0.05, 0.2553),
        (0.10, 0.2491),
        (0.15, 0.2432),
        (0.20, 0.2376),
        (0.25, 0.2324),
        (0.30, 0.2278),
        (0.35, 0.2238),
        (0.40, 0.2205),
        (0.45, 0.2177),
        (0.50, 0.2155),
        (0.55, 0.2138),
        (0.60, 0.2126),
        (0.65, 0.2121),
        (0.70, 0.2122),
        (0.75, 0.2132),
        (0.80, 0.2154),
        (0.85, 0.2194),
        (0.875, 0.2229),
        (0.90, 0.2297),
        (0.925, 0.2449),
        (0.95, 0.2732),
        (0.975, 0.3141),
        (1.0, 0.3597),
        (1.025, 0.3994),
        (1.05, 0.4261),
        (1.075, 0.4402),
        (1.10, 0.4465),
        (1.125, 0.4490),
        (1.15, 0.4497),
        (1.175, 0.4494),
        (1.20, 0.4482),
        (1.225, 0.4464),
        (1.25, 0.4441),
        (1.30, 0.4390),
        (1.35, 0.4336),
        (1.40, 0.4279),
        (1.45, 0.4221),
        (1.50, 0.4162),
        (1.55, 0.4102),
        (1.60, 0.4042),
        (1.65, 0.3981),
        (1.70, 0.3919),
        (1.75, 0.3855),
        (1.80, 0.3788),
        (1.85, 0.3721),
        (1.90, 0.3652),
        (1.95, 0.3583),
        (2.0, 0.3515),
        (2.05, 0.3447),
        (2.10, 0.3381),
        (2.15, 0.3314),
        (2.20, 0.3249),
        (2.25, 0.3185),
        (2.30, 0.3122),
        (2.35, 0.3060),
        (2.40, 0.3000),
        (2.45, 0.2941),
        (2.50, 0.2883),
        (2.60, 0.2772),
        (2.70, 0.2668),
        (2.80, 0.2574),
        (2.90, 0.2487),
        (3.0, 0.2407),
        (3.10, 0.2333),
        (3.20, 0.2265),
        (3.30, 0.2202),
        (3.40, 0.2144),
        (3.50, 0.2089),
        (3.60, 0.2039),
        (3.70, 0.1991),
        (3.80, 0.1947),
        (3.90, 0.1905),
        (4.0, 0.1866),
        (4.20, 0.1794),
        (4.40, 0.1730),
        (4.60, 0.1673),
        (4.80, 0.1621),
        (5.0, 0.1574),
    ];

    if let Some(drag_dir) = find_drag_tables_dir() {
        load_drag_table(&drag_dir, "g6", &fallback_data)
    } else {
        // Use fallback data if directory not found
        let mach_values: Vec<f64> = fallback_data.iter().map(|(m, _)| *m).collect();
        let cd_values: Vec<f64> = fallback_data.iter().map(|(_, cd)| *cd).collect();
        DragTable::new(mach_values, cd_values)
    }
});

/// G8 drag table - flat-base with 10 caliber secant ogive
/// MBA-156: Added for completeness with ballistics_rust
static G8_DRAG_TABLE: LazyLock<DragTable> = LazyLock::new(|| {
    let fallback_data = [
        (0.0, 0.2105),
        (0.05, 0.2105),
        (0.10, 0.2104),
        (0.15, 0.2104),
        (0.20, 0.2103),
        (0.25, 0.2103),
        (0.30, 0.2103),
        (0.35, 0.2103),
        (0.40, 0.2103),
        (0.45, 0.2102),
        (0.50, 0.2102),
        (0.55, 0.2102),
        (0.60, 0.2102),
        (0.65, 0.2102),
        (0.70, 0.2103),
        (0.75, 0.2103),
        (0.80, 0.2104),
        (0.825, 0.2104),
        (0.85, 0.2105),
        (0.875, 0.2106),
        (0.90, 0.2109),
        (0.925, 0.2183),
        (0.95, 0.2571),
        (0.975, 0.3358),
        (1.0, 0.4068),
        (1.025, 0.4378),
        (1.05, 0.4476),
        (1.075, 0.4493),
        (1.10, 0.4477),
        (1.125, 0.4450),
        (1.15, 0.4419),
        (1.20, 0.4353),
        (1.25, 0.4283),
        (1.30, 0.4208),
        (1.35, 0.4133),
        (1.40, 0.4059),
        (1.45, 0.3986),
        (1.50, 0.3915),
        (1.55, 0.3845),
        (1.60, 0.3777),
        (1.65, 0.3710),
        (1.70, 0.3645),
        (1.75, 0.3581),
        (1.80, 0.3519),
        (1.85, 0.3458),
        (1.90, 0.3400),
        (1.95, 0.3343),
        (2.0, 0.3288),
        (2.05, 0.3234),
        (2.10, 0.3182),
        (2.15, 0.3131),
        (2.20, 0.3081),
        (2.25, 0.3032),
        (2.30, 0.2983),
        (2.35, 0.2937),
        (2.40, 0.2891),
        (2.45, 0.2845),
        (2.50, 0.2802),
        (2.60, 0.2720),
        (2.70, 0.2642),
        (2.80, 0.2569),
        (2.90, 0.2499),
        (3.0, 0.2432),
        (3.10, 0.2368),
        (3.20, 0.2308),
        (3.30, 0.2251),
        (3.40, 0.2197),
        (3.50, 0.2147),
        (3.60, 0.2101),
        (3.70, 0.2058),
        (3.80, 0.2019),
        (3.90, 0.1983),
        (4.0, 0.1950),
        (4.20, 0.1890),
        (4.40, 0.1837),
        (4.60, 0.1791),
        (4.80, 0.1750),
        (5.0, 0.1713),
    ];

    if let Some(drag_dir) = find_drag_tables_dir() {
        load_drag_table(&drag_dir, "g8", &fallback_data)
    } else {
        // Use fallback data if directory not found
        let mach_values: Vec<f64> = fallback_data.iter().map(|(m, _)| *m).collect();
        let cd_values: Vec<f64> = fallback_data.iter().map(|(_, cd)| *cd).collect();
        DragTable::new(mach_values, cd_values)
    }
});

/// G2 drag table — banded-based projectile (Aberdeen/BRL standard, as tabulated in McCoy,
/// Modern Exterior Ballistics). High-resolution data baked in from data/g2.csv (MBA-1386);
/// see that file's header for provenance.
static G2_DRAG_TABLE: LazyLock<DragTable> = LazyLock::new(|| {
    // 3-point fallback for the impossible parse-failure path only.
    let fallback_data = [(0.0, 0.2303), (1.0, 0.3983), (5.0, 0.1648)];
    parse_embedded_drag_table(include_str!("../data/g2.csv"), &fallback_data)
});

/// G5 drag table — short 7.5 caliber tangent ogive boat-tail (Aberdeen/BRL standard, as
/// tabulated in McCoy, Modern Exterior Ballistics). High-resolution data baked in from
/// data/g5.csv (MBA-1386); see that file's header for provenance.
static G5_DRAG_TABLE: LazyLock<DragTable> = LazyLock::new(|| {
    // 3-point fallback for the impossible parse-failure path only.
    let fallback_data = [(0.0, 0.1710), (1.0, 0.3379), (5.0, 0.2280)];
    parse_embedded_drag_table(include_str!("../data/g5.csv"), &fallback_data)
});

/// GI drag table — flat-based, 5.5 caliber tangent ogive (Aberdeen/BRL standard, as
/// tabulated in McCoy, Modern Exterior Ballistics). High-resolution data baked in from
/// data/gi.csv (MBA-1386); see that file's header for provenance.
static GI_DRAG_TABLE: LazyLock<DragTable> = LazyLock::new(|| {
    // 3-point fallback for the impossible parse-failure path only.
    let fallback_data = [(0.0, 0.2282), (1.0, 0.4349), (5.0, 0.4082)];
    parse_embedded_drag_table(include_str!("../data/gi.csv"), &fallback_data)
});

/// GS drag table — spherical (round-ball) projectile (Aberdeen/BRL standard, as tabulated in
/// McCoy, Modern Exterior Ballistics). High-resolution data baked in from data/gs.csv
/// (MBA-1386); see that file's header for provenance. Note the source table only extends to
/// Mach 4.0 (not 5.0 like the other families).
static GS_DRAG_TABLE: LazyLock<DragTable> = LazyLock::new(|| {
    // 3-point fallback for the impossible parse-failure path only.
    let fallback_data = [(0.0, 0.4662), (1.0, 0.8140), (4.0, 0.9280)];
    parse_embedded_drag_table(include_str!("../data/gs.csv"), &fallback_data)
});

/// RA4 drag table — British RA 1929 reference function (as tabulated in McCoy, Modern
/// Exterior Ballistics). High-resolution data baked in from data/ra4.csv (MBA-1386); see that
/// file's header for provenance. Note the source table only extends to Mach 4.0 (not 5.0 like
/// most of the G-family tables).
static RA4_DRAG_TABLE: LazyLock<DragTable> = LazyLock::new(|| {
    // 3-point fallback for the impossible parse-failure path only.
    let fallback_data = [(0.0, 0.2283), (1.0, 0.3975), (4.0, 0.4969)];
    parse_embedded_drag_table(include_str!("../data/ra4.csv"), &fallback_data)
});

/// Get drag coefficient for given Mach number and drag model.
///
/// Every `DragModel` variant now has a dedicated high-resolution reference table (MBA-1386):
/// G1/G6/G7/G8 (Aberdeen/BRL, MBA-939/156), G2/G5/GI/GS (Aberdeen/BRL, as tabulated in McCoy),
/// and RA4 (British RA 1929 reference function).
pub fn get_drag_coefficient(mach: f64, drag_model: &DragModel) -> f64 {
    match drag_model {
        DragModel::G1 => G1_DRAG_TABLE.interpolate(mach),
        DragModel::G2 => G2_DRAG_TABLE.interpolate(mach),
        DragModel::G5 => G5_DRAG_TABLE.interpolate(mach),
        DragModel::G6 => G6_DRAG_TABLE.interpolate(mach),
        DragModel::G7 => G7_DRAG_TABLE.interpolate(mach),
        DragModel::G8 => G8_DRAG_TABLE.interpolate(mach),
        DragModel::GI => GI_DRAG_TABLE.interpolate(mach),
        DragModel::GS => GS_DRAG_TABLE.interpolate(mach),
        DragModel::RA4 => RA4_DRAG_TABLE.interpolate(mach),
    }
}

/// The built-in reference drag table for one model, as tabulated data (MBA-1424).
///
/// Exposes the same `(Mach, Cd)` pairs [`get_drag_coefficient`] interpolates, so a consumer can
/// plot or audit the standard curves without re-vendoring the numbers into their own code — and
/// without their copy drifting from the engine's as tables are refined.
///
/// These are the public-domain Aberdeen/BRL reference functions as tabulated in McCoy, *Modern
/// Exterior Ballistics*, plus the British RA 1929 function; each `data/*.csv` carries its own
/// provenance header. Nothing here is vendor drag data.
///
/// The tables do not share a Mach domain: most run to Mach 5, while GS and RA4 stop at Mach 4.
/// Read the returned table's `mach_values` rather than assuming a range.
///
/// This returns the *reference* curve for a standard projectile. For the effective Cd a specific
/// bullet actually flew — form-factor scaled, with segmented-BC band steps — see
/// [`crate::TrajectorySolver::effective_drag_coefficient`].
pub fn reference_drag_table(drag_model: &DragModel) -> &'static DragTable {
    match drag_model {
        DragModel::G1 => &G1_DRAG_TABLE,
        DragModel::G2 => &G2_DRAG_TABLE,
        DragModel::G5 => &G5_DRAG_TABLE,
        DragModel::G6 => &G6_DRAG_TABLE,
        DragModel::G7 => &G7_DRAG_TABLE,
        DragModel::G8 => &G8_DRAG_TABLE,
        DragModel::GI => &GI_DRAG_TABLE,
        DragModel::GS => &GS_DRAG_TABLE,
        DragModel::RA4 => &RA4_DRAG_TABLE,
    }
}

/// Get a standard G-table drag coefficient without double-counting transonic drag.
///
/// Standard G tables are total-drag curves that already contain the transonic
/// rise and wave drag. `apply_transonic_correction` and the shape inputs remain
/// in this public API for compatibility, but enabling the option does not stack
/// the separate empirical rise/wave model on top of a G-table coefficient.
pub fn get_drag_coefficient_with_transonic(
    mach: f64,
    drag_model: &DragModel,
    apply_transonic_correction: bool,
    projectile_shape: Option<ProjectileShape>,
    caliber: Option<f64>,
    weight_grains: Option<f64>,
) -> f64 {
    // Get base drag coefficient
    let base_cd = get_drag_coefficient(mach, drag_model);

    // Apply transonic corrections if requested and in transonic regime
    if apply_transonic_correction && (0.8..=1.3).contains(&mach) {
        // Determine projectile shape if not provided
        let shape = match projectile_shape {
            Some(s) => s,
            None => {
                if let (Some(cal), Some(weight)) = (caliber, weight_grains) {
                    get_projectile_shape(
                        cal,
                        weight,
                        match drag_model {
                            DragModel::G1 => "G1",
                            DragModel::G6 => "G6",
                            DragModel::G7 => "G7",
                            DragModel::G8 => "G8",
                            _ => "G1", // Default to G1
                        },
                    )
                } else {
                    ProjectileShape::Spitzer // Default
                }
            }
        };

        // Standard G-model tables are total-drag reference curves and already
        // contain their transonic rise and wave drag. Retain the public option
        // for API compatibility, but do not stack the empirical rise/wave model
        // on top of table Cd (MBA-1155).
        transonic_correction(mach, base_cd, shape, false)
    } else {
        base_cd
    }
}

/// Get drag coefficient with optional Reynolds correction.
///
/// The transonic option is retained for compatibility but, as documented by
/// [`get_drag_coefficient_with_transonic`], standard G tables are not corrected
/// a second time. Likewise, the Reynolds option only affects genuinely low-Re
/// (`Re < 10,000`) inputs; ordinary ballistic Reynolds numbers use the standard
/// table coefficient unchanged.
#[allow(clippy::too_many_arguments)] // Public compatibility API; grouping would be breaking.
pub fn get_drag_coefficient_full(
    mach: f64,
    drag_model: &DragModel,
    apply_transonic_correction: bool,
    apply_reynolds_correction: bool,
    projectile_shape: Option<ProjectileShape>,
    caliber: Option<f64>,
    weight_grains: Option<f64>,
    velocity_mps: Option<f64>,
    air_density_kg_m3: Option<f64>,
    temperature_c: Option<f64>,
) -> f64 {
    // Get base drag coefficient with transonic corrections if applicable
    let mut cd = get_drag_coefficient_with_transonic(
        mach,
        drag_model,
        apply_transonic_correction,
        projectile_shape,
        caliber,
        weight_grains,
    );

    // Route the opt-in low-Re helper for subsonic inputs. It leaves the ordinary
    // standard-table Reynolds-number range unchanged.
    if apply_reynolds_correction && mach < 1.0 {
        if let (Some(v), Some(cal), Some(rho), Some(temp)) =
            (velocity_mps, caliber, air_density_kg_m3, temperature_c)
        {
            use crate::reynolds::apply_reynolds_correction;
            cd = apply_reynolds_correction(cd, v, cal, rho, temp, mach);
        }
    }

    cd
}

#[cfg(test)]
#[allow(clippy::items_after_test_module)] // Keep the legacy public helper below in place.
mod tests {
    use super::*;

    #[test]
    fn test_g1_drag_coefficient_interpolation() {
        let cd = get_drag_coefficient(1.0, &DragModel::G1);
        // Should be close to the G1 standard value at Mach 1.0
        assert!(cd > 0.4 && cd < 0.6, "G1 CD at Mach 1.0: {cd}");
    }

    #[test]
    fn test_g7_drag_coefficient_interpolation() {
        let cd = get_drag_coefficient(1.0, &DragModel::G7);
        // Should be close to the G7 standard value at Mach 1.0
        assert!(cd > 0.3 && cd < 0.5, "G7 CD at Mach 1.0: {cd}");
    }

    #[test]
    fn standard_g_table_transonic_option_does_not_double_count_drag_rise() {
        let models = [
            DragModel::G1,
            DragModel::G2,
            DragModel::G5,
            DragModel::G6,
            DragModel::G7,
            DragModel::G8,
            DragModel::GI,
            DragModel::GS,
        ];
        for drag_model in models {
            for mach in [0.8, 0.95, 1.0, 1.1, 1.3] {
                let base_cd = get_drag_coefficient(mach, &drag_model);
                let corrected_cd = get_drag_coefficient_with_transonic(
                    mach,
                    &drag_model,
                    true,
                    Some(ProjectileShape::BoatTail),
                    Some(0.308),
                    Some(175.0),
                );
                assert_eq!(
                    corrected_cd.to_bits(),
                    base_cd.to_bits(),
                    "standard {drag_model:?} table already includes transonic drag at Mach \
                     {mach}: base={base_cd}, corrected={corrected_cd}"
                );

                let full_cd = get_drag_coefficient_full(
                    mach,
                    &drag_model,
                    true,
                    false,
                    Some(ProjectileShape::BoatTail),
                    Some(0.308),
                    Some(175.0),
                    None,
                    None,
                    None,
                );
                assert_eq!(full_cd.to_bits(), base_cd.to_bits());
            }
        }
    }

    #[test]
    fn test_drag_coefficient_continuity() {
        // Test that drag coefficient function is smooth
        for mach in [0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0] {
            let cd_before = get_drag_coefficient(mach - 0.01, &DragModel::G1);
            let cd_after = get_drag_coefficient(mach + 0.01, &DragModel::G1);
            let difference = (cd_after - cd_before).abs();
            assert!(
                difference < 0.05,
                "Large discontinuity at Mach {mach}: {cd_before} vs {cd_after}"
            );
        }
    }

    #[test]
    fn test_endpoint_bounds() {
        // Test endpoint hold below range
        let cd_low = get_drag_coefficient(0.0, &DragModel::G1);
        assert!(cd_low > 0.01 && cd_low < 0.5, "Low Mach G1: {cd_low}");

        // Test endpoint hold above range
        let cd_high = get_drag_coefficient(10.0, &DragModel::G1);
        assert!(cd_high > 0.01, "High Mach G1 should be positive: {cd_high}");

        // Same for G7
        let cd_low_g7 = get_drag_coefficient(0.0, &DragModel::G7);
        assert!(
            cd_low_g7 > 0.01,
            "Low Mach G7 should be positive: {cd_low_g7}"
        );

        let cd_high_g7 = get_drag_coefficient(20.0, &DragModel::G7);
        assert!(
            cd_high_g7 >= 0.01,
            "High Mach G7 should be positive: {cd_high_g7}"
        );
    }

    #[test]
    fn test_drag_table_creation() {
        let mach_vals = vec![0.5, 1.0, 1.5, 2.0];
        let cd_vals = vec![0.2, 0.5, 0.4, 0.3];
        let table = DragTable::new(mach_vals, cd_vals);

        // Test exact interpolation
        assert!((table.interpolate(1.0) - 0.5).abs() < 1e-10);

        // Test interpolation between points
        let cd_interp = table.interpolate(1.25);
        assert!(cd_interp > 0.4 && cd_interp < 0.5);
    }

    #[test]
    fn test_drag_table_empty() {
        let table = DragTable::new(vec![], vec![]);
        let result = table.interpolate(1.0);
        assert_eq!(result, 0.5); // Should return fallback value
    }

    #[test]
    fn test_drag_table_single_point() {
        let table = DragTable::new(vec![1.0], vec![0.4]);

        // Should return the single value for any Mach
        assert_eq!(table.interpolate(0.5), 0.4);
        assert_eq!(table.interpolate(1.0), 0.4);
        assert_eq!(table.interpolate(2.0), 0.4);
    }

    #[test]
    fn test_drag_table_two_points() {
        let table = DragTable::new(vec![1.0, 2.0], vec![0.4, 0.6]);

        // Exact matches
        assert!((table.interpolate(1.0) - 0.4).abs() < 1e-10);
        assert!((table.interpolate(2.0) - 0.6).abs() < 1e-10);

        // Linear interpolation
        let mid = table.interpolate(1.5);
        assert!((mid - 0.5).abs() < 1e-10);

        // Out-of-range values hold the nearest endpoint.
        let below = table.interpolate(0.5);
        assert_eq!(below.to_bits(), 0.4_f64.to_bits());

        let above = table.interpolate(3.0);
        assert_eq!(above.to_bits(), 0.6_f64.to_bits());
    }

    #[test]
    fn out_of_range_mach_holds_boundary_cd() {
        let table = DragTable::new(vec![0.5, 1.0, 2.0], vec![0.2, 0.5, 0.3]);

        for mach in [f64::NEG_INFINITY, -10.0, 0.49, 0.5] {
            assert_eq!(
                table.interpolate(mach).to_bits(),
                0.2_f64.to_bits(),
                "Mach {mach} must hold the first tabulated Cd"
            );
        }
        for mach in [2.0, 2.01, 100.0, f64::INFINITY] {
            assert_eq!(
                table.interpolate(mach).to_bits(),
                0.3_f64.to_bits(),
                "Mach {mach} must hold the last tabulated Cd"
            );
        }
    }

    #[test]
    fn test_linear_interpolation() {
        let table = DragTable::new(vec![0.0, 1.0, 2.0], vec![0.2, 0.5, 0.3]);

        // Test linear interpolation between first two points
        let result = table.linear_interpolate(0.5, 0);
        assert!((result - 0.35).abs() < 1e-10);

        // Test edge case with zero denominator
        let table_same = DragTable::new(vec![1.0, 1.0], vec![0.4, 0.6]);
        let result_same = table_same.linear_interpolate(1.0, 0);
        assert_eq!(result_same, 0.4); // Should return first value
    }

    #[test]
    fn test_cubic_interpolation() {
        // Create a table with enough points for cubic interpolation
        let table = DragTable::new(vec![0.5, 1.0, 1.5, 2.0, 2.5], vec![0.2, 0.4, 0.6, 0.5, 0.3]);

        // Test cubic interpolation in the middle
        let result = table.cubic_interpolate(1.25, 1);

        // Should be between the neighboring values
        assert!(result > 0.3 && result < 0.7);

        // Should be smooth (not exactly linear)
        let linear_result = table.linear_interpolate(1.25, 1);
        // Cubic and linear should be close but not identical for smooth curves
        assert!((result - linear_result).abs() < 0.2);
    }

    #[test]
    fn nonuniform_cubic_reproduces_affine_data() {
        let table = DragTable::new(
            vec![0.0, 1.0, 3.0, 4.0],
            vec![0.25, 0.3125, 0.4375, 0.5],
        );

        for mach in [1.5, 2.5] {
            let expected = 0.25 + mach / 16.0;
            let actual = table.interpolate(mach);
            assert_eq!(
                actual.to_bits(),
                expected.to_bits(),
                "non-uniform cubic bent affine data at Mach {mach}: {actual} vs {expected}"
            );
        }
    }

    #[test]
    fn nonuniform_cubic_is_c1_at_spacing_transition() {
        let table = DragTable::new(
            vec![0.0, 1.0, 3.0, 4.0, 7.0],
            vec![0.25, 0.265625, 0.390625, 0.5, 1.015625],
        );
        let knot = 3.0;
        let expected_at_knot = 0.390625_f64;
        let epsilon = 1e-6;
        let at_knot = table.interpolate(knot);
        let left_slope = (at_knot - table.interpolate(knot - epsilon)) / epsilon;
        let right_slope = (table.interpolate(knot + epsilon) - at_knot) / epsilon;

        assert_eq!(at_knot.to_bits(), expected_at_knot.to_bits());
        assert!(
            (left_slope - right_slope).abs() < 1e-5,
            "non-uniform cubic has a derivative kink: left={left_slope}, right={right_slope}"
        );
    }

    #[test]
    fn test_find_drag_tables_dir() {
        // This test may pass or fail depending on the environment
        // but should not panic
        let _dir = find_drag_tables_dir();
        // Just ensure the function doesn't panic
    }

    #[test]
    fn test_load_drag_table_fallback() {
        use std::path::Path;

        // Test with non-existent directory - should use fallback data
        let fake_dir = Path::new("/non/existent/directory");
        let fallback_data = [(0.5, 0.2), (1.0, 0.4), (1.5, 0.3)];

        let table = load_drag_table(fake_dir, "test", &fallback_data);

        // Should have fallback data
        assert_eq!(table.mach_values.len(), 3);
        assert_eq!(table.cd_values.len(), 3);
        assert_eq!(table.mach_values[0], 0.5);
        assert_eq!(table.cd_values[0], 0.2);
    }

    #[test]
    fn test_known_drag_values() {
        // Test against known ballistic standard values

        // G1 at Mach 1.0 should be around 0.4805
        let g1_mach1 = get_drag_coefficient(1.0, &DragModel::G1);
        assert!(
            (g1_mach1 - 0.4805).abs() < 0.01,
            "G1 at Mach 1.0: {g1_mach1}"
        );

        // G7 at Mach 1.0 should be around 0.3803
        let g7_mach1 = get_drag_coefficient(1.0, &DragModel::G7);
        assert!(
            (g7_mach1 - 0.3803).abs() < 0.01,
            "G7 at Mach 1.0: {g7_mach1}"
        );

        // G1 should generally be higher than G7 in transonic region
        assert!(g1_mach1 > g7_mach1, "G1 should be > G7 at Mach 1.0");
    }

    #[test]
    fn test_monotonicity_properties() {
        // Test general drag curve properties

        // G1 should peak somewhere in transonic region
        let mach_values: Vec<f64> = (8..20).map(|i| i as f64 * 0.1).collect(); // 0.8 to 1.9
        let g1_values: Vec<f64> = mach_values
            .iter()
            .map(|&m| get_drag_coefficient(m, &DragModel::G1))
            .collect();

        // Find maximum
        let max_value = g1_values.iter().copied().fold(0.0_f64, f64::max);
        let max_index = g1_values
            .iter()
            .position(|&x| x == max_value)
            .expect("Should find maximum in non-empty vector");
        let peak_mach = mach_values
            .get(max_index)
            .copied()
            .expect("Index should be valid");

        // Peak should be in reasonable range
        assert!(
            peak_mach > 1.0 && peak_mach < 1.6,
            "G1 peak at Mach {peak_mach}"
        );
        assert!(
            max_value > 0.5 && max_value < 1.0,
            "G1 peak value: {max_value}"
        );
    }

    #[test]
    fn test_physical_constraints() {
        let test_machs = [0.1, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0];

        for &mach in &test_machs {
            let g1_cd = get_drag_coefficient(mach, &DragModel::G1);
            let g7_cd = get_drag_coefficient(mach, &DragModel::G7);

            // All drag coefficients should be positive
            assert!(g1_cd > 0.0, "G1 CD negative at Mach {mach}: {g1_cd}");
            assert!(g7_cd > 0.0, "G7 CD negative at Mach {mach}: {g7_cd}");

            // Should be in reasonable physical ranges
            assert!(g1_cd < 2.0, "G1 CD too high at Mach {mach}: {g1_cd}");
            assert!(g7_cd < 1.5, "G7 CD too high at Mach {mach}: {g7_cd}");
        }
    }

    #[test]
    fn test_performance_characteristics() {
        // This test ensures the implementation is efficient
        use std::time::Instant;

        let start = Instant::now();

        // Perform many calculations
        for i in 0..1000 {
            let mach = 0.5 + (i as f64) * 0.004; // 0.5 to 4.5
            let _g1 = get_drag_coefficient(mach, &DragModel::G1);
            let _g7 = get_drag_coefficient(mach, &DragModel::G7);
        }

        let elapsed = start.elapsed();

        // Should complete 2000 calculations quickly (within 100ms)
        assert!(
            elapsed.as_millis() < 100,
            "Performance test took {}ms",
            elapsed.as_millis()
        );
    }

    #[test]
    fn try_new_accepts_valid_table() {
        let t = DragTable::try_new(vec![0.5, 1.0, 2.0], vec![0.20, 0.40, 0.30]).unwrap();
        assert_eq!(t.mach_values.len(), 3);
    }

    #[test]
    fn try_new_rejects_mismatched_lengths() {
        let e = DragTable::try_new(vec![0.5, 1.0, 2.0], vec![0.20, 0.40]).unwrap_err();
        assert!(e.contains("Mach") && e.contains("Cd"), "got: {e}");
    }

    #[test]
    fn try_new_rejects_too_few_points() {
        assert!(DragTable::try_new(vec![1.0], vec![0.3]).is_err());
    }

    #[test]
    fn try_new_rejects_non_ascending_mach() {
        assert!(DragTable::try_new(vec![1.0, 1.0, 2.0], vec![0.3, 0.3, 0.3]).is_err());
        assert!(DragTable::try_new(vec![2.0, 1.0], vec![0.3, 0.3]).is_err());
    }

    #[test]
    fn try_new_rejects_nonpositive_or_nonfinite_cd() {
        assert!(DragTable::try_new(vec![1.0, 2.0], vec![0.3, 0.0]).is_err());
        assert!(DragTable::try_new(vec![1.0, 2.0], vec![0.3, f64::NAN]).is_err());
    }

    #[test]
    fn interpolate_does_not_panic_on_mismatched_table() {
        // `new` is infallible; a caller could build a bad table. interpolate must not panic.
        let bad = DragTable::new(vec![0.5, 1.0, 2.0], vec![0.2]);
        let _ = bad.interpolate(0.1);
        let _ = bad.interpolate(5.0);
        let _ = bad.interpolate(1.0);
    }

    #[test]
    fn from_csv_str_parses_with_header_and_comments() {
        let csv = "# my deck\nmach,cd\n0.5, 0.230\n1.0,0.400\n2.0 , 0.300\n";
        let t = DragTable::from_csv_str(csv).unwrap();
        assert_eq!(t.mach_values, vec![0.5, 1.0, 2.0]);
        assert_eq!(t.cd_values, vec![0.230, 0.400, 0.300]);
    }

    #[test]
    fn from_csv_str_rejects_malformed_data_row() {
        // header skip is allowed once; a bad DATA row must error with a line number.
        let e = DragTable::from_csv_str("0.5,0.23\n1.0,notanumber\n").unwrap_err();
        assert!(e.contains("line 2"), "got: {e}");
    }

    #[test]
    fn from_csv_str_rejects_empty() {
        assert!(DragTable::from_csv_str("# only comments\n\n").is_err());
    }

    #[test]
    fn from_csv_str_rejects_malformed_first_data_row() {
        // first column is a valid number => it's data, not a header => must error, not vanish
        assert!(DragTable::from_csv_str("0.5\n1.0,0.4\n2.0,0.3\n").is_err());
        assert!(DragTable::from_csv_str("0.5,O.2\n1.0,0.4\n2.0,0.3\n").is_err());
    }

    #[test]
    fn from_csv_str_still_skips_textual_header() {
        // genuine header (first column non-numeric) is still tolerated
        let t = DragTable::from_csv_str("mach,cd\n0.5,0.2\n1.0,0.4\n").unwrap();
        assert_eq!(t.mach_values, vec![0.5, 1.0]);
    }

    #[test]
    fn from_csv_str_roundtrips_shipped_g7() {
        // The embedded G7 deck must load and validate through the public loader.
        let g7 = include_str!("../data/g7.csv");
        let t = DragTable::from_csv_str(g7).unwrap();
        assert!(t.mach_values.len() > 20);
    }

    // -- MBA-1386: real G2/G5/GI/GS reference tables + RA4 -----------------------------------

    #[test]
    fn test_g2_drag_coefficient_spot_values() {
        // Exact values copied from data/g2.csv (JBM mcg2.txt).
        for (mach, expected) in [(0.0, 0.2303), (0.70, 0.1702), (1.20, 0.4021), (5.0, 0.1648)] {
            let cd = get_drag_coefficient(mach, &DragModel::G2);
            assert!(
                (cd - expected).abs() < 1e-6,
                "G2 CD at Mach {mach}: expected {expected}, got {cd}"
            );
        }
    }

    #[test]
    fn test_g5_drag_coefficient_spot_values() {
        // Exact values copied from data/g5.csv (JBM mcg5.txt).
        for (mach, expected) in [(0.0, 0.1710), (0.75, 0.1463), (1.0, 0.3379), (5.0, 0.2280)] {
            let cd = get_drag_coefficient(mach, &DragModel::G5);
            assert!(
                (cd - expected).abs() < 1e-6,
                "G5 CD at Mach {mach}: expected {expected}, got {cd}"
            );
        }
    }

    #[test]
    fn test_gi_drag_coefficient_spot_values() {
        // Exact values copied from data/gi.csv (JBM mcgi.txt).
        for (mach, expected) in [(0.0, 0.2282), (0.90, 0.3170), (1.20, 0.6279), (5.0, 0.4082)] {
            let cd = get_drag_coefficient(mach, &DragModel::GI);
            assert!(
                (cd - expected).abs() < 1e-6,
                "GI CD at Mach {mach}: expected {expected}, got {cd}"
            );
        }
    }

    #[test]
    fn test_gs_drag_coefficient_spot_values() {
        // Exact values copied from data/gs.csv (JBM mcgs.txt). Table tops out at Mach 4.0.
        for (mach, expected) in [(0.0, 0.4662), (0.60, 0.5260), (1.60, 1.0090), (4.0, 0.9280)] {
            let cd = get_drag_coefficient(mach, &DragModel::GS);
            assert!(
                (cd - expected).abs() < 1e-6,
                "GS CD at Mach {mach}: expected {expected}, got {cd}"
            );
        }
    }

    #[test]
    fn test_ra4_drag_coefficient_spot_values() {
        // Exact values copied from data/ra4.csv (JBM ra4.txt). Table tops out at Mach 4.0.
        for (mach, expected) in [(0.0, 0.2283), (0.70, 0.2288), (1.15, 0.5943), (4.0, 0.4969)] {
            let cd = get_drag_coefficient(mach, &DragModel::RA4);
            assert!(
                (cd - expected).abs() < 1e-6,
                "RA4 CD at Mach {mach}: expected {expected}, got {cd}"
            );
        }
    }

    #[test]
    fn embedded_family_tables_have_ascending_mach_and_positive_cd() {
        // Structural check, via the public DragTable::from_csv_str loader (which enforces
        // strictly-ascending Mach and positive Cd through try_new): every embedded reference
        // table parses cleanly and has a sane point count.
        let decks: [(&str, &str); 7] = [
            ("G1", include_str!("../data/g1.csv")),
            ("G2", include_str!("../data/g2.csv")),
            ("G5", include_str!("../data/g5.csv")),
            ("G7", include_str!("../data/g7.csv")),
            ("GI", include_str!("../data/gi.csv")),
            ("GS", include_str!("../data/gs.csv")),
            ("RA4", include_str!("../data/ra4.csv")),
        ];
        for (name, csv) in decks {
            let table =
                DragTable::from_csv_str(csv).unwrap_or_else(|e| panic!("{name} failed to parse: {e}"));
            assert!(
                table.mach_values.len() > 20,
                "{name}: expected a high-resolution table, got {} points",
                table.mach_values.len()
            );
            assert_eq!(table.mach_values[0], 0.0, "{name}: expected Mach axis to start at 0.0");
        }

        // Also sweep the full public get_drag_coefficient API across all 9 families (including
        // the runtime-loaded G6/G8) and confirm Cd stays finite and positive everywhere.
        let models = [
            DragModel::G1,
            DragModel::G2,
            DragModel::G5,
            DragModel::G6,
            DragModel::G7,
            DragModel::G8,
            DragModel::GI,
            DragModel::GS,
            DragModel::RA4,
        ];
        for model in models {
            for i in 0..=50 {
                let mach = i as f64 * 0.1; // 0.0 ..= 5.0
                let cd = get_drag_coefficient(mach, &model);
                assert!(
                    cd.is_finite() && cd > 0.0 && cd < 2.0,
                    "{model:?} CD out of range at Mach {mach}: {cd}"
                );
            }
        }
    }

    /// Interpolate the trajectory height (position.y) at a given downrange distance,
    /// holding the nearest endpoint if the trajectory never reaches it. Mirrors the
    /// interpolation TrajectorySolver uses internally for zero-angle trials.
    fn interpolated_drop_at(points: &[crate::TrajectoryPoint], target_x_m: f64) -> f64 {
        for (i, point) in points.iter().enumerate() {
            if point.position.x >= target_x_m {
                if i == 0 {
                    return point.position.y;
                }
                let previous = &points[i - 1];
                let span = point.position.x - previous.position.x;
                if span.abs() < crate::constants::MIN_DIVISION_THRESHOLD {
                    return point.position.y;
                }
                let fraction = (target_x_m - previous.position.x) / span;
                return previous.position.y + fraction * (point.position.y - previous.position.y);
            }
        }
        points.last().map(|p| p.position.y).unwrap_or(f64::NAN)
    }

    fn drop_at_500yd_m(model: DragModel) -> f64 {
        let inputs = crate::BallisticInputs {
            bc_type: model,
            bc_value: 0.5,
            ground_threshold: -1000.0,
            ..Default::default()
        };
        let solver = crate::TrajectorySolver::new(
            inputs,
            crate::WindConditions::default(),
            crate::AtmosphericConditions::default(),
        );
        let result = solver.solve().expect("solve should succeed");
        assert!(result.max_range.is_finite());
        assert!(result.impact_velocity.is_finite() && result.impact_velocity > 0.0);
        assert!(!result.points.is_empty(), "{model:?}: solver produced no points");
        let drop = interpolated_drop_at(&result.points, 457.2); // 500 yards
        assert!(drop.is_finite(), "{model:?}: non-finite drop at 500yd");
        drop
    }

    #[test]
    fn mba1386_g2_g5_gi_gs_drop_differs_from_g1_now_that_fallback_is_gone() {
        let g1_drop = drop_at_500yd_m(DragModel::G1);
        for model in [DragModel::G2, DragModel::G5, DragModel::GI, DragModel::GS] {
            let drop = drop_at_500yd_m(model);
            assert!(
                (drop - g1_drop).abs() > 0.01,
                "{model:?} 500yd drop ({drop}) must differ from the G1 result ({g1_drop}) by \
                 more than solver noise now that the real table is wired in"
            );
        }
    }

    #[test]
    fn mba1386_ra4_solver_smoke() {
        // RA4 never fell back to G1 (it's a new variant), so just prove it solves sanely.
        let drop = drop_at_500yd_m(DragModel::RA4);
        assert!(drop < 0.0, "500yd drop should be a fall below the muzzle line: {drop}");
    }
}

/// Interpolate BC value for given Mach number from segments
pub fn interpolated_bc(mach: f64, segments: &[(f64, f64)]) -> f64 {
    if segments.is_empty() {
        return crate::constants::BC_FALLBACK_CONSERVATIVE; // Conservative fallback based on database analysis
    }

    // Get just the mach values
    let mach_values: Vec<f64> = segments.iter().map(|(m, _)| *m).collect();

    // Double-check we have values after collection
    if mach_values.is_empty() || segments.is_empty() {
        return crate::constants::BC_FALLBACK_CONSERVATIVE; // Conservative fallback based on database analysis
    }

    // Handle edge cases with safe indexing
    if let Some(first_mach) = mach_values.first() {
        if mach <= *first_mach {
            return segments.first().map(|(_, bc)| *bc).unwrap_or(0.5);
        }
    }

    if let Some(last_mach) = mach_values.last() {
        if mach >= *last_mach {
            return segments.last().map(|(_, bc)| *bc).unwrap_or(0.5);
        }
    }

    // Binary search to find the right segment with safe comparison
    let idx = match mach_values
        .binary_search_by(|&m| m.partial_cmp(&mach).unwrap_or(std::cmp::Ordering::Equal))
    {
        Ok(idx) => {
            // Exact match - safely get the BC value
            return segments.get(idx).map(|(_, bc)| *bc).unwrap_or(0.5);
        }
        Err(idx) => idx, // Insert position
    };

    // Ensure idx is valid for interpolation
    if idx == 0 || idx >= segments.len() {
        // Shouldn't happen given the edge case checks above, but be defensive
        // Use safe indexing
        let safe_idx = idx.saturating_sub(1).min(segments.len().saturating_sub(1));
        return segments.get(safe_idx).map(|(_, bc)| *bc).unwrap_or(0.5);
    }

    // Linear interpolation between the two closest points with safe indexing
    match (segments.get(idx - 1), segments.get(idx)) {
        (Some((lo_mach, lo_bc)), Some((hi_mach, hi_bc))) => {
            // Ensure denominator is not zero for safe interpolation
            let denominator = hi_mach - lo_mach;
            if denominator.abs() < crate::constants::MIN_DIVISION_THRESHOLD {
                return *lo_bc; // Return lower BC if Mach values are too close
            }
            let frac = (mach - lo_mach) / denominator;
            lo_bc + frac * (hi_bc - lo_bc)
        }
        _ => 0.5, // Fallback if indices are somehow invalid
    }
}

// Removed Python-specific function

/// MBA-1424: the built-in reference tables exposed as data.
#[cfg(test)]
mod reference_drag_table_tests {
    use super::*;

    const ALL_MODELS: [DragModel; 9] = [
        DragModel::G1,
        DragModel::G2,
        DragModel::G5,
        DragModel::G6,
        DragModel::G7,
        DragModel::G8,
        DragModel::GI,
        DragModel::GS,
        DragModel::RA4,
    ];

    /// The reason this accessor exists is so a consumer's chart cannot drift from the engine.
    /// That only holds if the data handed out is the data interpolated, so assert it directly at
    /// every tabulated point rather than trusting that both read the same static.
    #[test]
    fn the_exposed_table_is_the_one_the_solver_interpolates() {
        for model in ALL_MODELS {
            let table = reference_drag_table(&model);
            for (&mach, &cd) in table.mach_values.iter().zip(table.cd_values.iter()) {
                let interpolated = get_drag_coefficient(mach, &model);
                assert!(
                    (interpolated - cd).abs() < 1e-12,
                    "{model:?} disagrees at Mach {mach}: table {cd} vs solver {interpolated}"
                );
            }
        }
    }

    #[test]
    fn every_model_has_a_usable_table() {
        for model in ALL_MODELS {
            let table = reference_drag_table(&model);
            assert_eq!(
                table.mach_values.len(),
                table.cd_values.len(),
                "{model:?} axes differ in length"
            );
            assert!(
                table.mach_values.len() > 2,
                "{model:?} fell back to a stub table"
            );
            assert!(
                table.mach_values.windows(2).all(|w| w[1] > w[0]),
                "{model:?} Mach axis is not strictly ascending"
            );
            assert!(
                table.cd_values.iter().all(|cd| cd.is_finite() && *cd > 0.0),
                "{model:?} has a non-physical Cd"
            );
        }
    }

    /// The tables genuinely do not share a domain, which is why the accessor's contract is "read
    /// mach_values" rather than "assume 0 to 5". If this ever stops being true the docs and the
    /// CLI help that both say so should change with it.
    #[test]
    fn the_mach_domain_is_per_table_not_universal() {
        let g7_max = *reference_drag_table(&DragModel::G7).mach_values.last().unwrap();
        let gs_max = *reference_drag_table(&DragModel::GS).mach_values.last().unwrap();
        let ra4_max = *reference_drag_table(&DragModel::RA4).mach_values.last().unwrap();

        assert!(g7_max > gs_max, "G7 should extend past GS ({g7_max} vs {gs_max})");
        assert!(g7_max > ra4_max, "G7 should extend past RA4 ({g7_max} vs {ra4_max})");
    }
}
