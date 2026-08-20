// BC5D - 5-Dimensional BC Correction Table with Caliber-Specific Files
//
// This module provides offline BC corrections by loading precomputed tables
// of correction factors derived from ML model predictions. The tables are
// caliber-specific and indexed by:
//   - Weight (grains) - caliber-specific ranges
//   - Base BC (0.05-1.2)
//   - Muzzle Velocity (2000-4000 fps)
//   - Current Velocity (500-4000 fps, dense in transonic)
//   - Drag Model (G1, G7)
//
// Binary file format (BC5D v2):
//   Header (80 bytes):
//     - Magic: 4 bytes ('BC5D')
//     - Version: 4 bytes (uint32)
//     - Caliber: 4 bytes (float32)
//     - Flags: 4 bytes (uint32)
//     - Padding: 4 bytes
//     - dim_weight: 4 bytes (uint32)
//     - dim_bc: 4 bytes (uint32)
//     - dim_muzzle_vel: 4 bytes (uint32)
//     - dim_current_vel: 4 bytes (uint32)
//     - dim_drag_types: 4 bytes (uint32)
//     - timestamp: 8 bytes (uint64)
//     - checksum: 4 bytes (uint32, CRC32 of data section)
//     - api_version: 16 bytes (null-padded string)
//     - reserved: 12 bytes
//   Bin definitions:
//     - Weight bins: dim_weight * 4 bytes (float32)
//     - BC bins: dim_bc * 4 bytes (float32)
//     - Muzzle velocity bins: dim_muzzle_vel * 4 bytes (float32)
//     - Current velocity bins: dim_current_vel * 4 bytes (float32)
//   Data section:
//     - Correction factors: total_cells * 4 bytes (float32)
//     - Layout: [drag_type][weight][bc][muzzle_vel][current_vel]
//
// Correction factors are ratios: predicted_bc / base_bc
// Range: 0.5 to 1.5 (clipped during generation)

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

/// Magic bytes for BC5D format
const MAGIC: &[u8; 4] = b"BC5D";

/// Supported format version
const SUPPORTED_VERSION: u32 = 2;

/// BC5D table with 4D interpolation (drag type is discrete)
#[derive(Debug)]
pub struct Bc5dTable {
    /// Caliber this table is for
    caliber: f32,
    /// Correction data: [drag_type][weight][bc][muzzle_vel][current_vel]
    data: Vec<f32>,
    /// Weight bin values (grains)
    weight_bins: Vec<f32>,
    /// BC bin values
    bc_bins: Vec<f32>,
    /// Muzzle velocity bin values (fps)
    muzzle_vel_bins: Vec<f32>,
    /// Current velocity bin values (fps)
    current_vel_bins: Vec<f32>,
    /// Number of drag types (typically 2: G1=0, G7=1)
    num_drag_types: usize,
    /// Table version
    version: u32,
    /// API version used to generate the table
    api_version: String,
    /// Generation timestamp
    timestamp: u64,
}

/// A velocity-keyed BC schedule and the scalar BC used for any interior coverage gap.
#[cfg(any(test, target_arch = "wasm32"))]
pub(crate) struct Bc5dSegmentSchedule {
    pub(crate) segments: Vec<crate::BCSegmentData>,
    pub(crate) fallback_bc: f64,
}

/// Manager for loading caliber-specific BC5D tables
#[derive(Debug, Default)]
pub struct Bc5dTableManager {
    /// Directory containing BC5D table files
    table_dir: Option<PathBuf>,
    /// Loaded tables by caliber (rounded to 3 decimal places)
    tables: HashMap<i32, Bc5dTable>,
}

/// Error type for BC5D table operations
#[derive(Debug)]
pub enum Bc5dError {
    IoError(std::io::Error),
    InvalidMagic,
    UnsupportedVersion(u32),
    ChecksumMismatch { expected: u32, actual: u32 },
    InvalidDimensions,
    TableNotFound(f64),
    NoTableDirectory,
    /// The table's own header caliber is not the caliber of the shot it was handed to.
    /// See [`Bc5dTable::ensure_caliber_matches`] for why this is refused rather than
    /// applied or silently ignored. Both calibers are in inches.
    CaliberMismatch { table_caliber: f64, shot_caliber: f64 },
}

impl std::fmt::Display for Bc5dError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Bc5dError::IoError(e) => write!(f, "IO error: {}", e),
            Bc5dError::InvalidMagic => write!(f, "Invalid file magic (expected 'BC5D')"),
            Bc5dError::UnsupportedVersion(v) => write!(f, "Unsupported table version: {}", v),
            Bc5dError::ChecksumMismatch { expected, actual } => {
                write!(f, "Checksum mismatch: expected {:08x}, got {:08x}", expected, actual)
            }
            Bc5dError::InvalidDimensions => write!(f, "Invalid table dimensions"),
            Bc5dError::TableNotFound(cal) => write!(f, "No BC5D table found for caliber {:.3}", cal),
            Bc5dError::NoTableDirectory => write!(f, "No BC table directory configured"),
            Bc5dError::CaliberMismatch {
                table_caliber,
                shot_caliber,
            } => write!(
                f,
                "BC5D table caliber does not match the shot: table is for {:.3}, shot is {:.3} \
                 (BC5D tables are keyed to the nearest 0.001 in: {} vs {})",
                table_caliber,
                shot_caliber,
                caliber_to_key(*table_caliber),
                caliber_to_key(*shot_caliber),
            ),
        }
    }
}

impl std::error::Error for Bc5dError {}

impl From<std::io::Error> for Bc5dError {
    fn from(e: std::io::Error) -> Self {
        Bc5dError::IoError(e)
    }
}

impl Bc5dTable {
    /// Load a BC5D table from a binary file
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, Bc5dError> {
        let file = File::open(&path)?;
        Self::from_reader(BufReader::new(file))
    }

    /// Parse a BC5D table directly from an in-memory byte slice.
    ///
    /// Behaves identically to [`Bc5dTable::load`] but performs no filesystem
    /// access, which makes it usable from WASM (`wasm32-unknown-unknown`, where
    /// there is no `std::fs`). The host JS/Node layer fetches the `.bin` and
    /// hands the raw bytes across the boundary.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Bc5dError> {
        Self::from_reader(std::io::Cursor::new(bytes))
    }

    fn from_reader<R: Read>(mut reader: R) -> Result<Self, Bc5dError> {
        // Read and validate magic
        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)?;
        if &magic != MAGIC {
            return Err(Bc5dError::InvalidMagic);
        }

        // Read header fields
        let version = read_u32(&mut reader)?;
        if version != SUPPORTED_VERSION {
            return Err(Bc5dError::UnsupportedVersion(version));
        }

        let caliber = read_f32(&mut reader)?;
        let _flags = read_u32(&mut reader)?;
        let _padding = read_u32(&mut reader)?;

        let dim_weight = read_u32(&mut reader)? as usize;
        let dim_bc = read_u32(&mut reader)? as usize;
        let dim_muzzle_vel = read_u32(&mut reader)? as usize;
        let dim_current_vel = read_u32(&mut reader)? as usize;
        let dim_drag_types = read_u32(&mut reader)? as usize;

        let timestamp = read_u64(&mut reader)?;
        let stored_checksum = read_u32(&mut reader)?;

        // Read API version (16 bytes, null-terminated)
        let mut api_version_bytes = [0u8; 16];
        reader.read_exact(&mut api_version_bytes)?;
        let api_version = String::from_utf8_lossy(&api_version_bytes)
            .trim_end_matches('\0')
            .to_string();

        // Skip reserved bytes
        let mut reserved = [0u8; 12];
        reader.read_exact(&mut reserved)?;

        // Validate dimensions
        if dim_weight == 0 || dim_bc == 0 || dim_muzzle_vel == 0 || dim_current_vel == 0 || dim_drag_types == 0 {
            return Err(Bc5dError::InvalidDimensions);
        }

        // Read bin definitions
        let weight_bins = read_f32_array(&mut reader, dim_weight)?;
        let bc_bins = read_f32_array(&mut reader, dim_bc)?;
        let muzzle_vel_bins = read_f32_array(&mut reader, dim_muzzle_vel)?;
        let current_vel_bins = read_f32_array(&mut reader, dim_current_vel)?;

        // Read data section. Bound the product with checked arithmetic so a corrupt or
        // hostile file cannot overflow (debug panic / release wrap) or trigger a huge OOM
        // allocation before the trailing CRC check can reject it.
        const MAX_TOTAL_CELLS: usize = 64_000_000; // ~256 MB of f32; far above any real table
        let total_cells = dim_drag_types
            .checked_mul(dim_weight)
            .and_then(|x| x.checked_mul(dim_bc))
            .and_then(|x| x.checked_mul(dim_muzzle_vel))
            .and_then(|x| x.checked_mul(dim_current_vel))
            .filter(|&n| n <= MAX_TOTAL_CELLS)
            .ok_or(Bc5dError::InvalidDimensions)?;
        let data = read_f32_array(&mut reader, total_cells)?;

        // Verify checksum (CRC32 of bins + data)
        let mut checksum_data = Vec::new();
        for &v in &weight_bins {
            checksum_data.extend_from_slice(&v.to_le_bytes());
        }
        for &v in &bc_bins {
            checksum_data.extend_from_slice(&v.to_le_bytes());
        }
        for &v in &muzzle_vel_bins {
            checksum_data.extend_from_slice(&v.to_le_bytes());
        }
        for &v in &current_vel_bins {
            checksum_data.extend_from_slice(&v.to_le_bytes());
        }
        for &v in &data {
            checksum_data.extend_from_slice(&v.to_le_bytes());
        }

        let calculated_checksum = crc32_ieee(&checksum_data);
        if calculated_checksum != stored_checksum {
            return Err(Bc5dError::ChecksumMismatch {
                expected: stored_checksum,
                actual: calculated_checksum,
            });
        }

        Ok(Bc5dTable {
            caliber,
            data,
            weight_bins,
            bc_bins,
            muzzle_vel_bins,
            current_vel_bins,
            num_drag_types: dim_drag_types,
            version,
            api_version,
            timestamp,
        })
    }

    /// Look up a BC correction factor with 4D linear interpolation
    /// (drag type is discrete, not interpolated)
    ///
    /// # Arguments
    /// * `weight_grains` - Bullet weight in grains
    /// * `base_bc` - Published BC value
    /// * `muzzle_velocity` - Initial muzzle velocity in fps
    /// * `current_velocity` - Current bullet velocity in fps
    /// * `drag_type` - "G1" or "G7"
    ///
    /// # Returns
    /// Correction factor (multiply published BC by this value)
    pub fn lookup(
        &self,
        weight_grains: f64,
        base_bc: f64,
        muzzle_velocity: f64,
        current_velocity: f64,
        drag_type: &str,
    ) -> f64 {
        // Get drag type index (0 = G1, 1 = G7)
        let drag_idx = if drag_type.eq_ignore_ascii_case("G7") { 1 } else { 0 };

        // Clamp drag_idx to valid range
        let drag_idx = drag_idx.min(self.num_drag_types - 1);

        // Find interpolation indices and weights for each continuous dimension
        let (weight_idx, weight_w) = self.interp_idx(weight_grains as f32, &self.weight_bins);
        let (bc_idx, bc_w) = self.interp_idx(base_bc as f32, &self.bc_bins);
        let (muzzle_idx, muzzle_w) = self.interp_idx(muzzle_velocity as f32, &self.muzzle_vel_bins);
        let (current_idx, current_w) = self.interp_idx(current_velocity as f32, &self.current_vel_bins);

        // 4D linear interpolation (16 corners of a hypercube)
        let mut result = 0.0f64;

        for dw in 0..2 {
            for db in 0..2 {
                for dm in 0..2 {
                    for dc in 0..2 {
                        // Calculate weight for this corner
                        let weight = (if dw == 0 { 1.0 - weight_w } else { weight_w })
                            * (if db == 0 { 1.0 - bc_w } else { bc_w })
                            * (if dm == 0 { 1.0 - muzzle_w } else { muzzle_w })
                            * (if dc == 0 { 1.0 - current_w } else { current_w });

                        // Get clamped indices
                        let wi = (weight_idx + dw).min(self.weight_bins.len() - 1);
                        let bi = (bc_idx + db).min(self.bc_bins.len() - 1);
                        let mi = (muzzle_idx + dm).min(self.muzzle_vel_bins.len() - 1);
                        let ci = (current_idx + dc).min(self.current_vel_bins.len() - 1);

                        // Calculate flat index
                        let idx = self.flat_index(drag_idx, wi, bi, mi, ci);
                        result += weight * self.data[idx] as f64;
                    }
                }
            }
        }

        // A correction is multiplicative, so an undefined table result must be neutral rather
        // than silently becoming the most aggressive allowed degradation.
        if !result.is_finite() {
            return 1.0;
        }

        result.clamp(0.5, 1.5)
    }

    /// Get the effective BC at a given velocity
    ///
    /// This multiplies the base BC by the correction factor from the table.
    pub fn get_effective_bc(
        &self,
        weight_grains: f64,
        base_bc: f64,
        muzzle_velocity: f64,
        current_velocity: f64,
        drag_type: &str,
    ) -> f64 {
        let correction = self.lookup(weight_grains, base_bc, muzzle_velocity, current_velocity, drag_type);
        base_bc * correction
    }

    /// Generate velocity-dependent BC segments for a bullet from this table.
    ///
    /// This mirrors the CLI's `--bc-table-dir` segment synthesis: the 4D
    /// correction surface is sampled at a fixed ladder of velocity breakpoints
    /// (from 500 fps up through the muzzle velocity) and each adjacent pair
    /// becomes a [`crate::BCSegmentData`] carrying the corrected BC over that
    /// band. The solver consumes these via `inputs.bc_segments_data` +
    /// `use_bc_segments`, giving the same velocity-dependent BC degradation
    /// offline that the online solver produces.
    ///
    /// All velocities are in fps and `weight_grains` in grains, matching the
    /// table's native units. Returns `None` when the table carries no
    /// meaningful correction for this bullet (every sampled cell ~= 1.0), so
    /// callers can leave the constant published BC in place.
    pub fn generate_segments(
        &self,
        base_bc: f64,
        drag_type: &str,
        weight_grains: f64,
        muzzle_velocity_fps: Option<f64>,
    ) -> Option<Vec<crate::BCSegmentData>> {
        let mut breakpoints: Vec<f64> = vec![
            4000.0, 3500.0, 3000.0, 2700.0, 2500.0, 2300.0, 2100.0, 2000.0, 1900.0, 1800.0, 1700.0,
            1600.0, 1500.0, 1400.0, 1350.0, 1300.0, 1250.0, 1200.0, 1150.0, 1100.0, 1050.0, 1000.0,
            950.0, 900.0, 850.0, 800.0, 700.0, 600.0, 500.0,
        ];
        if let Some(mv) = muzzle_velocity_fps {
            breakpoints.push(mv);
        }

        let mut velocities: Vec<f64> = breakpoints
            .into_iter()
            .filter(|&v| v >= 500.0 && muzzle_velocity_fps.is_none_or(|mv| v <= mv))
            .collect();
        velocities.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
        velocities.dedup();

        // Correction factors were generated relative to the highest (muzzle)
        // velocity, so anchor every lookup to that reference.
        let reference_mv = velocities.first().copied().unwrap_or(3000.0);

        let mut segments: Vec<crate::BCSegmentData> = Vec::new();
        let mut any_correction = false;
        for i in 0..velocities.len().saturating_sub(1) {
            let vel_max = velocities[i];
            let vel_min = velocities[i + 1];
            let vel_mid = (vel_max + vel_min) / 2.0;

            let correction =
                self.lookup(weight_grains, base_bc, reference_mv, vel_mid, drag_type);
            if (correction - 1.0).abs() > 1e-6 {
                any_correction = true;
            }
            segments.push(crate::BCSegmentData {
                velocity_min: vel_min,
                velocity_max: vel_max,
                bc_value: base_bc * correction,
            });
        }

        if any_correction && !segments.is_empty() {
            Some(segments)
        } else {
            None
        }
    }

    /// Generate the BC5D schedule consumed by the WASM frontend.
    #[cfg(any(test, target_arch = "wasm32"))]
    pub(crate) fn generate_segment_schedule(
        &self,
        base_bc: f64,
        drag_type: &str,
        weight_grains: f64,
        muzzle_velocity_fps: f64,
    ) -> Option<Bc5dSegmentSchedule> {
        let segments =
            self.generate_segments(base_bc, drag_type, weight_grains, Some(muzzle_velocity_fps))?;
        let fallback_bc = self.get_effective_bc(
            weight_grains,
            base_bc,
            muzzle_velocity_fps,
            muzzle_velocity_fps,
            drag_type,
        );

        Some(Bc5dSegmentSchedule {
            segments,
            fallback_bc,
        })
    }

    /// Find interpolation index and weight for a value in bins
    fn interp_idx(&self, value: f32, bins: &[f32]) -> (usize, f64) {
        if bins.len() < 2 || value.is_nan() {
            return (0, 0.0);
        }

        // Handle out of range (clamp to edges)
        if value <= bins[0] {
            return (0, 0.0);
        }
        if value >= bins[bins.len() - 1] {
            return (bins.len().saturating_sub(2), 1.0);
        }

        // Binary search for interval containing value
        let last_interval = bins.len().saturating_sub(2);
        let idx = match bins.binary_search_by(|probe| {
            probe
                .partial_cmp(&value)
                .unwrap_or(std::cmp::Ordering::Equal)
        }) {
            Ok(i) => i.saturating_sub(1).min(last_interval),
            Err(i) => i.saturating_sub(1).min(last_interval),
        };

        // Calculate interpolation weight
        let low = bins[idx];
        let high = bins[idx + 1];
        let weight = if high > low {
            ((value - low) / (high - low)) as f64
        } else {
            0.0
        };

        (idx, weight)
    }

    /// Calculate flat array index from 5D indices
    fn flat_index(&self, drag_idx: usize, weight_idx: usize, bc_idx: usize, muzzle_idx: usize, current_idx: usize) -> usize {
        let n_weight = self.weight_bins.len();
        let n_bc = self.bc_bins.len();
        let n_muzzle = self.muzzle_vel_bins.len();
        let n_current = self.current_vel_bins.len();

        drag_idx * (n_weight * n_bc * n_muzzle * n_current)
            + weight_idx * (n_bc * n_muzzle * n_current)
            + bc_idx * (n_muzzle * n_current)
            + muzzle_idx * n_current
            + current_idx
    }

    /// Get caliber this table is for
    pub fn caliber(&self) -> f32 {
        self.caliber
    }

    /// The table's own caliber as a BC5D key ([`caliber_to_key`]) — exactly the value
    /// [`Self::ensure_caliber_matches`] compares against, and the value `bc5d.info`
    /// reports as `caliber_key` so an app can pre-check a downloaded table itself.
    pub fn caliber_key(&self) -> i32 {
        caliber_to_key(self.caliber as f64)
    }

    /// Refuse this table for a shot of a different caliber.
    ///
    /// `shot_caliber_in` is the shot's bullet diameter in INCHES.
    ///
    /// # Why this is an error and not a warning or a silent skip
    ///
    /// Nothing in the lookup path fails on a foreign table: [`Self::lookup`] clamps
    /// out-of-range values to the edge bins, so a table for another caliber still
    /// returns a plausible-looking correction for every cell and
    /// [`Self::generate_segments`] still emits a full ladder. Measured: the published
    /// `bc5d_224.bin` handed a 175 gr / G1 0.505 / 2600 fps .308 shot yields a
    /// 25-segment ladder whose segment BCs are 0.4710..0.5114 where the .308 table
    /// gives 0.4989..0.5072 — ~6.7 % off in the low bands, with no diagnostic
    /// anywhere; a .308 table on a .243 shot measures -17.9 % effective BC. A
    /// wrong-caliber table is therefore WORSE than no table at all: no table leaves
    /// the published BC intact, while a wrong one silently biases every row. So the
    /// caller is refused outright — never corrected with foreign data, and never
    /// quietly downgraded to an uncorrected solve either, because a caller that asked
    /// for a table and got an unannotated uncorrected answer cannot tell.
    ///
    /// # Matching rule
    ///
    /// Equality of the 3-digit BC5D caliber key ([`caliber_to_key`]) — the same key
    /// `find_table_file` uses to choose `bc5d_<key>.bin` for `--bc-table-dir`. The
    /// guard therefore accepts exactly the diameters for which the CLI would have
    /// selected this table's own filename: precedent decides the rule, so the CLI and
    /// the path/bytes consumers (bridge cards, bridge solve, solve-json, WASM) cannot
    /// disagree about what "this table is for my shot" means.
    ///
    /// In practice that is a half-thousandth tolerance — a `308` table accepts
    /// `[0.3075, 0.3085)` — expressed as a bucket rather than a centered epsilon.
    /// Bucketing matters: the header caliber is an `f32`, so `0.308` arrives as
    /// 0.30799998, and rounding to the nearest thousandth absorbs that representation
    /// error, whereas a centered `|a - b| <= 0.0005` comparison would have to carry a
    /// separate fudge for it.
    pub fn ensure_caliber_matches(&self, shot_caliber_in: f64) -> Result<(), Bc5dError> {
        if self.caliber_key() == caliber_to_key(shot_caliber_in) {
            return Ok(());
        }
        Err(Bc5dError::CaliberMismatch {
            table_caliber: self.caliber as f64,
            shot_caliber: shot_caliber_in,
        })
    }

    /// Get table version
    pub fn version(&self) -> u32 {
        self.version
    }

    /// Get API version used to generate the table
    pub fn api_version(&self) -> &str {
        &self.api_version
    }

    /// Get generation timestamp
    pub fn timestamp(&self) -> u64 {
        self.timestamp
    }

    /// Get total number of cells in the table
    pub fn total_cells(&self) -> usize {
        self.data.len()
    }

    /// Bin counts per axis: `(weight, bc, muzzle_vel, current_vel, drag_types)`.
    ///
    /// The structured counterpart of [`Self::dimensions_str`], for callers that
    /// report table metadata over a wire (e.g. the bridge's `bc5d.info`).
    pub fn bin_counts(&self) -> (usize, usize, usize, usize, usize) {
        (
            self.weight_bins.len(),
            self.bc_bins.len(),
            self.muzzle_vel_bins.len(),
            self.current_vel_bins.len(),
            self.num_drag_types,
        )
    }

    /// Get table dimensions as a string
    pub fn dimensions_str(&self) -> String {
        format!(
            "{}x{}x{}x{}x{} (weight x bc x muzzle_vel x current_vel x drag_types)",
            self.weight_bins.len(),
            self.bc_bins.len(),
            self.muzzle_vel_bins.len(),
            self.current_vel_bins.len(),
            self.num_drag_types
        )
    }

    /// Get weight range
    pub fn weight_range(&self) -> (f32, f32) {
        (*self.weight_bins.first().unwrap_or(&0.0), *self.weight_bins.last().unwrap_or(&0.0))
    }

    /// Get velocity range
    pub fn velocity_range(&self) -> (f32, f32) {
        (*self.current_vel_bins.first().unwrap_or(&0.0), *self.current_vel_bins.last().unwrap_or(&0.0))
    }
}

impl Bc5dTableManager {
    /// Create a new table manager with a directory path
    pub fn new<P: AsRef<Path>>(table_dir: P) -> Self {
        Bc5dTableManager {
            table_dir: Some(table_dir.as_ref().to_path_buf()),
            tables: HashMap::new(),
        }
    }

    /// Create an empty manager (no table directory)
    pub fn empty() -> Self {
        Bc5dTableManager {
            table_dir: None,
            tables: HashMap::new(),
        }
    }

    /// Get or load the table for a caliber
    ///
    /// Tables are cached after first load.
    ///
    /// The file is selected by NAME (`bc5d_<key>.bin`) but verified by CONTENT: a file
    /// whose header caliber is not `caliber` is rejected with
    /// [`Bc5dError::CaliberMismatch`] and never cached, so a rotated manifest, a
    /// hand-copied `.bin`, or a future generator bug cannot silently bias every lookup
    /// (see [`Bc5dTable::ensure_caliber_matches`]).
    pub fn get_table(&mut self, caliber: f64) -> Result<&Bc5dTable, Bc5dError> {
        let caliber_key = caliber_to_key(caliber);

        // Check if already loaded
        if self.tables.contains_key(&caliber_key) {
            return Ok(self.tables.get(&caliber_key).unwrap());
        }

        // Need to load
        let table_dir = self.table_dir.as_ref().ok_or(Bc5dError::NoTableDirectory)?;
        let table_path = find_table_file(table_dir, caliber)?;
        let table = Bc5dTable::load(&table_path)?;
        table.ensure_caliber_matches(caliber)?;
        self.tables.insert(caliber_key, table);
        Ok(self.tables.get(&caliber_key).unwrap())
    }

    /// Look up BC correction for a bullet
    pub fn lookup(
        &mut self,
        caliber: f64,
        weight_grains: f64,
        base_bc: f64,
        muzzle_velocity: f64,
        current_velocity: f64,
        drag_type: &str,
    ) -> Result<f64, Bc5dError> {
        let table = self.get_table(caliber)?;
        Ok(table.lookup(weight_grains, base_bc, muzzle_velocity, current_velocity, drag_type))
    }

    /// Get effective BC with correction applied
    pub fn get_effective_bc(
        &mut self,
        caliber: f64,
        weight_grains: f64,
        base_bc: f64,
        muzzle_velocity: f64,
        current_velocity: f64,
        drag_type: &str,
    ) -> Result<f64, Bc5dError> {
        let table = self.get_table(caliber)?;
        Ok(table.get_effective_bc(weight_grains, base_bc, muzzle_velocity, current_velocity, drag_type))
    }

    /// Check if a table is available for a caliber
    pub fn has_table(&self, caliber: f64) -> bool {
        if let Some(ref table_dir) = self.table_dir {
            find_table_file(table_dir, caliber).is_ok()
        } else {
            false
        }
    }

    /// List available calibers in the table directory
    pub fn available_calibers(&self) -> Vec<f64> {
        let mut calibers = Vec::new();
        if let Some(ref table_dir) = self.table_dir {
            if let Ok(entries) = std::fs::read_dir(table_dir) {
                for entry in entries.flatten() {
                    let path = entry.path();
                    if let Some(ext) = path.extension() {
                        if ext == "bin" {
                            if let Some(stem) = path.file_stem() {
                                let name = stem.to_string_lossy();
                                if let Some(caliber) = name.strip_prefix("bc5d_") {
                                    // Parse caliber from filename (e.g., bc5d_308.bin -> 0.308)
                                    if let Ok(cal_int) = caliber.parse::<i32>() {
                                        calibers.push(cal_int as f64 / 1000.0);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        calibers.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        calibers
    }
}

/// The canonical BC5D caliber key: a caliber in inches rounded to the nearest
/// thousandth, as an integer (0.308 -> 308, 0.224 -> 224).
///
/// This is the ONE key BC5D identity is expressed in. `find_table_file` builds
/// `bc5d_<key>.bin` from it (so it decides which file `--bc-table-dir` picks),
/// [`Bc5dTableManager`] caches by it, and [`Bc5dTable::ensure_caliber_matches`]
/// compares by it — one rule, so the CLI's file selection and every path/bytes
/// consumer's identity check cannot drift apart.
pub fn caliber_to_key(caliber: f64) -> i32 {
    // A non-finite caliber saturates to 0 here (Rust's float->int cast), which is not
    // any real table's key, so it falls out as a mismatch rather than matching anything.
    (caliber * 1000.0).round() as i32
}

/// Find the table file for a caliber
fn find_table_file(table_dir: &Path, caliber: f64) -> Result<PathBuf, Bc5dError> {
    let caliber_int = (caliber * 1000.0).round() as i32;
    let filename = format!("bc5d_{}.bin", caliber_int);
    let path = table_dir.join(&filename);

    if path.exists() {
        return Ok(path);
    }

    // Try common variations
    let variations = [
        format!("bc5d_{:03}.bin", caliber_int),
        format!("bc5d_0{}.bin", caliber_int),
    ];

    for var in &variations {
        let var_path = table_dir.join(var);
        if var_path.exists() {
            return Ok(var_path);
        }
    }

    Err(Bc5dError::TableNotFound(caliber))
}

// Helper functions for reading binary data

fn read_u32<R: Read>(reader: &mut R) -> Result<u32, std::io::Error> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_u64<R: Read>(reader: &mut R) -> Result<u64, std::io::Error> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn read_f32<R: Read>(reader: &mut R) -> Result<f32, std::io::Error> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(f32::from_le_bytes(buf))
}

// A newer clippy stable added `chunks_exact_to_as_chunks`, which flags a constant chunk
// size and suggests `as_chunks`. Suppressed rather than restructured: this is binary-format
// parsing, `as_chunks` changes the element type from `&[u8]` to `&[u8; N]` and so ripples
// into every use inside the loop, and the change would land in the middle of a 13-platform
// release. `unknown_lints` is allowed alongside it so toolchains predating the lint do not
// warn on the name. Adopting `as_chunks` properly is a follow-up.
#[allow(unknown_lints, clippy::chunks_exact_to_as_chunks)]
fn read_f32_array<R: Read>(reader: &mut R, count: usize) -> Result<Vec<f32>, std::io::Error> {
    // Defensive bounds: reject absurd lengths from corrupt/hostile files before
    // allocating, and guard the byte-count multiply against overflow.
    const MAX_ELEMS: usize = 64_000_000; // 256 MB of f32
    if count > MAX_ELEMS {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "f32 array length too large",
        ));
    }
    let byte_len = count.checked_mul(4).ok_or_else(|| {
        std::io::Error::new(std::io::ErrorKind::InvalidData, "f32 array length overflow")
    })?;
    let mut data = vec![0f32; count];
    let mut buf = vec![0u8; byte_len];
    reader.read_exact(&mut buf)?;

    for (i, chunk) in buf.chunks_exact(4).enumerate() {
        data[i] = f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
    }

    Ok(data)
}

/// Simple CRC32 (IEEE polynomial) implementation
pub(crate) fn crc32_ieee(data: &[u8]) -> u32 {
    const TABLE: [u32; 256] = make_crc32_table();
    let mut crc = 0xFFFFFFFFu32;
    for &byte in data {
        let idx = ((crc ^ byte as u32) & 0xFF) as usize;
        crc = (crc >> 8) ^ TABLE[idx];
    }
    !crc
}

const fn make_crc32_table() -> [u32; 256] {
    const POLY: u32 = 0xEDB88320;
    let mut table = [0u32; 256];
    let mut i = 0;
    while i < 256 {
        let mut crc = i as u32;
        let mut j = 0;
        while j < 8 {
            if crc & 1 != 0 {
                crc = (crc >> 1) ^ POLY;
            } else {
                crc >>= 1;
            }
            j += 1;
        }
        table[i] = crc;
        i += 1;
    }
    table
}

/// Process-wide cache of BC5D tables loaded from explicit filesystem paths.
///
/// The bridge/solve-json surfaces accept a caller-supplied table PATH (mobile apps
/// download the `.bin` themselves and hand the engine a file path), and PARSING a
/// several-MB table on every card or solve call would dominate the request. What the
/// cache saves is the parse; the read and CRC are the price of knowing what is
/// actually on disk.
///
/// Entries are keyed by `(canonical path, file size, CRC32 of the file's bytes)` —
/// i.e. by CONTENT. An earlier version keyed on `(canonical path, file size, mtime)`
/// and could serve a stale parsed table when a file was replaced in place by
/// same-size content within one filesystem mtime tick: the key was unchanged, so the
/// new bytes were never read. That is a live scenario here, not a theoretical one —
/// a table-set refresh overwrites `bc5d_<caliber>.bin` in place, and a regenerated
/// table with identical dimensions has identical size. It reached a release because
/// mtime granularity is fine enough on macOS and Linux to hide it, and only surfaced
/// on an OpenBSD guest whose granularity is coarse enough to collide.
///
/// Size is retained alongside the CRC purely as a second, free discriminator.
/// The cache is bounded (oldest entry evicted at capacity).
///
/// Filesystem-only by construction, so the whole module is compiled out on
/// `wasm32` (where WASM callers pass table BYTES via `loadBc5dTable` instead).
#[cfg(not(target_arch = "wasm32"))]
pub mod path_cache {
    use super::{Bc5dError, Bc5dTable};
    use std::path::{Path, PathBuf};
    use std::sync::{Arc, OnceLock, RwLock};

    /// Small on purpose: a mobile app realistically has one or two calibers live.
    const CACHE_CAPACITY: usize = 4;

    #[derive(Debug, Clone, PartialEq, Eq)]
    struct CacheKey {
        canonical_path: PathBuf,
        file_size: u64,
        /// CRC32 of the file's ENTIRE byte content. Not the checksum field stored
        /// inside the table (that describes only the data section, and a corrupted
        /// data byte leaves it untouched — exactly the case that must invalidate).
        content_crc: u32,
    }

    type CacheEntries = Vec<(CacheKey, Arc<Bc5dTable>)>;

    fn cache() -> &'static RwLock<CacheEntries> {
        static CACHE: OnceLock<RwLock<CacheEntries>> = OnceLock::new();
        CACHE.get_or_init(|| RwLock::new(Vec::new()))
    }

    /// Load a BC5D table from `path`, verifying the header (magic, version,
    /// dimensions) and the stored CRC32 exactly as [`Bc5dTable::load`] does, with
    /// the parsed result cached process-wide.
    ///
    /// A cache hit requires the canonical path, file size, AND a CRC32 over the
    /// file's bytes to match, so ANY change to the content invalidates — including a
    /// same-size in-place replacement and single-byte corruption, neither of which a
    /// timestamp reliably distinguishes. Corrupt, truncated, or missing files are
    /// never cached.
    ///
    /// The file is therefore read on every call; only the parse is cached. At the
    /// once-per-request call sites (bridge cards, bridge `solve`, `solve-json`,
    /// `bc5d.info`) that trade is not measurable against the solve itself.
    pub fn load_verified(path: &Path) -> Result<Arc<Bc5dTable>, Bc5dError> {
        let canonical_path = std::fs::canonicalize(path)?;
        // Read BEFORE consulting the cache: the bytes are the identity, so there is
        // nothing trustworthy to look up until we have them.
        let bytes = std::fs::read(&canonical_path)?;
        let key = CacheKey {
            canonical_path,
            file_size: bytes.len() as u64,
            content_crc: super::crc32_ieee(&bytes),
        };

        if let Ok(entries) = cache().read() {
            if let Some((_, table)) = entries.iter().find(|(k, _)| *k == key) {
                return Ok(Arc::clone(table));
            }
        }

        let table = Arc::new(Bc5dTable::from_bytes(&bytes)?);

        if let Ok(mut entries) = cache().write() {
            // Re-check under the write lock: another thread may have inserted the
            // same key between our read miss and here.
            if let Some((_, existing)) = entries.iter().find(|(k, _)| *k == key) {
                return Ok(Arc::clone(existing));
            }
            if entries.len() >= CACHE_CAPACITY {
                entries.remove(0);
            }
            entries.push((key, Arc::clone(&table)));
        }
        Ok(table)
    }

    /// [`load_verified`] plus the caliber-identity guard: the loaded table must be for
    /// the caliber of the shot that is about to use it, per
    /// [`Bc5dTable::ensure_caliber_matches`] (`shot_caliber_in` in INCHES).
    ///
    /// This is the entry point every consumer that HAS a shot must use — bridge cards,
    /// bridge `solve`, and `solve-json` all take a caller-supplied table path, and a
    /// path says nothing about content. Only surfaces with no shot in hand (e.g. the
    /// bridge's `bc5d.info`, which just describes a file) call [`load_verified`]
    /// directly. Keeping the guard inside the loader is deliberate: adding a fourth
    /// path-based consumer cannot forget it.
    pub fn load_verified_for_caliber(
        path: &Path,
        shot_caliber_in: f64,
    ) -> Result<Arc<Bc5dTable>, Bc5dError> {
        let table = load_verified(path)?;
        table.ensure_caliber_matches(shot_caliber_in)?;
        Ok(table)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_table() -> Bc5dTable {
        // Create a small test table with known values
        let weight_bins = vec![100.0, 150.0, 200.0];
        let bc_bins = vec![0.3, 0.4, 0.5];
        let muzzle_vel_bins = vec![2500.0, 3000.0];
        let current_vel_bins = vec![1000.0, 2000.0, 3000.0];
        let num_drag_types = 2;

        // Total cells: 2 * 3 * 3 * 2 * 3 = 108
        let total = num_drag_types * weight_bins.len() * bc_bins.len() * muzzle_vel_bins.len() * current_vel_bins.len();
        let mut data = vec![1.0f32; total];

        // Set some non-uniform values for testing interpolation
        // At weight=150, bc=0.4, muzzle=2750 (interpolated), current=2000, G1
        // We'll set corners to test 4D interpolation
        data[0] = 0.95; // First corner
        data[total - 1] = 1.05; // Last corner

        Bc5dTable {
            caliber: 0.308,
            data,
            weight_bins,
            bc_bins,
            muzzle_vel_bins,
            current_vel_bins,
            num_drag_types,
            version: 2,
            api_version: "test".to_string(),
            timestamp: 0,
        }
    }

    /// The defect that reached 0.33.3: the cache keyed on `(path, size, mtime)`, so a
    /// file replaced IN PLACE by same-size content within one mtime tick kept its key
    /// and the stale parsed table was served. The sibling test above only caught it on
    /// filesystems whose timestamp granularity happens to collide (it failed on an
    /// OpenBSD guest and passed everywhere else), so this one removes the luck: it
    /// pins both writes to an IDENTICAL mtime and asserts content still decides.
    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn path_cache_detects_same_size_replacement_under_an_identical_mtime() {
        use super::path_cache;
        use std::time::{Duration, SystemTime};

        let dir = std::env::temp_dir().join(format!(
            "bc5d_same_mtime_{}_{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("bc5d_308.bin");

        // A fixed, explicitly-set mtime that BOTH writes will carry, emulating a
        // filesystem that cannot separate them. `File::set_modified` is std, so this
        // stays portable instead of shelling out to `touch`.
        let pinned = SystemTime::UNIX_EPOCH + Duration::from_secs(1_600_000_000);
        let pin = |p: &std::path::Path| {
            let f = std::fs::OpenOptions::new().write(true).open(p).unwrap();
            f.set_modified(pinned).unwrap();
            f.sync_all().unwrap();
        };

        let original = create_test_table();
        let good = serialize_test_table(&original);
        std::fs::write(&path, &good).unwrap();
        pin(&path);
        let first = path_cache::load_verified(&path).expect("valid table loads");

        // 1. A same-size VALID replacement must be parsed, not served from cache.
        let mut replacement = create_test_table();
        let last = replacement.data.len() - 1;
        replacement.data[last] = 0.5; // different content, identical dimensions => identical size
        let replaced = serialize_test_table(&replacement);
        assert_eq!(replaced.len(), good.len(), "test requires an identical size");
        std::fs::write(&path, &replaced).unwrap();
        pin(&path);
        assert_eq!(
            std::fs::metadata(&path).unwrap().modified().unwrap(),
            pinned,
            "both writes must share one mtime for this test to mean anything"
        );
        let second = path_cache::load_verified(&path).expect("replacement loads");
        assert!(
            !std::sync::Arc::ptr_eq(&first, &second),
            "a same-size replacement under an identical mtime must NOT be served from cache"
        );

        // 2. Same-size CORRUPTION must be reported, never served from cache. The stored
        //    checksum field is untouched here, so only hashing the real bytes catches it.
        let mut corrupt = good.clone();
        *corrupt.last_mut().unwrap() ^= 0xFF;
        assert_eq!(corrupt.len(), good.len());
        std::fs::write(&path, &corrupt).unwrap();
        pin(&path);
        assert!(
            matches!(
                path_cache::load_verified(&path),
                Err(Bc5dError::ChecksumMismatch { .. })
            ),
            "corruption under an identical mtime must be detected, not cached"
        );

        // 3. Restoring the original bytes is a cache HIT again — content, not history.
        std::fs::write(&path, &good).unwrap();
        pin(&path);
        let restored = path_cache::load_verified(&path).expect("original loads again");
        assert!(
            std::sync::Arc::ptr_eq(&first, &restored),
            "identical content must still be served from the cache"
        );

        let _ = std::fs::remove_dir_all(&dir);
    }

    fn create_single_cell_test_table() -> Bc5dTable {
        Bc5dTable {
            caliber: 0.308,
            data: vec![0.875],
            weight_bins: vec![168.0],
            bc_bins: vec![0.4],
            muzzle_vel_bins: vec![2500.0],
            current_vel_bins: vec![2000.0],
            num_drag_types: 1,
            version: 2,
            api_version: "test".to_string(),
            timestamp: 0,
        }
    }

    /// Serialize a table into the BC5D v2 `.bin` byte layout so we can exercise
    /// the `from_bytes` parser without depending on an external file.
    fn serialize_test_table(t: &Bc5dTable) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(MAGIC);
        out.extend_from_slice(&t.version.to_le_bytes());
        out.extend_from_slice(&t.caliber.to_le_bytes());
        out.extend_from_slice(&0u32.to_le_bytes()); // flags
        out.extend_from_slice(&0u32.to_le_bytes()); // padding
        out.extend_from_slice(&(t.weight_bins.len() as u32).to_le_bytes());
        out.extend_from_slice(&(t.bc_bins.len() as u32).to_le_bytes());
        out.extend_from_slice(&(t.muzzle_vel_bins.len() as u32).to_le_bytes());
        out.extend_from_slice(&(t.current_vel_bins.len() as u32).to_le_bytes());
        out.extend_from_slice(&(t.num_drag_types as u32).to_le_bytes());
        out.extend_from_slice(&t.timestamp.to_le_bytes());

        // Checksum is CRC32 of bins + data, in declaration order.
        let mut checksum_data = Vec::new();
        for v in t.weight_bins.iter().chain(&t.bc_bins).chain(&t.muzzle_vel_bins)
            .chain(&t.current_vel_bins).chain(&t.data) {
            checksum_data.extend_from_slice(&v.to_le_bytes());
        }
        out.extend_from_slice(&crc32_ieee(&checksum_data).to_le_bytes());

        let mut api = [0u8; 16];
        let bytes = t.api_version.as_bytes();
        api[..bytes.len().min(16)].copy_from_slice(&bytes[..bytes.len().min(16)]);
        out.extend_from_slice(&api);
        out.extend_from_slice(&[0u8; 12]); // reserved

        for v in t.weight_bins.iter().chain(&t.bc_bins).chain(&t.muzzle_vel_bins)
            .chain(&t.current_vel_bins).chain(&t.data) {
            out.extend_from_slice(&v.to_le_bytes());
        }
        out
    }

    #[test]
    fn test_from_bytes_roundtrip() {
        let original = create_test_table();
        let bytes = serialize_test_table(&original);
        let parsed = Bc5dTable::from_bytes(&bytes).expect("from_bytes should parse");

        assert_eq!(parsed.caliber, original.caliber);
        assert_eq!(parsed.num_drag_types, original.num_drag_types);
        assert_eq!(parsed.weight_bins, original.weight_bins);
        assert_eq!(parsed.current_vel_bins, original.current_vel_bins);
        assert_eq!(parsed.data, original.data);
        assert_eq!(parsed.api_version, original.api_version);

        // A corrupted body must be rejected by the CRC check.
        let mut bad = bytes.clone();
        *bad.last_mut().unwrap() ^= 0xFF;
        assert!(Bc5dTable::from_bytes(&bad).is_err());
    }

    #[test]
    fn test_generate_segments() {
        // A table whose corrections are all exactly 1.0 carries no useful
        // correction, so generate_segments returns None (leave published BC).
        let mut uniform = create_test_table();
        uniform.data.iter_mut().for_each(|v| *v = 1.0);
        assert!(uniform
            .generate_segments(0.4, "G1", 150.0, Some(2700.0))
            .is_none());

        // A table with a real (0.9) correction across the sampled slice must
        // produce contiguous, descending velocity segments carrying bc*corr.
        let mut corrected = create_test_table();
        corrected.data.iter_mut().for_each(|v| *v = 0.9);
        let segments = corrected
            .generate_segments(0.4, "G1", 150.0, Some(2700.0))
            .expect("segments expected for a table with corrections");
        assert!(!segments.is_empty());
        for w in segments.windows(2) {
            // Bands are contiguous and descend in velocity.
            assert!((segments[0].velocity_max - w[0].velocity_max).abs() >= 0.0);
            assert!(w[0].velocity_min >= w[1].velocity_max - 1e-6);
        }
        for s in &segments {
            assert!((s.bc_value - 0.4 * 0.9).abs() < 1e-6); // base_bc * correction
            assert!(s.velocity_max > s.velocity_min);
        }
    }

    #[test]
    fn segment_schedule_carries_muzzle_corrected_fallback_bc() {
        let table = create_single_cell_test_table();
        let base_bc = 0.4;
        let schedule = table
            .generate_segment_schedule(base_bc, "G1", 168.0, 2500.0)
            .expect("uniform non-neutral correction should produce a schedule");
        let expected_fallback = table.get_effective_bc(168.0, base_bc, 2500.0, 2500.0, "G1");

        assert!(!schedule.segments.is_empty());
        assert_eq!(expected_fallback.to_bits(), (base_bc * 0.875).to_bits());
        assert_eq!(schedule.fallback_bc.to_bits(), expected_fallback.to_bits());
    }

    #[test]
    fn test_interp_idx_in_range() {
        let table = create_test_table();

        // Test middle of range
        let (idx, weight) = table.interp_idx(125.0, &table.weight_bins);
        assert_eq!(idx, 0);
        assert!((weight - 0.5).abs() < 0.01);

        // Test at bin boundary
        let (idx, weight) = table.interp_idx(150.0, &table.weight_bins);
        assert_eq!(idx, 0);
        assert!((weight - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_interp_idx_out_of_range() {
        let table = create_test_table();

        // Test below range
        let (idx, weight) = table.interp_idx(50.0, &table.weight_bins);
        assert_eq!(idx, 0);
        assert_eq!(weight, 0.0);

        // Test above range
        let (idx, weight) = table.interp_idx(250.0, &table.weight_bins);
        assert_eq!(idx, 1); // len - 2
        assert_eq!(weight, 1.0);
    }

    #[test]
    fn test_interp_idx_nan_defaults_to_first_bin() {
        let table = create_test_table();

        assert_eq!(table.interp_idx(f32::NAN, &table.weight_bins), (0, 0.0));
        assert_eq!(
            table.interp_idx(f32::NEG_INFINITY, &table.weight_bins),
            (0, 0.0)
        );
        assert_eq!(
            table.interp_idx(f32::INFINITY, &table.weight_bins),
            (1, 1.0)
        );
    }

    #[test]
    fn test_lookup_nan_with_single_bin_axes_uses_only_cell() {
        let table = create_single_cell_test_table();

        assert_eq!(table.lookup(f64::NAN, 0.4, 2500.0, 2000.0, "G1"), 0.875);
    }

    #[test]
    fn test_lookup_non_finite_table_cells_are_neutral() {
        for cell in [f32::NAN, f32::INFINITY, f32::NEG_INFINITY] {
            let mut source = create_test_table();
            source.data.fill(cell);
            let bytes = serialize_test_table(&source);
            let table = Bc5dTable::from_bytes(&bytes).expect("CRC-valid table should load");

            assert_eq!(
                table.lookup(125.0, 0.35, 2750.0, 1500.0, "G1"),
                1.0,
                "non-finite cell {cell:?} must produce a neutral correction"
            );
            assert_eq!(
                table.get_effective_bc(125.0, 0.35, 2750.0, 1500.0, "G1"),
                0.35,
                "neutral correction must preserve base BC for {cell:?}"
            );
        }
    }

    #[test]
    fn test_lookup_returns_valid_range() {
        let table = create_test_table();

        let correction = table.lookup(150.0, 0.4, 2750.0, 2000.0, "G1");
        assert!((0.5..=1.5).contains(&correction));

        let correction = table.lookup(150.0, 0.4, 2750.0, 2000.0, "G7");
        assert!((0.5..=1.5).contains(&correction));
    }

    #[test]
    fn test_effective_bc() {
        let table = create_test_table();

        let base_bc = 0.4;
        let effective = table.get_effective_bc(150.0, base_bc, 2750.0, 2000.0, "G1");

        // Effective BC should be base_bc * correction
        assert!(effective >= base_bc * 0.5 && effective <= base_bc * 1.5);
    }

    #[test]
    fn test_caliber_to_key() {
        assert_eq!(caliber_to_key(0.308), 308);
        assert_eq!(caliber_to_key(0.224), 224);
        assert_eq!(caliber_to_key(0.338), 338);
    }

    /// The matching rule's boundaries, pinned. A `308` table is the bucket of diameters
    /// that round to 0.308 at the thousandth — precisely the diameters for which
    /// `find_table_file` would have chosen this table's own `bc5d_308.bin`.
    #[test]
    fn ensure_caliber_matches_accepts_the_rounding_bucket_and_refuses_outside_it() {
        let table = create_test_table(); // header caliber 0.308 (as f32)

        // The f32 header (0.30799998) must not cost us the exact match.
        assert_eq!(table.caliber_key(), 308);
        assert!(table.ensure_caliber_matches(0.308).is_ok());

        // Inclusive at the bottom edge: 0.3075 * 1000 is exactly 307.5, which rounds
        // half-away-from-zero to 308.
        assert!(table.ensure_caliber_matches(0.3075).is_ok());
        assert!(table.ensure_caliber_matches(0.3084).is_ok());

        // Exclusive at the top edge: 0.3085 * 1000 is exactly 308.5 -> key 309.
        assert!(matches!(
            table.ensure_caliber_matches(0.3085),
            Err(Bc5dError::CaliberMismatch { .. })
        ));
        assert!(matches!(
            table.ensure_caliber_matches(0.3074),
            Err(Bc5dError::CaliberMismatch { .. })
        ));

        // The real-world failure mode: a whole different caliber, named in the message.
        let err = table
            .ensure_caliber_matches(0.224)
            .expect_err("a .224 shot must not be served a .308 table");
        let message = err.to_string();
        assert!(
            message.contains("table is for 0.308, shot is 0.224"),
            "the error must name both calibers: {message}"
        );

        // A garbage diameter cannot accidentally match a real table either.
        assert!(table.ensure_caliber_matches(f64::NAN).is_err());
        assert!(table.ensure_caliber_matches(0.0).is_err());
    }

    /// A file named `bc5d_308.bin` whose CONTENT is a .224 table must be refused by the
    /// CLI's manager, and must not be cached (so a retry cannot resurrect it).
    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn manager_refuses_a_file_whose_header_caliber_is_foreign() {
        let dir = std::env::temp_dir().join(format!(
            "bc5d-mislabeled-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();

        let mut foreign = create_test_table();
        foreign.caliber = 0.224;
        foreign.data.fill(0.9);
        std::fs::write(dir.join("bc5d_308.bin"), serialize_test_table(&foreign)).unwrap();

        let mut manager = Bc5dTableManager::new(&dir);
        let err = manager
            .get_table(0.308)
            .expect_err("a mislabeled table must be refused, not applied");
        assert!(matches!(err, Bc5dError::CaliberMismatch { .. }), "{err}");
        // Not cached: the second attempt fails the same way rather than succeeding.
        assert!(matches!(
            manager.get_table(0.308),
            Err(Bc5dError::CaliberMismatch { .. })
        ));
        // And no correction leaks out through the convenience wrappers.
        assert!(manager
            .lookup(0.308, 168.0, 0.4, 2500.0, 2000.0, "G1")
            .is_err());

        // The same bytes ARE usable for the caliber they actually describe.
        std::fs::write(dir.join("bc5d_224.bin"), serialize_test_table(&foreign)).unwrap();
        assert!(manager.get_table(0.224).is_ok());

        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// `load_verified_for_caliber` is the guarded loader every shot-bearing consumer
    /// uses: same bytes, accepted for their own caliber and refused for another.
    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn path_cache_guarded_loader_refuses_a_foreign_caliber() {
        let dir = std::env::temp_dir().join(format!(
            "bc5d-guarded-load-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("bc5d_308.bin");
        std::fs::write(&path, serialize_test_table(&create_test_table())).unwrap();

        assert!(path_cache::load_verified_for_caliber(&path, 0.308).is_ok());
        let err = path_cache::load_verified_for_caliber(&path, 0.243)
            .expect_err("a .243 shot must be refused a .308 table");
        assert!(matches!(err, Bc5dError::CaliberMismatch { .. }), "{err}");
        assert!(
            err.to_string().contains("table is for 0.308, shot is 0.243"),
            "{err}"
        );

        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn test_table_metadata() {
        let table = create_test_table();
        assert!((table.caliber() - 0.308).abs() < 0.001);
        assert_eq!(table.version(), 2);
        assert_eq!(table.api_version(), "test");
    }

    #[test]
    fn test_crc32() {
        // Test with known CRC32 value
        let data = b"123456789";
        let crc = crc32_ieee(data);
        assert_eq!(crc, 0xCBF43926);
    }

    #[test]
    fn test_bin_counts_matches_dimensions() {
        let table = create_test_table();
        assert_eq!(table.bin_counts(), (3, 3, 2, 3, 2));
    }

    /// The path cache must hand back the SAME parsed table for repeated loads of an
    /// unchanged file, reject corruption instead of caching it, and pick up an
    /// in-place replacement (size change breaks the key).
    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn path_cache_reuses_parsed_tables_and_detects_replacement() {
        let dir = std::env::temp_dir().join(format!(
            "bc5d-path-cache-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("bc5d_308.bin");

        let table = create_test_table();
        let bytes = serialize_test_table(&table);
        std::fs::write(&path, &bytes).unwrap();

        let first = path_cache::load_verified(&path).expect("valid table loads");
        let second = path_cache::load_verified(&path).expect("cached table loads");
        assert!(
            std::sync::Arc::ptr_eq(&first, &second),
            "an unchanged file must be served from the cache"
        );

        // Replace the file with a differently sized (still valid) table: the key
        // changes, so the next load parses the new content.
        let mut replacement = create_test_table();
        replacement.weight_bins.push(250.0);
        let extra_cells = replacement.num_drag_types
            * replacement.bc_bins.len()
            * replacement.muzzle_vel_bins.len()
            * replacement.current_vel_bins.len();
        replacement
            .data
            .resize(replacement.data.len() + extra_cells, 0.9f32);
        std::fs::write(&path, serialize_test_table(&replacement)).unwrap();
        let third = path_cache::load_verified(&path).expect("replacement loads");
        assert!(!std::sync::Arc::ptr_eq(&first, &third));
        assert_eq!(third.bin_counts().0, 4, "replacement content must be parsed");

        // Corruption is a clean error, not a cached table.
        let mut corrupt = serialize_test_table(&table);
        *corrupt.last_mut().unwrap() ^= 0xFF;
        std::fs::write(&path, corrupt).unwrap();
        assert!(matches!(
            path_cache::load_verified(&path),
            Err(Bc5dError::ChecksumMismatch { .. })
        ));

        std::fs::remove_dir_all(&dir).unwrap();
    }
}
