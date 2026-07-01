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
        let mut reader = BufReader::new(file);

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

        // Clamp result to valid range
        result.max(0.5).min(1.5)
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

    /// Find interpolation index and weight for a value in bins
    fn interp_idx(&self, value: f32, bins: &[f32]) -> (usize, f64) {
        if bins.is_empty() {
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
        let idx = match bins.binary_search_by(|probe| {
            probe.partial_cmp(&value).unwrap_or(std::cmp::Ordering::Equal)
        }) {
            Ok(i) => i.saturating_sub(1).min(bins.len() - 2),
            Err(i) => i.saturating_sub(1).min(bins.len() - 2),
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
                                if name.starts_with("bc5d_") {
                                    // Parse caliber from filename (e.g., bc5d_308.bin -> 0.308)
                                    if let Ok(cal_int) = name[5..].parse::<i32>() {
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

/// Convert caliber to integer key (multiply by 1000)
fn caliber_to_key(caliber: f64) -> i32 {
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
    fn test_lookup_returns_valid_range() {
        let table = create_test_table();

        let correction = table.lookup(150.0, 0.4, 2750.0, 2000.0, "G1");
        assert!(correction >= 0.5 && correction <= 1.5);

        let correction = table.lookup(150.0, 0.4, 2750.0, 2000.0, "G7");
        assert!(correction >= 0.5 && correction <= 1.5);
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
}
