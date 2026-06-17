// BC Correction Table - 5D Lookup Table with Linear Interpolation
//
// This module provides offline BC corrections by loading a precomputed table
// of correction factors. The table is indexed by:
//   - BC value (published BC)
//   - BC type (G1 or G7)
//   - Bullet mass (grains)
//   - Bullet length (inches)
//   - Velocity (fps)
//
// Binary file format (BCCR):
//   Header (64 bytes):
//     - Magic: 4 bytes ('BCCR')
//     - Version: 4 bytes (uint32)
//     - Flags: 4 bytes (uint32)
//     - num_bc: 4 bytes (uint32)
//     - num_mass: 4 bytes (uint32)
//     - num_length: 4 bytes (uint32)
//     - num_velocity: 4 bytes (uint32)
//     - num_bc_types: 4 bytes (uint32)
//     - timestamp: 8 bytes (uint64)
//     - checksum: 4 bytes (uint32, CRC32)
//     - reserved: 16 bytes
//   Bin definitions:
//     - BC bins: num_bc * 4 bytes (float32)
//     - Mass bins: num_mass * 4 bytes (float32)
//     - Length bins: num_length * 4 bytes (float32)
//     - Velocity bins: num_velocity * 4 bytes (float32)
//   Data section:
//     - Correction factors: total_cells * 4 bytes (float32)
//     - Layout: [bc_type][bc][mass][length][velocity]

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

/// BC correction table with 5D interpolation
#[derive(Debug)]
pub struct BcCorrectionTable {
    /// Table data: [bc_type][bc][mass][length][velocity]
    data: Vec<f32>,
    /// BC bin values
    bc_bins: Vec<f32>,
    /// Mass bin values (grains)
    mass_bins: Vec<f32>,
    /// Length bin values (inches)
    length_bins: Vec<f32>,
    /// Velocity bin values (fps)
    velocity_bins: Vec<f32>,
    /// Number of BC types (typically 2: G1, G7)
    num_types: usize,
    /// Table version
    version: u32,
}

/// Error type for BC table operations
#[derive(Debug)]
pub enum BcTableError {
    IoError(std::io::Error),
    InvalidMagic,
    UnsupportedVersion(u32),
    ChecksumMismatch,
    InvalidDimensions,
}

impl std::fmt::Display for BcTableError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BcTableError::IoError(e) => write!(f, "IO error: {}", e),
            BcTableError::InvalidMagic => write!(f, "Invalid file magic (expected 'BCCR')"),
            BcTableError::UnsupportedVersion(v) => write!(f, "Unsupported table version: {}", v),
            BcTableError::ChecksumMismatch => write!(f, "Data checksum mismatch"),
            BcTableError::InvalidDimensions => write!(f, "Invalid table dimensions"),
        }
    }
}

impl std::error::Error for BcTableError {}

impl From<std::io::Error> for BcTableError {
    fn from(e: std::io::Error) -> Self {
        BcTableError::IoError(e)
    }
}

impl BcCorrectionTable {
    /// Load a BC correction table from a binary file
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, BcTableError> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        // Read header (64 bytes)
        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)?;
        if &magic != b"BCCR" {
            return Err(BcTableError::InvalidMagic);
        }

        let version = read_u32(&mut reader)?;
        if version != 1 {
            return Err(BcTableError::UnsupportedVersion(version));
        }

        let _flags = read_u32(&mut reader)?;
        let num_bc = read_u32(&mut reader)? as usize;
        let num_mass = read_u32(&mut reader)? as usize;
        let num_length = read_u32(&mut reader)? as usize;
        let num_velocity = read_u32(&mut reader)? as usize;
        let num_types = read_u32(&mut reader)? as usize;
        let _timestamp = read_u64(&mut reader)?;
        let _checksum = read_u32(&mut reader)?;

        // Skip reserved bytes
        let mut reserved = [0u8; 16];
        reader.read_exact(&mut reserved)?;

        // Validate dimensions
        if num_bc == 0 || num_mass == 0 || num_length == 0 || num_velocity == 0 || num_types == 0 {
            return Err(BcTableError::InvalidDimensions);
        }

        // Read bin definitions
        let bc_bins = read_f32_array(&mut reader, num_bc)?;
        let mass_bins = read_f32_array(&mut reader, num_mass)?;
        let length_bins = read_f32_array(&mut reader, num_length)?;
        let velocity_bins = read_f32_array(&mut reader, num_velocity)?;

        // Read data section. Bound the product with checked arithmetic so a corrupt/hostile
        // file cannot overflow (debug panic / release wrap to a too-short `data`, which then
        // panics on the unchecked index in lookup()) or trigger a huge OOM allocation.
        // Mirrors the bc_table_5d.rs hardening.
        const MAX_TOTAL_CELLS: usize = 64_000_000; // ~256 MB of f32; far above any real table
        let total_cells = num_types
            .checked_mul(num_bc)
            .and_then(|x| x.checked_mul(num_mass))
            .and_then(|x| x.checked_mul(num_length))
            .and_then(|x| x.checked_mul(num_velocity))
            .filter(|&n| n <= MAX_TOTAL_CELLS)
            .ok_or(BcTableError::InvalidDimensions)?;
        let data = read_f32_array(&mut reader, total_cells)?;

        Ok(BcCorrectionTable {
            data,
            bc_bins,
            mass_bins,
            length_bins,
            velocity_bins,
            num_types,
            version,
        })
    }

    /// Look up a BC correction factor with 5D linear interpolation
    ///
    /// # Arguments
    /// * `bc` - Published BC value
    /// * `bc_type` - "G1" or "G7"
    /// * `mass` - Bullet mass in grains
    /// * `length` - Bullet length in inches
    /// * `velocity` - Current velocity in fps
    ///
    /// # Returns
    /// Correction factor (multiply published BC by this value)
    pub fn lookup(&self, bc: f64, bc_type: &str, mass: f64, length: f64, velocity: f64) -> f64 {
        // Get type index (0 = G1, 1 = G7), clamped to the table's type axis so a
        // G1-only table (num_types == 1) queried with "G7" stays in bounds instead of
        // producing an out-of-range flat_index (mirrors Bc5dTable::lookup).
        let type_idx =
            (if bc_type.to_uppercase() == "G1" { 0 } else { 1 }).min(self.num_types.saturating_sub(1));

        // Find interpolation indices and weights for each dimension
        let (bc_idx, bc_weight) = self.interp_idx(bc as f32, &self.bc_bins);
        let (mass_idx, mass_weight) = self.interp_idx(mass as f32, &self.mass_bins);
        let (length_idx, length_weight) = self.interp_idx(length as f32, &self.length_bins);
        let (vel_idx, vel_weight) = self.interp_idx(velocity as f32, &self.velocity_bins);

        // 4D linear interpolation (type is discrete, not interpolated)
        let mut result = 0.0f64;

        for di in 0..2 {
            for dj in 0..2 {
                for dk in 0..2 {
                    for dl in 0..2 {
                        // Calculate weight for this corner
                        let w = (if di == 0 { 1.0 - bc_weight } else { bc_weight })
                            * (if dj == 0 { 1.0 - mass_weight } else { mass_weight })
                            * (if dk == 0 { 1.0 - length_weight } else { length_weight })
                            * (if dl == 0 { 1.0 - vel_weight } else { vel_weight });

                        // Get clamped indices
                        let i = (bc_idx + di).min(self.bc_bins.len() - 1);
                        let j = (mass_idx + dj).min(self.mass_bins.len() - 1);
                        let k = (length_idx + dk).min(self.length_bins.len() - 1);
                        let l = (vel_idx + dl).min(self.velocity_bins.len() - 1);

                        // Calculate flat index
                        let idx = self.flat_index(type_idx, i, j, k, l);
                        // Bounds-safe access (true defense-in-depth): the continuous axes are
                        // clamped above and type_idx is clamped to num_types, so idx is always
                        // in bounds for any valid table; the unwrap_or only guards corrupt data.
                        result += w * self.data.get(idx).copied().unwrap_or(0.0) as f64;
                    }
                }
            }
        }

        result
    }

    /// Get the BC correction factor for a given bullet and velocity
    /// This is a convenience method that estimates bullet length from mass and caliber
    ///
    /// # Arguments
    /// * `bc` - Published BC value
    /// * `bc_type` - "G1" or "G7"
    /// * `mass_grains` - Bullet mass in grains
    /// * `caliber_inches` - Bullet caliber in inches
    /// * `velocity_fps` - Current velocity in fps
    pub fn lookup_with_caliber(
        &self,
        bc: f64,
        bc_type: &str,
        mass_grains: f64,
        caliber_inches: f64,
        velocity_fps: f64,
    ) -> f64 {
        // Estimate bullet length from caliber (typical rifle bullets are 3-4 calibers long)
        let estimated_length = caliber_inches * 3.5;
        self.lookup(bc, bc_type, mass_grains, estimated_length, velocity_fps)
    }

    /// Find interpolation index and weight for a value in bins
    fn interp_idx(&self, value: f32, bins: &[f32]) -> (usize, f64) {
        if bins.is_empty() {
            return (0, 0.0);
        }

        // Handle out of range
        if value <= bins[0] {
            return (0, 0.0);
        }
        if value >= bins[bins.len() - 1] {
            return (bins.len().saturating_sub(2), 1.0);
        }

        // Binary search for interval
        let idx = match bins
            .binary_search_by(|probe| probe.partial_cmp(&value).unwrap_or(std::cmp::Ordering::Equal))
        {
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
    fn flat_index(&self, type_idx: usize, bc_idx: usize, mass_idx: usize, length_idx: usize, vel_idx: usize) -> usize {
        let n_bc = self.bc_bins.len();
        let n_mass = self.mass_bins.len();
        let n_length = self.length_bins.len();
        let n_velocity = self.velocity_bins.len();

        type_idx * (n_bc * n_mass * n_length * n_velocity)
            + bc_idx * (n_mass * n_length * n_velocity)
            + mass_idx * (n_length * n_velocity)
            + length_idx * n_velocity
            + vel_idx
    }

    /// Get table version
    pub fn version(&self) -> u32 {
        self.version
    }

    /// Get total number of cells in the table
    pub fn total_cells(&self) -> usize {
        self.data.len()
    }

    /// Get table dimensions as a string
    pub fn dimensions_str(&self) -> String {
        format!(
            "{}x{}x{}x{}x{}",
            self.bc_bins.len(),
            self.mass_bins.len(),
            self.length_bins.len(),
            self.velocity_bins.len(),
            self.num_types
        )
    }
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

fn read_f32_array<R: Read>(reader: &mut R, count: usize) -> Result<Vec<f32>, std::io::Error> {
    // Defensive bounds: reject absurd lengths from corrupt/hostile files before allocating,
    // and guard the byte-count multiply against overflow (mirrors bc_table_5d.rs).
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interp_idx_in_range() {
        let table = BcCorrectionTable {
            data: vec![1.0; 100],
            bc_bins: vec![0.1, 0.2, 0.3, 0.4, 0.5],
            mass_bins: vec![100.0, 150.0, 200.0],
            length_bins: vec![1.0, 1.2, 1.4],
            velocity_bins: vec![3000.0, 2500.0, 2000.0],
            num_types: 2,
            version: 1,
        };

        // Test middle of range
        let (idx, weight) = table.interp_idx(0.25, &table.bc_bins);
        assert_eq!(idx, 1);
        assert!((weight - 0.5).abs() < 0.01);

        // Test at bin boundary
        let (idx, weight) = table.interp_idx(0.2, &table.bc_bins);
        assert_eq!(idx, 0);
        assert!((weight - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_interp_idx_out_of_range() {
        let table = BcCorrectionTable {
            data: vec![1.0; 100],
            bc_bins: vec![0.1, 0.2, 0.3, 0.4, 0.5],
            mass_bins: vec![100.0, 150.0, 200.0],
            length_bins: vec![1.0, 1.2, 1.4],
            velocity_bins: vec![3000.0, 2500.0, 2000.0],
            num_types: 2,
            version: 1,
        };

        // Test below range
        let (idx, weight) = table.interp_idx(0.05, &table.bc_bins);
        assert_eq!(idx, 0);
        assert_eq!(weight, 0.0);

        // Test above range
        let (idx, weight) = table.interp_idx(0.6, &table.bc_bins);
        assert_eq!(idx, 3); // len - 2
        assert_eq!(weight, 1.0);
    }
}
