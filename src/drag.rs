/// Drag coefficient calculations for ballistics using actual drag table data
use std::path::Path;
use ndarray::ArrayD;
use once_cell::sync::Lazy;
use pyo3::prelude::*;
use crate::DragModel;
use crate::transonic_drag::{transonic_correction, get_projectile_shape, ProjectileShape};

/// Drag table data structure
#[derive(Clone)]
pub struct DragTable {
    pub mach_values: Vec<f64>,
    pub cd_values: Vec<f64>,
}

impl DragTable {
    /// Create a new drag table from mach and cd arrays
    pub fn new(mach_values: Vec<f64>, cd_values: Vec<f64>) -> Self {
        Self { mach_values, cd_values }
    }

    /// Interpolate drag coefficient for given Mach number using cubic spline-like interpolation
    pub fn interpolate(&self, mach: f64) -> f64 {
        let n = self.mach_values.len();
        
        if n == 0 {
            return 0.5; // Fallback
        }
        
        if n == 1 {
            return self.cd_values[0];
        }

        // Handle out-of-bounds cases with extrapolation
        if mach <= self.mach_values[0] {
            // Linear extrapolation below range, but prevent negative values
            if n >= 2 {
                let slope = (self.cd_values[1] - self.cd_values[0]) / 
                           (self.mach_values[1] - self.mach_values[0]);
                let extrapolated = self.cd_values[0] + slope * (mach - self.mach_values[0]);
                // Clamp to minimum reasonable value to prevent negative drag coefficients
                return extrapolated.max(0.01);
            }
            return self.cd_values[0];
        }
        
        if mach >= self.mach_values[n - 1] {
            // Linear extrapolation above range, but prevent negative values
            if n >= 2 {
                let slope = (self.cd_values[n - 1] - self.cd_values[n - 2]) / 
                           (self.mach_values[n - 1] - self.mach_values[n - 2]);
                let extrapolated = self.cd_values[n - 1] + slope * (mach - self.mach_values[n - 1]);
                // Clamp to minimum reasonable value to prevent negative drag coefficients
                return extrapolated.max(0.01);
            }
            return self.cd_values[n - 1];
        }

        // Find the segment containing the mach value
        let mut idx = 0;
        for i in 0..n - 1 {
            if mach >= self.mach_values[i] && mach <= self.mach_values[i + 1] {
                idx = i;
                break;
            }
        }

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

    /// Cubic interpolation using 4 points (similar to cubic spline)
    pub fn cubic_interpolate(&self, mach: f64, idx: usize) -> f64 {
        // Ensure we have enough points for cubic interpolation
        if idx == 0 || idx + 1 >= self.mach_values.len() || 
           idx + 1 >= self.cd_values.len() {
            // Fall back to linear interpolation if not enough points
            return self.linear_interpolate(mach, idx);
        }
        
        // Use points at idx-1, idx, idx+1, idx+2
        let x = [
            self.mach_values[idx - 1],
            self.mach_values[idx],
            self.mach_values[idx + 1],
            if idx + 2 < self.mach_values.len() { self.mach_values[idx + 2] } else { self.mach_values[idx + 1] }
        ];
        let y = [
            self.cd_values[idx - 1],
            self.cd_values[idx],
            self.cd_values[idx + 1],
            if idx + 2 < self.cd_values.len() { self.cd_values[idx + 2] } else { self.cd_values[idx + 1] }
        ];

        // Catmull-Rom spline interpolation
        // Ensure denominator is not zero
        let denominator = x[2] - x[1];
        if denominator.abs() < crate::constants::MIN_DIVISION_THRESHOLD {
            return y[1]; // Return the value at the current point if points are too close
        }
        let t = (mach - x[1]) / denominator;
        let t2 = t * t;
        let t3 = t2 * t;

        let a0 = -0.5 * y[0] + 1.5 * y[1] - 1.5 * y[2] + 0.5 * y[3];
        let a1 = y[0] - 2.5 * y[1] + 2.0 * y[2] - 0.5 * y[3];
        let a2 = -0.5 * y[0] + 0.5 * y[2];
        let a3 = y[1];

        a0 * t3 + a1 * t2 + a2 * t + a3
    }
}

/// Load drag table from NumPy binary file or CSV fallback
pub fn load_drag_table(drag_tables_dir: &Path, filename: &str, fallback_data: &[(f64, f64)]) -> DragTable {
    // Try to load NumPy binary file first
    let npy_path = drag_tables_dir.join(format!("{filename}.npy"));
    if let Ok(array) = ndarray_npy::read_npy::<_, ArrayD<f64>>(&npy_path) {
        if let Ok(array_2d) = array.into_dimensionality::<ndarray::Ix2>() {
            let mach_values: Vec<f64> = array_2d.column(0).to_vec();
            let cd_values: Vec<f64> = array_2d.column(1).to_vec();
            return DragTable::new(mach_values, cd_values);
        }
    }

    // Fallback to CSV file
    let csv_path = drag_tables_dir.join(format!("{filename}.csv"));
    if let Ok(mut reader) = csv::Reader::from_path(&csv_path) {
        let mut mach_values = Vec::new();
        let mut cd_values = Vec::new();
        
        for result in reader.records() {
            if let Ok(record) = result {
                if record.len() >= 2 {
                    if let (Ok(mach), Ok(cd)) = (record[0].parse::<f64>(), record[1].parse::<f64>()) {
                        mach_values.push(mach);
                        cd_values.push(cd);
                    }
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

/// G1 drag table with lazy loading
static G1_DRAG_TABLE: Lazy<DragTable> = Lazy::new(|| {
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

    if let Some(drag_dir) = find_drag_tables_dir() {
        load_drag_table(&drag_dir, "g1", &fallback_data)
    } else {
        // Use fallback data if directory not found
        let mach_values: Vec<f64> = fallback_data.iter().map(|(m, _)| *m).collect();
        let cd_values: Vec<f64> = fallback_data.iter().map(|(_, cd)| *cd).collect();
        DragTable::new(mach_values, cd_values)
    }
});

/// G7 drag table with lazy loading
static G7_DRAG_TABLE: Lazy<DragTable> = Lazy::new(|| {
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

    if let Some(drag_dir) = find_drag_tables_dir() {
        load_drag_table(&drag_dir, "g7", &fallback_data)
    } else {
        // Use fallback data if directory not found
        let mach_values: Vec<f64> = fallback_data.iter().map(|(m, _)| *m).collect();
        let cd_values: Vec<f64> = fallback_data.iter().map(|(_, cd)| *cd).collect();
        DragTable::new(mach_values, cd_values)
    }
});

/// Get drag coefficient for given Mach number and drag model
pub fn get_drag_coefficient(mach: f64, drag_model: &DragModel) -> f64 {
    match drag_model {
        DragModel::G1 => G1_DRAG_TABLE.interpolate(mach),
        DragModel::G7 => G7_DRAG_TABLE.interpolate(mach),
        _ => G1_DRAG_TABLE.interpolate(mach), // Default to G1 for other models
    }
}

/// Get drag coefficient with optional transonic corrections
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
                    get_projectile_shape(cal, weight, match drag_model {
                        DragModel::G1 => "G1",
                        DragModel::G7 => "G7",
                        _ => "G1", // Default to G1
                    })
                } else {
                    ProjectileShape::Spitzer // Default
                }
            }
        };
        
        // Apply correction
        transonic_correction(mach, base_cd, shape, true)
    } else {
        base_cd
    }
}

/// Get drag coefficient with optional transonic and Reynolds corrections
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
    
    // Apply Reynolds corrections for low velocities (subsonic only)
    if apply_reynolds_correction && mach < 1.0 {
        if let (Some(v), Some(cal), Some(rho), Some(temp)) = 
            (velocity_mps, caliber, air_density_kg_m3, temperature_c) {
            use crate::reynolds::apply_reynolds_correction;
            cd = apply_reynolds_correction(cd, v, cal, rho, temp, mach);
        }
    }
    
    cd
}

#[cfg(test)]
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
    fn test_drag_coefficient_continuity() {
        // Test that drag coefficient function is smooth
        for mach in [0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0] {
            let cd_before = get_drag_coefficient(mach - 0.01, &DragModel::G1);
            let cd_after = get_drag_coefficient(mach + 0.01, &DragModel::G1);
            let difference = (cd_after - cd_before).abs();
            assert!(difference < 0.05, "Large discontinuity at Mach {mach}: {cd_before} vs {cd_after}");
        }
    }

    #[test]
    fn test_extrapolation_bounds() {
        // Test extrapolation below range
        let cd_low = get_drag_coefficient(0.0, &DragModel::G1);
        assert!(cd_low > 0.01 && cd_low < 0.5, "Low Mach G1: {cd_low}");
        
        // Test extrapolation above range - should be clamped to positive values
        let cd_high = get_drag_coefficient(10.0, &DragModel::G1);
        assert!(cd_high > 0.01, "High Mach G1 should be positive: {cd_high}");
        
        // Same for G7
        let cd_low_g7 = get_drag_coefficient(0.0, &DragModel::G7);
        assert!(cd_low_g7 > 0.01, "Low Mach G7 should be positive: {cd_low_g7}");
        
        let cd_high_g7 = get_drag_coefficient(20.0, &DragModel::G7);
        assert!(cd_high_g7 > 0.01, "High Mach G7 should be positive: {cd_high_g7}");
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
        
        // Extrapolation
        let below = table.interpolate(0.5);
        assert!(below.abs() > 1e-10); // Should have some value
        
        let above = table.interpolate(3.0);
        assert!(above.abs() > 1e-10); // Should have some value
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
        let table = DragTable::new(
            vec![0.5, 1.0, 1.5, 2.0, 2.5],
            vec![0.2, 0.4, 0.6, 0.5, 0.3]
        );
        
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
        assert!((g1_mach1 - 0.4805).abs() < 0.01, "G1 at Mach 1.0: {g1_mach1}");
        
        // G7 at Mach 1.0 should be around 0.3803
        let g7_mach1 = get_drag_coefficient(1.0, &DragModel::G7);
        assert!((g7_mach1 - 0.3803).abs() < 0.01, "G7 at Mach 1.0: {g7_mach1}");
        
        // G1 should generally be higher than G7 in transonic region
        assert!(g1_mach1 > g7_mach1, "G1 should be > G7 at Mach 1.0");
    }

    #[test]
    fn test_monotonicity_properties() {
        // Test general drag curve properties
        
        // G1 should peak somewhere in transonic region
        let mach_values: Vec<f64> = (8..20).map(|i| i as f64 * 0.1).collect(); // 0.8 to 1.9
        let g1_values: Vec<f64> = mach_values.iter()
            .map(|&m| get_drag_coefficient(m, &DragModel::G1))
            .collect();
        
        // Find maximum
        let max_value = g1_values.iter().copied().fold(0.0_f64, f64::max);
        let max_index = g1_values.iter().position(|&x| x == max_value)
            .expect("Should find maximum in non-empty vector");
        let peak_mach = mach_values.get(max_index).copied()
            .expect("Index should be valid");
        
        // Peak should be in reasonable range
        assert!(peak_mach > 1.0 && peak_mach < 1.6, "G1 peak at Mach {peak_mach}");
        assert!(max_value > 0.5 && max_value < 1.0, "G1 peak value: {max_value}");
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
        assert!(elapsed.as_millis() < 100, "Performance test took {}ms", elapsed.as_millis());
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
    let idx = match mach_values.binary_search_by(|&m| {
        m.partial_cmp(&mach).unwrap_or(std::cmp::Ordering::Equal)
    }) {
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

/// Python-exposed function for BC interpolation
#[pyfunction]
#[pyo3(name = "interpolated_bc_rust")]
pub fn interpolated_bc_py(mach: f64, segments: Vec<(f64, f64)>) -> PyResult<f64> {
    Ok(interpolated_bc(mach, &segments))
}