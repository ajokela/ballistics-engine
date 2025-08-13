/// Pure Rust tests for drag coefficient functionality that don't require PyO3
#[cfg(test)]
mod pure_rust_tests {
    use crate::drag::{DragTable, load_drag_table};
    use std::path::Path;

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
        
        // Extrapolation should not be negative
        let below = table.interpolate(0.5);
        assert!(below > 0.0);
        
        let above = table.interpolate(3.0);
        assert!(above > 0.0);
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
        
        // Should be between reasonable bounds
        assert!(result > 0.1 && result < 1.0);
        
        // Should be smooth (not exactly linear)
        let linear_result = table.linear_interpolate(1.25, 1);
        // Results should be related but not identical
        assert!((result - linear_result).abs() < 0.5);
    }

    #[test]
    fn test_load_drag_table_fallback() {
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
    fn test_extrapolation_bounds() {
        let table = DragTable::new(vec![1.0, 2.0, 3.0], vec![0.4, 0.6, 0.3]);
        
        // Test extrapolation below range - should not be negative
        let cd_low = table.interpolate(0.0);
        assert!(cd_low > 0.01, "Low extrapolation should be positive: {cd_low}");
        
        // Test extrapolation above range - should not be negative
        let cd_high = table.interpolate(10.0);
        assert!(cd_high > 0.01, "High extrapolation should be positive: {cd_high}");
    }

    #[test]
    fn test_interpolation_continuity() {
        let table = DragTable::new(
            vec![0.5, 1.0, 1.5, 2.0, 2.5, 3.0], 
            vec![0.2, 0.4, 0.6, 0.5, 0.3, 0.25]
        );
        
        // Test that interpolation is reasonably continuous
        for i in 1..=25 {
            let mach = 0.5 + i as f64 * 0.1;
            let cd_before = table.interpolate(mach - 0.01);
            let cd_after = table.interpolate(mach + 0.01);
            let difference = (cd_after - cd_before).abs();
            
            // Should not have large jumps
            assert!(difference < 0.1, "Large discontinuity at Mach {mach}: {cd_before} vs {cd_after}");
        }
    }

    #[test]
    fn test_physical_constraints() {
        let table = DragTable::new(
            vec![0.5, 1.0, 1.5, 2.0, 3.0, 5.0], 
            vec![0.2, 0.5, 0.6, 0.4, 0.3, 0.2]
        );
        
        let test_machs = [0.1, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0];
        
        for &mach in &test_machs {
            let cd = table.interpolate(mach);
            
            // All drag coefficients should be positive
            assert!(cd > 0.0, "CD negative at Mach {mach}: {cd}");
            
            // Should be in reasonable physical ranges
            assert!(cd < 5.0, "CD too high at Mach {mach}: {cd}");
        }
    }

    #[test]
    fn test_performance_basic() {
        let table = DragTable::new(
            vec![0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0], 
            vec![0.2, 0.5, 0.6, 0.4, 0.3, 0.25, 0.22, 0.2]
        );
        
        // Perform many interpolations to ensure no performance issues
        for i in 0..1000 {
            let mach = 0.5 + (i as f64) * 0.004; // 0.5 to 4.5
            let _cd = table.interpolate(mach);
        }
        
        // If we get here without timing out, performance is acceptable
        assert!(true);
    }
}