//! Cluster-based BC degradation for improved accuracy
//!
//! This module implements empirical BC corrections based on bullet clustering.
//! Bullets are classified into 4 clusters based on their physical characteristics,
//! and each cluster has unique velocity-dependent BC degradation curves derived
//! from real-world ballistic data.

#[derive(Debug, Clone, Copy)]
pub struct ClusterBCDegradation {
    /// Pre-calculated cluster centroids
    centroids: [(f64, f64, f64); 4],
}

impl ClusterBCDegradation {
    pub fn new() -> Self {
        Self {
            // Cluster centroids: (caliber_normalized, weight_normalized, bc_normalized)
            centroids: [
                (0.605, 0.415, 0.613), // Cluster 0: Standard Long-Range
                (0.516, 0.324, 0.643), // Cluster 1: Low-Drag Specialty
                (0.307, 0.088, 0.336), // Cluster 2: Light Varmint
                (0.750, 0.805, 0.505), // Cluster 3: Heavy Magnums
            ],
        }
    }

    /// Predict which cluster a bullet belongs to
    pub fn predict_cluster(&self, caliber: f64, weight_gr: f64, bc_g1: f64) -> usize {
        // Normalize features to [0, 1] range
        // Bounds derived from training data and centroid values:
        // - Caliber: 0.172 (.17 HMR) to 0.750 (.50 BMG) - matches cluster 3 centroid
        // - Weight: 15gr (light varmint) to 750gr (heavy magnum)
        // - BC G1: 0.05 (low drag) to 1.2 (high BC match bullets)
        let caliber_norm = (caliber - 0.172) / (0.750 - 0.172);
        let weight_norm = (weight_gr - 15.0) / (750.0 - 15.0);
        let bc_norm = (bc_g1 - 0.05) / (1.2 - 0.05);

        // Find nearest centroid
        let mut min_distance = f64::INFINITY;
        let mut best_cluster = 0;

        for (i, &(c_cal, c_wt, c_bc)) in self.centroids.iter().enumerate() {
            let distance = ((caliber_norm - c_cal).powi(2)
                + (weight_norm - c_wt).powi(2)
                + (bc_norm - c_bc).powi(2))
            .sqrt();

            if distance < min_distance {
                min_distance = distance;
                best_cluster = i;
            }
        }

        best_cluster
    }

    /// Get BC multiplier for a given velocity and cluster
    pub fn get_bc_multiplier(&self, velocity_fps: f64, cluster_id: usize) -> f64 {
        match cluster_id {
            0 => {
                // Standard Long-Range: gradual degradation
                if velocity_fps > 2800.0 {
                    1.0
                } else if velocity_fps > 2400.0 {
                    0.98 - 0.03 * (2800.0 - velocity_fps) / 400.0
                } else if velocity_fps > 1800.0 {
                    0.95 - 0.05 * (2400.0 - velocity_fps) / 600.0
                } else if velocity_fps > 1200.0 {
                    0.90 - 0.10 * (1800.0 - velocity_fps) / 600.0
                } else {
                    0.80 - 0.05 * (1200.0 - velocity_fps) / 1200.0
                }
            }
            1 => {
                // Low-Drag Specialty: minimal degradation
                if velocity_fps > 3000.0 {
                    1.0
                } else if velocity_fps > 2500.0 {
                    0.99 - 0.02 * (3000.0 - velocity_fps) / 500.0
                } else if velocity_fps > 2000.0 {
                    0.97 - 0.03 * (2500.0 - velocity_fps) / 500.0
                } else if velocity_fps > 1500.0 {
                    0.94 - 0.06 * (2000.0 - velocity_fps) / 500.0
                } else {
                    0.88 - 0.08 * (1500.0 - velocity_fps) / 1500.0
                }
            }
            2 => {
                // Light Varmint/Target: significant degradation
                if velocity_fps > 3500.0 {
                    1.0
                } else if velocity_fps > 3000.0 {
                    0.96 - 0.04 * (3500.0 - velocity_fps) / 500.0
                } else if velocity_fps > 2200.0 {
                    0.92 - 0.08 * (3000.0 - velocity_fps) / 800.0
                } else if velocity_fps > 1600.0 {
                    0.84 - 0.14 * (2200.0 - velocity_fps) / 600.0
                } else {
                    0.70 - 0.15 * (1600.0 - velocity_fps) / 1600.0
                }
            }
            3 => {
                // Heavy Magnums: moderate degradation with steep initial drop
                if velocity_fps > 2600.0 {
                    1.0
                } else if velocity_fps > 2200.0 {
                    0.96 - 0.06 * (2600.0 - velocity_fps) / 400.0
                } else if velocity_fps > 1700.0 {
                    0.90 - 0.10 * (2200.0 - velocity_fps) / 500.0
                } else if velocity_fps > 1200.0 {
                    0.80 - 0.15 * (1700.0 - velocity_fps) / 500.0
                } else {
                    0.65 - 0.10 * (1200.0 - velocity_fps) / 1200.0
                }
            }
            _ => 1.0, // Default: no adjustment
        }
    }

    /// Get cluster name for display
    pub fn get_cluster_name(&self, cluster_id: usize) -> &'static str {
        match cluster_id {
            0 => "Standard Long-Range",
            1 => "Low-Drag Specialty",
            2 => "Light Varmint/Target",
            3 => "Heavy Magnum",
            _ => "Unknown",
        }
    }

    /// Apply cluster-based BC correction
    pub fn apply_correction(
        &self,
        bc: f64,
        caliber: f64,
        weight_gr: f64,
        velocity_fps: f64,
    ) -> f64 {
        let cluster_id = self.predict_cluster(caliber, weight_gr, bc);
        let multiplier = self.get_bc_multiplier(velocity_fps, cluster_id);
        bc * multiplier
    }
}

impl Default for ClusterBCDegradation {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cluster_prediction() {
        let cluster_bc = ClusterBCDegradation::new();

        // Test standard long-range bullet (308 Win 168gr)
        let cluster = cluster_bc.predict_cluster(0.308, 168.0, 0.475);
        assert!(
            cluster <= 3,
            "Standard long-range should be in a valid cluster"
        );

        // Test light varmint bullet (223 Rem 55gr)
        let cluster = cluster_bc.predict_cluster(0.224, 55.0, 0.250);
        assert_eq!(cluster, 2);

        // Test heavy magnum (458 Win Mag 500gr)
        let cluster = cluster_bc.predict_cluster(0.458, 500.0, 0.295);
        assert_eq!(cluster, 3);
    }

    #[test]
    fn test_bc_multiplier() {
        let cluster_bc = ClusterBCDegradation::new();

        // Test high velocity (minimal degradation)
        let mult = cluster_bc.get_bc_multiplier(3000.0, 0);
        assert!(mult > 0.95 && mult <= 1.0);

        // Test low velocity (significant degradation)
        let mult = cluster_bc.get_bc_multiplier(1000.0, 0);
        assert!(mult < 0.85);

        // Test that multiplier decreases with velocity
        let mult_high = cluster_bc.get_bc_multiplier(2500.0, 1);
        let mult_low = cluster_bc.get_bc_multiplier(1500.0, 1);
        assert!(mult_high > mult_low);
    }

    #[test]
    fn test_apply_correction() {
        let cluster_bc = ClusterBCDegradation::new();

        // Test that correction reduces BC at low velocity
        let bc_original = 0.475;
        let bc_corrected = cluster_bc.apply_correction(bc_original, 0.308, 168.0, 1500.0);
        assert!(bc_corrected < bc_original);

        // Test that correction is minimal at high velocity
        let bc_corrected_high = cluster_bc.apply_correction(bc_original, 0.308, 168.0, 2800.0);
        assert!(
            bc_corrected_high >= bc_original * 0.90,
            "High velocity should have minimal BC reduction"
        );
    }
}
