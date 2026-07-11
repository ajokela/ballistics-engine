//! Cluster-based BC degradation for improved accuracy
//!
//! This module implements empirical velocity-dependent BC corrections based on
//! bullet clustering. Bullets are classified into 4 clusters from their physical
//! characteristics, and each cluster carries a velocity-keyed multiplier ladder.
//!
//! MBA-1127 (2026-07): this is a faithful port of the Python reference model
//! `ballistics.ml.cluster_bc_degradation.ClusterBCDegradation`
//! (version tag `2026-07-08-radar-recal-mba1126`), which is the single source of
//! truth. The ladders and the classifier below MUST stay byte-for-byte equivalent
//! to that model; `tests::test_parity_with_python_reference` locks the values with
//! goldens generated from the Python implementation.
//!
//! MBA-1126 recalibration notes (why these values differ from the old MBA-645
//! ladders): the multipliers are ratios r(v) = BC_G7(v) / BC_G7_banded_avg computed
//! in form-factor space against the 84-point standard G7 table, from 50 Lapua radar
//! .drg curves + 281 doppler-derived curves. Because the G7 reference curve already
//! carries the transonic drag rise, a G7 BC is nearly FLAT through the supersonic
//! band and actually RISES slightly (~+4-8%) around 1000-1200 fps for match bullets.
//! The old deep transonic dips (0.82-0.90 at 1200-1500) were G1-band shapes
//! misapplied to G7 semantics and over-degraded every transonic trajectory (and the
//! BC5D tables sampled from this model). Multipliers may therefore exceed 1.0.

use crate::DragModel;

/// One velocity band: (velocity_threshold_fps, multiplier). Bands are ordered by
/// DESCENDING threshold, matching the Python reference; `get_bc_multiplier`
/// interpolates linearly between adjacent bands.
type VelocityBand = (f64, f64);

/// Population-average conversion used by the Python source of truth when a G7
/// BC must be classified by the legacy G1-keyed cluster gates. The converted BC
/// is used only to select a cluster; correction multipliers scale the original BC.
const G7_TO_G1_CLASSIFICATION_RATIO: f64 = 1.98;

/// Per-cluster velocity ladders (MBA-1126 radar-recalibrated, G7-referenced).
/// Index = cluster id (0..=3). Kept identical to the Python `velocity_bands`.
const VELOCITY_BANDS: [[VelocityBand; 8]; 4] = [
    // Cluster 0: Standard Long-Range (Sierra MatchKing, etc.; n~248)
    [
        (3500.0, 1.00),
        (3000.0, 1.01),
        (2500.0, 1.01),
        (2000.0, 0.99),
        (1500.0, 0.99),
        (1200.0, 1.03),
        (1000.0, 1.05),
        (0.0, 0.94),
    ],
    // Cluster 1: Low-Drag Specialty (blunt/monolithic + stubby varmint; n~12,
    // shrunk ~25% toward 1.0). G7 BC runs below the banded avg high, rises
    // steeply through transonic.
    [
        (3500.0, 0.90),
        (3000.0, 0.91),
        (2500.0, 0.94),
        (2000.0, 1.02),
        (1500.0, 1.06),
        (1200.0, 1.17),
        (1000.0, 1.40),
        (0.0, 1.14),
    ],
    // Cluster 2: Light Varmint/Target (n~54)
    [
        (3500.0, 1.00),
        (3000.0, 1.00),
        (2500.0, 1.01),
        (2000.0, 0.99),
        (1500.0, 0.985),
        (1200.0, 1.01),
        (1000.0, 1.03),
        (0.0, 0.90),
    ],
    // Cluster 3: Heavy Magnums (n~6, shrunk 50% toward 1.0)
    [
        (3500.0, 1.00),
        (3000.0, 1.01),
        (2500.0, 1.00),
        (2000.0, 0.99),
        (1500.0, 0.99),
        (1200.0, 1.02),
        (1000.0, 1.04),
        (0.0, 0.94),
    ],
];

#[derive(Debug, Clone, Copy, Default)]
pub struct ClusterBCDegradation;

impl ClusterBCDegradation {
    pub fn new() -> Self {
        Self
    }

    /// Predict which cluster a bullet belongs to.
    ///
    /// Rule-based classifier ported verbatim from the Python reference
    /// `predict_cluster`. `caliber` must be in INCHES (the sectional-density and
    /// form-factor thresholds are inch-referenced); passing meters misclassifies.
    pub fn predict_cluster(&self, caliber: f64, weight_gr: f64, bc_g1: f64) -> usize {
        // Sectional density (grains / (7000 * caliber_in^2))
        let sd = if caliber > 0.0 {
            weight_gr / (7000.0 * caliber * caliber)
        } else {
            0.0
        };
        // Approximate form factor
        let form_factor = if sd > 0.0 { bc_g1 / sd } else { 1.5 };

        // Heavy magnum check
        if weight_gr > 400.0 && caliber > 0.35 {
            return 3;
        }
        // Low-drag specialty check (low form factor)
        if form_factor < 1.1 && bc_g1 < 0.2 {
            return 1;
        }
        // Light varmint check
        if weight_gr < 100.0 && caliber < 0.3 {
            return 2;
        }
        // Default to standard long-range
        0
    }

    /// Predict a cluster from a BC expressed in its active reference-model space.
    ///
    /// The legacy classifier is G1-keyed. G7 values are converted to the same
    /// approximate G1 space used by the Python source of truth before applying
    /// those unchanged gates. Other drag models retain the legacy behavior.
    pub fn predict_cluster_for_drag_model(
        &self,
        caliber: f64,
        weight_gr: f64,
        bc: f64,
        drag_model: DragModel,
    ) -> usize {
        let classification_bc_g1 = if drag_model == DragModel::G7 {
            bc * G7_TO_G1_CLASSIFICATION_RATIO
        } else {
            bc
        };
        self.predict_cluster(caliber, weight_gr, classification_bc_g1)
    }

    /// Get the BC multiplier for a given velocity and cluster.
    ///
    /// Mirrors the Python reference: iterate the cluster's bands (descending
    /// threshold); at the first band whose threshold `<=` velocity, return the
    /// band multiplier if it is the top band, otherwise linearly interpolate
    /// between it and the band above. Below the lowest threshold, return the last
    /// band's multiplier. Unknown clusters fall back to cluster 0 (matching Python,
    /// which coerces an unknown id to the standard profile — NOT to 1.0).
    pub fn get_bc_multiplier(&self, velocity_fps: f64, cluster_id: usize) -> f64 {
        let cid = if cluster_id < VELOCITY_BANDS.len() {
            cluster_id
        } else {
            0
        };
        let bands = &VELOCITY_BANDS[cid];

        for i in 0..bands.len() {
            let (v_threshold, multiplier) = bands[i];
            if velocity_fps >= v_threshold {
                if i == 0 {
                    // Above the highest threshold
                    return multiplier;
                }
                // Interpolate between bands[i-1] (high) and bands[i] (low)
                let (v_high, mult_high) = bands[i - 1];
                let (v_low, mult_low) = bands[i];
                let t = (velocity_fps - v_low) / (v_high - v_low);
                return mult_low + t * (mult_high - mult_low);
            }
        }
        // Below the lowest threshold
        bands[bands.len() - 1].1
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

    /// Apply cluster-based BC correction with the legacy G1 classification contract.
    ///
    /// Retained for API compatibility. Use [`Self::apply_correction_for_drag_model`]
    /// when the BC's reference drag model is known.
    pub fn apply_correction(
        &self,
        bc: f64,
        caliber: f64,
        weight_gr: f64,
        velocity_fps: f64,
    ) -> f64 {
        self.apply_correction_for_drag_model(bc, caliber, weight_gr, velocity_fps, DragModel::G1)
    }

    /// Apply cluster-based BC correction in the BC's reference-model space.
    ///
    /// The classifier's thresholds are calibrated for G1 BCs, while the velocity
    /// ladders are ratios that scale the caller's original BC. G7 values are
    /// therefore converted to an approximate G1 equivalent for classification
    /// only; the selected multiplier is still applied to the original G7 value.
    pub fn apply_correction_for_drag_model(
        &self,
        bc: f64,
        caliber: f64,
        weight_gr: f64,
        velocity_fps: f64,
        drag_model: DragModel,
    ) -> f64 {
        let cluster_id = self.predict_cluster_for_drag_model(caliber, weight_gr, bc, drag_model);
        let multiplier = self.get_bc_multiplier(velocity_fps, cluster_id);
        bc * multiplier
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cluster_prediction() {
        let cluster_bc = ClusterBCDegradation::new();

        // Standard long-range bullet (.308 Win 168gr) -> cluster 0
        assert_eq!(cluster_bc.predict_cluster(0.308, 168.0, 0.475), 0);

        // Light varmint bullet (.223 Rem 55gr) -> cluster 2
        assert_eq!(cluster_bc.predict_cluster(0.224, 55.0, 0.250), 2);

        // Heavy magnum (.458 Win Mag 500gr) -> cluster 3
        assert_eq!(cluster_bc.predict_cluster(0.458, 500.0, 0.295), 3);
    }

    #[test]
    fn test_caliber_must_be_inches_not_meters() {
        // predict_cluster's SD/form-factor thresholds are inch-referenced.
        // Passing the caliber in meters (the old cli_api bug, caliber_inches *
        // 0.0254) collapses the sectional density and misclassifies. Guard it.
        let cluster_bc = ClusterBCDegradation::new();
        let inches = cluster_bc.predict_cluster(0.458, 500.0, 0.295);
        let meters = cluster_bc.predict_cluster(0.458 * 0.0254, 500.0, 0.295);
        assert_eq!(inches, 3, ".458/500gr in inches -> Heavy Magnum cluster");
        assert_ne!(
            meters, inches,
            "caliber in meters must misclassify (caller must pass inches)"
        );
    }

    #[test]
    fn test_bc_multiplier_g7_semantics() {
        let cluster_bc = ClusterBCDegradation::new();

        // MBA-1126: G7-referenced multipliers hover around 1.0 and may exceed it.
        // Cluster 0 at 3000 fps is 1.01 (slight rise), not a sub-1.0 degradation.
        let mult_hi = cluster_bc.get_bc_multiplier(3000.0, 0);
        assert!(
            (mult_hi - 1.01).abs() < 1e-9,
            "cluster 0 @3000fps should be 1.01, got {mult_hi}"
        );

        // Transonic RISE for match bullets (the old deep dip is gone): the
        // 1000-1200 fps region is above the mid-supersonic (1500-2000) values.
        let mult_transonic = cluster_bc.get_bc_multiplier(1100.0, 0);
        let mult_mid = cluster_bc.get_bc_multiplier(1750.0, 0);
        assert!(
            mult_transonic > mult_mid,
            "G7 transonic ({mult_transonic}) should now rise above mid-supersonic ({mult_mid})"
        );
    }

    #[test]
    fn test_apply_correction_g7() {
        let cluster_bc = ClusterBCDegradation::new();
        let bc = 0.190;

        // A representative .224 77gr match bullet belongs to light-target
        // cluster 2 when its G7 BC is classified in G1-equivalent space.
        assert_eq!(cluster_bc.predict_cluster(0.224, 77.0, bc), 1);
        let classification_cases = [
            (0.224, 55.0, 0.120, 2),
            (0.224, 69.0, 0.169, 2),
            (0.224, 77.0, bc, 2),
            (0.243, 87.0, 0.190, 2),
            (0.20, 40.0, 0.075, 1),
        ];
        for (caliber, weight, g7_bc, expected) in classification_cases {
            let got = cluster_bc.predict_cluster_for_drag_model(
                caliber,
                weight,
                g7_bc,
                DragModel::G7,
            );
            assert_eq!(
                got,
                expected,
                "G7 classification for ({caliber}, {weight}, {g7_bc})"
            );
        }

        let hi = cluster_bc.apply_correction_for_drag_model(bc, 0.224, 77.0, 2800.0, DragModel::G7);
        assert!(
            (hi / bc - 1.004).abs() < 1e-12,
            "wrong cluster: {}",
            hi / bc
        );

        let transonic =
            cluster_bc.apply_correction_for_drag_model(bc, 0.224, 77.0, 1100.0, DragModel::G7);
        assert!(
            (transonic / bc - 1.02).abs() < 1e-12,
            "wrong transonic cluster: {}",
            transonic / bc
        );

        // The existing API remains byte-identical to explicit G1 classification.
        let legacy = cluster_bc.apply_correction(bc, 0.224, 77.0, 2800.0);
        let explicit_g1 = cluster_bc.apply_correction_for_drag_model(
            bc,
            0.224,
            77.0,
            2800.0,
            DragModel::G1,
        );
        assert_eq!(legacy.to_bits(), explicit_g1.to_bits());
    }

    /// Golden parity test — values generated directly from the Python reference
    /// `ClusterBCDegradation` (tag 2026-07-08-radar-recal-mba1126). If this fails,
    /// the Rust port has drifted from the source of truth; regenerate and reconcile.
    #[test]
    fn test_parity_with_python_reference() {
        let m = ClusterBCDegradation::new();

        // (cluster, velocity_fps, expected_multiplier)
        let mult_goldens: &[(usize, f64, f64)] = &[
            (0, 3600.0, 1.000000), (0, 3500.0, 1.000000), (0, 3250.0, 1.005000),
            (0, 3000.0, 1.010000), (0, 2800.0, 1.010000), (0, 2500.0, 1.010000),
            (0, 2250.0, 1.000000), (0, 2000.0, 0.990000), (0, 1750.0, 0.990000),
            (0, 1500.0, 0.990000), (0, 1350.0, 1.010000), (0, 1200.0, 1.030000),
            (0, 1100.0, 1.040000), (0, 1000.0, 1.050000), (0, 900.0, 1.039000),
            (0, 800.0, 1.028000), (0, 500.0, 0.995000), (0, 0.0, 0.940000),
            (1, 3600.0, 0.900000), (1, 3500.0, 0.900000), (1, 3250.0, 0.905000),
            (1, 3000.0, 0.910000), (1, 2800.0, 0.922000), (1, 2500.0, 0.940000),
            (1, 2250.0, 0.980000), (1, 2000.0, 1.020000), (1, 1750.0, 1.040000),
            (1, 1500.0, 1.060000), (1, 1350.0, 1.115000), (1, 1200.0, 1.170000),
            (1, 1100.0, 1.285000), (1, 1000.0, 1.400000), (1, 900.0, 1.374000),
            (1, 800.0, 1.348000), (1, 500.0, 1.270000), (1, 0.0, 1.140000),
            (2, 3600.0, 1.000000), (2, 3500.0, 1.000000), (2, 3250.0, 1.000000),
            (2, 3000.0, 1.000000), (2, 2800.0, 1.004000), (2, 2500.0, 1.010000),
            (2, 2250.0, 1.000000), (2, 2000.0, 0.990000), (2, 1750.0, 0.987500),
            (2, 1500.0, 0.985000), (2, 1350.0, 0.997500), (2, 1200.0, 1.010000),
            (2, 1100.0, 1.020000), (2, 1000.0, 1.030000), (2, 900.0, 1.017000),
            (2, 800.0, 1.004000), (2, 500.0, 0.965000), (2, 0.0, 0.900000),
            (3, 3600.0, 1.000000), (3, 3500.0, 1.000000), (3, 3250.0, 1.005000),
            (3, 3000.0, 1.010000), (3, 2800.0, 1.006000), (3, 2500.0, 1.000000),
            (3, 2250.0, 0.995000), (3, 2000.0, 0.990000), (3, 1750.0, 0.990000),
            (3, 1500.0, 0.990000), (3, 1350.0, 1.005000), (3, 1200.0, 1.020000),
            (3, 1100.0, 1.030000), (3, 1000.0, 1.040000), (3, 900.0, 1.030000),
            (3, 800.0, 1.020000), (3, 500.0, 0.990000), (3, 0.0, 0.940000),
        ];
        for &(c, v, expected) in mult_goldens {
            let got = m.get_bc_multiplier(v, c);
            assert!(
                (got - expected).abs() < 1e-6,
                "multiplier parity: cluster {c} @{v}fps expected {expected}, got {got}"
            );
        }

        // (caliber_in, weight_gr, bc_g1, expected_cluster)
        let cluster_goldens: &[(f64, f64, f64, usize)] = &[
            (0.308, 168.0, 0.475, 0), (0.308, 168.0, 0.223, 0),
            (0.224, 55.0, 0.250, 2), (0.224, 77.0, 0.372, 2),
            (0.458, 500.0, 0.295, 3), (0.375, 350.0, 0.310, 0),
            (0.264, 140.0, 0.315, 0), (0.284, 180.0, 0.360, 0),
            (0.338, 300.0, 0.400, 0), (0.243, 105.0, 0.278, 0),
            (0.20, 40.0, 0.15, 1), (0.172, 20.0, 0.10, 1),
            (0.510, 750.0, 0.500, 3), (0.30, 150.0, 0.19, 1),
            (0.36, 410.0, 0.40, 3),
        ];
        for &(cal, wt, bc, expected) in cluster_goldens {
            let got = m.predict_cluster(cal, wt, bc);
            assert_eq!(
                got, expected,
                "cluster parity: ({cal}, {wt}, {bc}) expected {expected}, got {got}"
            );
        }
    }
}
