use crate::BCSegmentData;

/// Bullet type classification based on model name
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BulletType {
    MatchBoatTail,
    MatchFlatBase,
    HuntingBoatTail,
    HuntingFlatBase,
    VldHighBc,
    Hybrid,
    FMJ,
    RoundNose,
    Unknown,
}

/// BC degradation factors for different bullet types
pub struct BulletTypeFactors {
    pub drop: f64,
    pub transition_curve: f64,
}

impl BulletType {
    /// Get degradation factors for this bullet type
    pub fn get_factors(&self) -> BulletTypeFactors {
        match self {
            BulletType::MatchBoatTail => BulletTypeFactors {
                drop: 0.075,  // 7.5% total drop for match boat tail
                transition_curve: 0.3,
            },
            BulletType::MatchFlatBase => BulletTypeFactors {
                drop: 0.10,  // 10% for match flat base
                transition_curve: 0.35,
            },
            BulletType::HuntingBoatTail => BulletTypeFactors {
                drop: 0.15,  // 15% for hunting boat tail
                transition_curve: 0.45,
            },
            BulletType::HuntingFlatBase => BulletTypeFactors {
                drop: 0.20,  // 20% for hunting flat base
                transition_curve: 0.5,
            },
            BulletType::VldHighBc => BulletTypeFactors {
                drop: 0.05,  // 5% for VLD (very low drag)
                transition_curve: 0.25,
            },
            BulletType::Hybrid => BulletTypeFactors {
                drop: 0.06,  // 6% for hybrid designs
                transition_curve: 0.28,
            },
            BulletType::FMJ => BulletTypeFactors {
                drop: 0.12,  // 12% for military ball
                transition_curve: 0.4,
            },
            BulletType::RoundNose => BulletTypeFactors {
                drop: 0.35,  // 35% for round nose
                transition_curve: 0.7,
            },
            BulletType::Unknown => BulletTypeFactors {
                drop: 0.15,  // Conservative 15%
                transition_curve: 0.5,
            },
        }
    }
}

/// BC segment estimator based on physics and known patterns
pub struct BCSegmentEstimator;

impl BCSegmentEstimator {
    /// Identify bullet type from model name and characteristics
    pub fn identify_bullet_type(model: &str, weight: f64, caliber: f64, bc_value: Option<f64>) -> BulletType {
        let model_lower = model.to_lowercase();
        
        // VLD/High BC bullets
        if model_lower.contains("vld") || model_lower.contains("berger") ||
           model_lower.contains("hybrid") || model_lower.contains("elite") {
            if model_lower.contains("hybrid") {
                return BulletType::Hybrid;
            }
            return BulletType::VldHighBc;
        }
        
        // Match bullets (competition/target)
        if model_lower.contains("smk") || model_lower.contains("matchking") ||
           model_lower.contains("match") || model_lower.contains("bthp") ||
           model_lower.contains("competition") || model_lower.contains("target") ||
           model_lower.contains("a-max") || model_lower.contains("eld-m") ||
           model_lower.contains("scenar") || model_lower.contains("x-ring") {
            // Check for boat tail
            if model_lower.contains("bt") || model_lower.contains("boat") {
                return BulletType::MatchBoatTail;
            }
            // Check if high BC indicates boat tail
            if let Some(bc) = bc_value {
                let sd = Self::calculate_sectional_density(weight, caliber);
                if bc / sd > 1.6 {
                    return BulletType::MatchBoatTail;
                }
            }
            return BulletType::MatchFlatBase;
        }
        
        // Hunting bullets (expanding)
        if model_lower.contains("gameking") || model_lower.contains("hunting") ||
           model_lower.contains("sst") || model_lower.contains("eld-x") ||
           model_lower.contains("partition") || model_lower.contains("accubond") ||
           model_lower.contains("core-lokt") || model_lower.contains("ballistic tip") ||
           model_lower.contains("v-max") || model_lower.contains("hornady sp") ||
           model_lower.contains("interlock") || model_lower.contains("tsx") {
            // Check for boat tail
            if model_lower.contains("bt") || model_lower.contains("boat") ||
               model_lower.contains("sst") || model_lower.contains("accubond") {
                return BulletType::HuntingBoatTail;
            }
            return BulletType::HuntingFlatBase;
        }
        
        // FMJ/Military
        if model_lower.contains("fmj") || model_lower.contains("ball") ||
           model_lower.contains("m80") || model_lower.contains("m855") ||
           model_lower.contains("tracer") {
            return BulletType::FMJ;
        }
        
        // Round nose
        if model_lower.contains("rn") || model_lower.contains("round nose") ||
           model_lower.contains("rnsp") {
            return BulletType::RoundNose;
        }
        
        // Use BC value as hint if available
        if let Some(bc) = bc_value {
            let sd = Self::calculate_sectional_density(weight, caliber);
            let bc_to_sd_ratio = bc / sd;
            
            if bc_to_sd_ratio > 1.8 {
                return BulletType::VldHighBc;
            } else if bc_to_sd_ratio > 1.5 {
                return BulletType::MatchBoatTail;
            } else if bc_to_sd_ratio < 1.2 {
                return BulletType::HuntingFlatBase;
            }
        }
        
        BulletType::Unknown
    }
    
    /// Calculate sectional density (SD) from weight and caliber
    pub fn calculate_sectional_density(weight_grains: f64, caliber_inches: f64) -> f64 {
        // SD = weight / (7000 * caliber^2)
        // Protect against division by zero or negative caliber
        if caliber_inches <= 0.0 {
            return 0.0;
        }
        weight_grains / (7000.0 * caliber_inches * caliber_inches)
    }
    
    /// Estimate BC segments based on bullet characteristics
    pub fn estimate_bc_segments(
        base_bc: f64,
        caliber: f64,
        weight: f64,
        model: &str,
        drag_model: &str,
    ) -> Vec<BCSegmentData> {
        // Identify bullet type
        let bullet_type = Self::identify_bullet_type(model, weight, caliber, Some(base_bc));
        let type_factors = bullet_type.get_factors();
        
        // Calculate sectional density
        let sd = Self::calculate_sectional_density(weight, caliber);
        
        // Adjust BC drop based on sectional density
        // Higher SD = more stable BC
        let sd_factor = (sd / 0.25).max(0.7).min(1.3);
        let adjusted_drop = type_factors.drop / sd_factor;
        
        // Adjust transition curve based on drag model
        let transition_adjustment = if drag_model == "G7" { 0.8 } else { 1.0 };
        let _adjusted_curve = type_factors.transition_curve * transition_adjustment;
        
        // Generate segments based on bullet type
        let mut segments = Vec::new();
        
        // Determine velocity ranges and BC retention factors
        match bullet_type {
            BulletType::MatchBoatTail => {
                // Match boat tail - minimal BC degradation
                segments.push(BCSegmentData {
                    velocity_min: 2800.0,
                    velocity_max: 5000.0,
                    bc_value: base_bc * 1.000,
                });
                segments.push(BCSegmentData {
                    velocity_min: 2400.0,
                    velocity_max: 2800.0,
                    bc_value: base_bc * 0.985,
                });
                segments.push(BCSegmentData {
                    velocity_min: 2000.0,
                    velocity_max: 2400.0,
                    bc_value: base_bc * 0.965,
                });
                segments.push(BCSegmentData {
                    velocity_min: 1600.0,
                    velocity_max: 2000.0,
                    bc_value: base_bc * 0.945,
                });
                segments.push(BCSegmentData {
                    velocity_min: 0.0,
                    velocity_max: 1600.0,
                    bc_value: base_bc * 0.925,
                });
            },
            BulletType::VldHighBc | BulletType::Hybrid => {
                // VLD/Hybrid - very stable BC
                segments.push(BCSegmentData {
                    velocity_min: 2800.0,
                    velocity_max: 5000.0,
                    bc_value: base_bc * 1.000,
                });
                segments.push(BCSegmentData {
                    velocity_min: 2200.0,
                    velocity_max: 2800.0,
                    bc_value: base_bc * 0.990,
                });
                segments.push(BCSegmentData {
                    velocity_min: 1600.0,
                    velocity_max: 2200.0,
                    bc_value: base_bc * 0.970,
                });
                segments.push(BCSegmentData {
                    velocity_min: 0.0,
                    velocity_max: 1600.0,
                    bc_value: base_bc * 0.950,
                });
            },
            BulletType::HuntingBoatTail => {
                // Hunting boat tail - moderate degradation
                segments.push(BCSegmentData {
                    velocity_min: 2600.0,
                    velocity_max: 5000.0,
                    bc_value: base_bc * 1.000,
                });
                segments.push(BCSegmentData {
                    velocity_min: 2200.0,
                    velocity_max: 2600.0,
                    bc_value: base_bc * 0.960,
                });
                segments.push(BCSegmentData {
                    velocity_min: 1800.0,
                    velocity_max: 2200.0,
                    bc_value: base_bc * 0.900,
                });
                segments.push(BCSegmentData {
                    velocity_min: 0.0,
                    velocity_max: 1800.0,
                    bc_value: base_bc * 0.850,
                });
            },
            _ => {
                // Default degradation profile
                segments.push(BCSegmentData {
                    velocity_min: 2800.0,
                    velocity_max: 5000.0,
                    bc_value: base_bc,
                });
                
                let transonic_bc = base_bc * (1.0 - adjusted_drop * 0.3);
                segments.push(BCSegmentData {
                    velocity_min: 1800.0,
                    velocity_max: 2800.0,
                    bc_value: transonic_bc,
                });
                
                let subsonic_bc = base_bc * (1.0 - adjusted_drop);
                segments.push(BCSegmentData {
                    velocity_min: 0.0,
                    velocity_max: 1800.0,
                    bc_value: subsonic_bc,
                });
            }
        }
        
        // Apply sectional density adjustment
        for segment in &mut segments {
            segment.bc_value *= sd_factor.powf(0.5);
            // Ensure we don't exceed nominal BC
            if segment.bc_value > base_bc {
                segment.bc_value = base_bc;
            }
        }
        
        segments
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_bullet_type_identification() {
        assert_eq!(BCSegmentEstimator::identify_bullet_type("168gr SMK", 168.0, 0.308, None), BulletType::MatchFlatBase);
        assert_eq!(BCSegmentEstimator::identify_bullet_type("168gr SMK BT", 168.0, 0.308, None), BulletType::MatchBoatTail);
        assert_eq!(BCSegmentEstimator::identify_bullet_type("150gr SST", 150.0, 0.308, None), BulletType::HuntingBoatTail);
        assert_eq!(BCSegmentEstimator::identify_bullet_type("147gr FMJ", 147.0, 0.308, None), BulletType::FMJ);
        assert_eq!(BCSegmentEstimator::identify_bullet_type("180gr RN", 180.0, 0.308, None), BulletType::RoundNose);
        assert_eq!(BCSegmentEstimator::identify_bullet_type("168gr VLD", 168.0, 0.308, None), BulletType::VldHighBc);
        assert_eq!(BCSegmentEstimator::identify_bullet_type("Some bullet", 150.0, 0.308, None), BulletType::Unknown);
    }
    
    #[test]
    fn test_sectional_density() {
        let sd = BCSegmentEstimator::calculate_sectional_density(168.0, 0.308);
        assert!((sd - 0.253).abs() < 0.001);
    }
    
    #[test]
    fn test_bc_estimation() {
        let segments = BCSegmentEstimator::estimate_bc_segments(
            0.450, 0.308, 168.0, "168gr SMK", "G1"
        );
        
        assert_eq!(segments.len(), 3);
        assert_eq!(segments[0].bc_value, 0.450);  // High velocity
        assert!(segments[1].bc_value < segments[0].bc_value);  // Transonic
        assert!(segments[2].bc_value < segments[1].bc_value);  // Subsonic
    }
}