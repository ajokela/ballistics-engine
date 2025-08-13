use nalgebra::Vector3;
use crate::InternalBallisticInputs as BallisticInputs;
use crate::constants::*;
use crate::atmosphere::{get_local_atmosphere, get_direct_atmosphere};
use crate::drag::get_drag_coefficient_full;
use crate::spin_drift::{calculate_enhanced_spin_drift, apply_enhanced_spin_drift};
// use crate::bc_estimation::BCSegmentEstimator;
use crate::form_factor::apply_form_factor_to_drag;
// use crate::cluster_bc::ClusterBCDegradation;

// Physics constants
const INCHES_PER_FOOT: f64 = 12.0;
const INCHES_TO_METERS: f64 = 0.0254;
const STANDARD_AIR_DENSITY_METRIC: f64 = 1.225; // kg/m³ at sea level

// Magnus Effect Constants
//
// The Magnus effect causes spinning projectiles to deflect perpendicular to both
// their velocity vector and spin axis due to asymmetric pressure distribution.
// These constants define the Magnus moment coefficient (C_Lα) for different flight regimes.

/// Magnus coefficient for subsonic flow (M < 0.8)
/// 
/// Value: 4.0 (dimensionless coefficient)
/// Physical basis: Fully developed boundary layer circulation around spinning projectile
/// Regime: Subsonic flow where boundary layer remains attached
/// Source: McCoy's "Modern Exterior Ballistics", validated against wind tunnel data
const MAGNUS_COEFF_SUBSONIC: f64 = 4.0;

/// Magnus coefficient reduction factor for transonic regime (0.8 < M < 1.2)
/// 
/// Value: 2.0 (50% reduction from subsonic value)
/// Physical basis: Shock waves disrupt circulation patterns, reducing Magnus effect
/// Effect: Spin drift significantly reduced in transonic flight
/// Source: Experimental spinning projectile studies
const MAGNUS_COEFF_TRANSONIC_REDUCTION: f64 = 2.0;

/// Base Magnus coefficient for supersonic flow (M > 1.2)
/// 
/// Value: 2.0 (dimensionless coefficient)
/// Physical basis: Shock-dominated flow with reduced but persistent circulation
/// Effect: Lower Magnus effect than subsonic, but higher than transonic minimum
const MAGNUS_COEFF_SUPERSONIC_BASE: f64 = 2.0;

/// Magnus coefficient scaling factor for high supersonic speeds
/// 
/// Value: 1.5 (additional scaling with Mach number)
/// Formula: Magnus_coeff = BASE + SCALE * (M - 1.2) for M > 1.2
/// Physical basis: Partial recovery of circulation effects at higher Mach numbers
const MAGNUS_COEFF_SUPERSONIC_SCALE: f64 = 1.5;

/// Transonic regime boundaries for Magnus effect calculations
const MAGNUS_TRANSONIC_LOWER: f64 = 0.8;  // Lower bound of transonic regime
const MAGNUS_TRANSONIC_UPPER: f64 = 1.2;  // Upper bound of transonic regime
const MAGNUS_TRANSONIC_RANGE: f64 = 0.4;  // Range width (1.2 - 0.8)
const MAGNUS_SUPERSONIC_RANGE: f64 = 1.8; // Scaling range for supersonic recovery

// Note: Previous implementation had Magnus coefficients ~100x too small, which 
// caused compensating factors elsewhere to make the Magnus force ~50x too large.
// These corrected values are calibrated against real-world spin drift measurements.

// Atmosphere detection thresholds
const MAX_REALISTIC_DENSITY: f64 = 2.0; // kg/m³
const MIN_REALISTIC_SPEED_OF_SOUND: f64 = 200.0; // m/s

/// Calculate spin rate from twist rate and velocity
fn calculate_spin_rate(twist_rate: f64, velocity_mps: f64) -> f64 {
    if twist_rate <= 0.0 {
        return 0.0;
    }
    
    // Convert velocity to ft/s and twist rate to ft/turn
    let velocity_fps = velocity_mps * MPS_TO_FPS;
    let twist_rate_ft = twist_rate / INCHES_PER_FOOT;
    
    // Calculate spin rate: revolutions per second = velocity_fps / twist_rate_ft
    // Convert to rad/s: rad/s = (revolutions/s) * 2π
    let revolutions_per_second = velocity_fps / twist_rate_ft;
    
    
    revolutions_per_second * 2.0 * std::f64::consts::PI
}

/// Calculate Magnus moment coefficient C_Lα based on Mach number
/// Based on McCoy's 'Modern Exterior Ballistics' and empirical data
fn calculate_magnus_moment_coefficient(mach: f64) -> f64 {
    // Magnus moment coefficient varies with Mach number
    // Values based on empirical data for spitzer bullets
    
    if mach < MAGNUS_TRANSONIC_LOWER {
        // Subsonic: relatively constant
        MAGNUS_COEFF_SUBSONIC
    } else if mach < MAGNUS_TRANSONIC_UPPER {
        // Transonic: reduced due to shock formation
        // Linear interpolation through transonic region
        MAGNUS_COEFF_SUBSONIC - MAGNUS_COEFF_TRANSONIC_REDUCTION * (mach - MAGNUS_TRANSONIC_LOWER) / MAGNUS_TRANSONIC_RANGE
    } else {
        // Supersonic: gradually recovers
        MAGNUS_COEFF_SUPERSONIC_BASE + MAGNUS_COEFF_SUPERSONIC_SCALE * ((mach - MAGNUS_TRANSONIC_UPPER) / MAGNUS_SUPERSONIC_RANGE).min(1.0)
    }
}

/// Compute ballistic derivatives for trajectory integration
pub fn compute_derivatives(
    pos: Vector3<f64>,
    vel: Vector3<f64>,
    inputs: &BallisticInputs,
    wind_vector: Vector3<f64>,
    atmos_params: (f64, f64, f64, f64),
    bc_used: f64,
    omega_vector: Option<Vector3<f64>>,
    time: f64,
) -> [f64; 6] {
    // Gravity acceleration vector
    let accel_gravity = Vector3::new(0.0, -G_ACCEL_MPS2, 0.0);
    
    // Wind-adjusted velocity
    let velocity_adjusted = vel - wind_vector;
    let speed_air = velocity_adjusted.norm();
    
    // Initialize drag acceleration
    let mut accel_drag = Vector3::zeros();
    let mut accel_magnus = Vector3::zeros();
    
    // Calculate drag if velocity is significant
    if speed_air > crate::constants::MIN_VELOCITY_THRESHOLD {
        let v_rel_fps = speed_air * MPS_TO_FPS;
        
        // Get atmospheric conditions
        let altitude_at_pos = inputs.altitude + pos[1];
        
        // Check if we have direct atmosphere values
        // Direct atmosphere is indicated by having only 2 parameters where:
        // params[0] = air density, params[1] = speed of sound
        // params[2] and params[3] would be 0.0
        // BUT: we need to check if params[0] is a reasonable density value (< 2.0 kg/m³)
        let (air_density, speed_of_sound) = if atmos_params.0 < MAX_REALISTIC_DENSITY && atmos_params.1 > MIN_REALISTIC_SPEED_OF_SOUND && atmos_params.2 == 0.0 && atmos_params.3 == 0.0 {
            // Direct atmosphere values
            get_direct_atmosphere(atmos_params.0, atmos_params.1)
        } else {
            // Calculate from base parameters
            get_local_atmosphere(
                altitude_at_pos,
                atmos_params.0,  // base_alt
                atmos_params.1,  // base_temp_c
                atmos_params.2,  // base_press_hpa
                atmos_params.3,  // base_ratio
            )
        };
        
        // Calculate Mach number with safe division
        let mach = if speed_of_sound > 1e-9 {
            speed_air / speed_of_sound
        } else {
            0.0 // No meaningful Mach number at zero speed of sound
        };
        
        // Get drag coefficient with transonic and Reynolds corrections
        let mut drag_factor = get_drag_coefficient_full(
            mach, 
            &inputs.bc_type,
            true, // apply transonic correction
            true, // apply Reynolds correction
            None, // let it determine shape
            if inputs.caliber_inches > 0.0 { Some(inputs.caliber_inches) } else { Some(inputs.bullet_diameter) },
            if inputs.weight_grains > 0.0 { Some(inputs.weight_grains) } else { Some(inputs.bullet_mass * 15.432358) },
            Some(speed_air),
            Some(air_density),
            Some(atmos_params.1), // temperature in Celsius
        );
        
        // Apply form factor if enabled
        if inputs.use_form_factor {
            drag_factor = apply_form_factor_to_drag(
                drag_factor,
                inputs.bullet_model.as_deref(),
                &inputs.bc_type,
                true,
            );
        }
        
        // Get BC value - check if cluster BC should override segments
        let mut bc_val = bc_used;
        
        // If cluster BC is enabled, skip BC segments to avoid double reduction
        if inputs.use_cluster_bc {
            // Apply cluster-based BC degradation
            let cluster_bc = ClusterBCDegradation::new();
            bc_val = cluster_bc.adjust_bc(
                bc_used, 
                v_rel_fps, 
                inputs.bullet_diameter, 
                inputs.bullet_mass,
                inputs.bullet_cluster
            );
        } else if inputs.use_bc_segments {
            // Use velocity-based segments if cluster BC is disabled
            bc_val = get_bc_for_velocity(v_rel_fps, inputs, bc_used);
        } else if let Some(ref segments) = inputs.bc_segments {
            // Fall back to Mach-based segments
            bc_val = interpolated_bc(mach, segments, Some(inputs));
        }
        
        // Calculate yaw effect with safe division
        let yaw_deg = if inputs.tipoff_decay_distance.abs() > 1e-9 {
            inputs.tipoff_yaw * (-pos[0] / inputs.tipoff_decay_distance).exp()
        } else {
            inputs.tipoff_yaw // No decay if distance is zero
        };
        let yaw_rad = yaw_deg.to_radians();
        let yaw_multiplier = 1.0 + yaw_rad.powi(2);
        
        // Calculate density scaling
        let density_scale = air_density / STANDARD_AIR_DENSITY;
        
        // Apply transonic correction if in transonic regime
        let drag_factor = if mach > 0.7 && mach < 1.3 {
            // Estimate projectile shape from parameters
            let shape = crate::transonic_drag::get_projectile_shape(
                inputs.bullet_diameter,
                inputs.bullet_mass,
                &inputs.bc_type.to_string(),
            );
            
            // Apply transonic correction
            let corrected_cd = crate::transonic_drag::transonic_correction(
                mach,
                drag_factor,
                shape,
                true, // include_wave_drag
            );
            
            // The drag_factor is already a coefficient, so we need to calculate the correction ratio
            corrected_cd / drag_factor
        } else {
            1.0
        } * drag_factor;
        
        // Apply Reynolds correction for low velocities
        let drag_factor = if mach < 1.0 && speed_air < 200.0 {
            // Get temperature from atmospheric parameters
            let temperature_c = atmos_params.1;  // base_temp_c from atmos_params
            
            // Apply Reynolds number correction
            
            crate::reynolds::apply_reynolds_correction(
                drag_factor,
                speed_air,
                inputs.bullet_diameter,
                air_density,
                temperature_c,
                mach,
            )
        } else {
            drag_factor
        };
        
        // Calculate drag acceleration
        let standard_factor = drag_factor * CD_TO_RETARD;
        let a_drag_ft_s2 = (v_rel_fps.powi(2) * standard_factor * yaw_multiplier * density_scale) / bc_val;
        let a_drag_m_s2 = a_drag_ft_s2 * FPS_TO_MPS;
        
        // Apply drag in opposite direction of relative velocity
        accel_drag = -a_drag_m_s2 * (velocity_adjusted / speed_air);
        
        // Magnus Effect calculation
        if inputs.enable_advanced_effects && inputs.bullet_diameter > 0.0 && inputs.twist_rate > 0.0 {
            
            // Calculate spin rate from twist rate and velocity
            let spin_rate_rad_s = calculate_spin_rate(inputs.twist_rate, speed_air);
            
            // Calculate Magnus moment coefficient
            let c_la = calculate_magnus_moment_coefficient(mach);
            
            // Convert diameter to meters
            let diameter_m = inputs.bullet_diameter * INCHES_TO_METERS;
            
            // Magnus force formula for spinning projectiles
            // Based on McCoy's Modern Exterior Ballistics
            // F_Magnus = C_L × ½ × ρ × V² × A
            // where C_L = C_Lα × spin_parameter
            
            // Calculate spin parameter (dimensionless) with safe division
            let spin_param = if speed_air > 1e-9 {
                spin_rate_rad_s * diameter_m / (2.0 * speed_air)
            } else {
                0.0 // No spin effect at zero speed
            };
            
            // Calculate lift coefficient
            let c_l = spin_param * c_la;
            
            // Calculate reference area
            let area = std::f64::consts::PI * (diameter_m / 2.0).powi(2);
            
            // Calculate Magnus force using standard lift equation
            // F = 0.5 * ρ * V² * A * C_L
            // Apply empirical calibration factor to match real-world data
            // This accounts for the fact that bullets are not perfect cylinders
            // and have complex flow patterns
            const MAGNUS_CALIBRATION_FACTOR: f64 = 1.8; // Calibrated to produce 4-6 inches drift at 200 yards
            let magnus_force_magnitude = MAGNUS_CALIBRATION_FACTOR * 0.5 * air_density * speed_air.powi(2) * area * c_l;
            
            // Magnus force is perpendicular to both velocity and spin axis
            // For a bullet spinning around its axis of travel, the spin vector is aligned with velocity
            let velocity_unit = velocity_adjusted / speed_air;
            
            // The Magnus force creates lift perpendicular to velocity
            // For right-hand twist, force is to the right when looking downrange
            // We need a vector perpendicular to velocity in the horizontal plane
            
            // Simplified approach: Magnus primarily causes horizontal drift
            // The force is perpendicular to both spin axis (velocity) and gravity
            let vertical = Vector3::new(0.0, 1.0, 0.0); // Up direction
            
            // Magnus force direction: velocity × vertical (for right-hand twist)
            let magnus_direction = velocity_unit.cross(&vertical);
            let magnus_norm = magnus_direction.norm();
            
            if magnus_norm > 1e-12 && magnus_force_magnitude > 1e-12 {
                let magnus_direction = magnus_direction / magnus_norm;
                
                // Reverse direction for left-hand twist
                let magnus_direction = if inputs.is_twist_right {
                    magnus_direction
                } else {
                    -magnus_direction
                };
                
                // Convert bullet mass to kg
                let bullet_mass_kg = inputs.bullet_mass * GRAINS_TO_KG;
                
                // Calculate acceleration
                accel_magnus = (magnus_force_magnitude / bullet_mass_kg) * magnus_direction;
            }
        }
    }
    
    // Total acceleration
    let mut accel = accel_gravity + accel_drag + accel_magnus;
    
    // Add Coriolis acceleration if omega vector is provided
    if let Some(omega) = omega_vector {
        let accel_coriolis = -2.0 * omega.cross(&vel);
        accel += accel_coriolis;
    }
    
    // Apply enhanced spin drift if enabled
    let mut derivatives = [vel[0], vel[1], vel[2], accel[0], accel[1], accel[2]];
    
    if inputs.use_enhanced_spin_drift && inputs.enable_advanced_effects && time > 0.0 {
        // Calculate crosswind component
        let velocity_adjusted = vel - wind_vector;
        let crosswind_speed = if velocity_adjusted.norm() > crate::constants::MIN_VELOCITY_THRESHOLD {
            let trajectory_unit = velocity_adjusted / velocity_adjusted.norm();
            let crosswind = wind_vector - wind_vector.dot(&trajectory_unit) * trajectory_unit;
            crosswind.norm()
        } else {
            0.0
        };
        
        // Get air density (already calculated above)
        let air_density = if speed_air > crate::constants::MIN_VELOCITY_THRESHOLD {
            let altitude_at_pos = inputs.altitude + pos[1];
            let (density, _) = if atmos_params.0 < MAX_REALISTIC_DENSITY && atmos_params.1 > MIN_REALISTIC_SPEED_OF_SOUND && atmos_params.2 == 0.0 && atmos_params.3 == 0.0 {
                get_direct_atmosphere(atmos_params.0, atmos_params.1)
            } else {
                get_local_atmosphere(
                    altitude_at_pos,
                    atmos_params.0,
                    atmos_params.1,
                    atmos_params.2,
                    atmos_params.3,
                )
            };
            density
        } else {
            STANDARD_AIR_DENSITY_METRIC // Standard air density
        };
        
        // Calculate enhanced spin drift components
        let spin_components = calculate_enhanced_spin_drift(
            inputs.bullet_mass,
            vel.norm(),
            inputs.twist_rate,
            inputs.bullet_diameter,
            inputs.bullet_length,
            inputs.is_twist_right,
            time,
            air_density,
            crosswind_speed,
            0.0,   // pitch_rate_rad_s - we don't track angular rates yet
            false, // use_pitch_damping - disabled for now
        );
        
        // Apply enhanced spin drift acceleration
        apply_enhanced_spin_drift(&mut derivatives, &spin_components, time, inputs.is_twist_right);
    }
    
    // Return state derivatives: [velocity, acceleration]
    derivatives
}

/// Calculate appropriate BC fallback based on available bullet parameters
fn calculate_bc_fallback(
    bullet_mass: Option<f64>,      // grains
    bullet_diameter: Option<f64>,  // inches  
    bc_type: Option<&str>          // "G1" or "G7"
) -> f64 {
    use crate::constants::*;
    
    // Weight-based fallback (most reliable predictor)
    if let Some(weight) = bullet_mass {
        let base_bc = if weight < 50.0 {
            BC_FALLBACK_ULTRA_LIGHT
        } else if weight < 100.0 {
            BC_FALLBACK_LIGHT
        } else if weight < 150.0 {
            BC_FALLBACK_MEDIUM
        } else if weight < 200.0 {
            BC_FALLBACK_HEAVY
        } else {
            BC_FALLBACK_VERY_HEAVY
        };
        
        // G7 vs G1 adjustment
        return if let Some(drag_model) = bc_type {
            if drag_model == "G7" {
                base_bc * 0.85 // G7 BCs are typically lower than G1
            } else {
                base_bc
            }
        } else {
            base_bc
        };
    }
    
    // Caliber-based fallback (second most reliable)
    if let Some(caliber) = bullet_diameter {
        let base_bc = if caliber <= 0.224 {
            BC_FALLBACK_SMALL_CALIBER
        } else if caliber <= 0.243 {
            BC_FALLBACK_MEDIUM_CALIBER
        } else if caliber <= 0.284 {
            BC_FALLBACK_LARGE_CALIBER
        } else {
            BC_FALLBACK_XLARGE_CALIBER
        };
        
        // G7 vs G1 adjustment
        return if let Some(drag_model) = bc_type {
            if drag_model == "G7" {
                base_bc * 0.85 // G7 BCs are typically lower than G1
            } else {
                base_bc
            }
        } else {
            base_bc
        };
    }
    
    // Final fallback - conservative overall
    let base_fallback = BC_FALLBACK_CONSERVATIVE;
    if let Some(drag_model) = bc_type {
        if drag_model == "G7" {
            return base_fallback * 0.85;
        }
    }
    
    base_fallback
}

/// Interpolate ballistic coefficient from segments with dynamic fallback
pub fn interpolated_bc(mach: f64, segments: &[(f64, f64)], inputs: Option<&BallisticInputs>) -> f64 {
    if segments.is_empty() {
        // Use dynamic fallback based on bullet characteristics if available
        if let Some(inputs) = inputs {
            let bc_type_str = match inputs.bc_type {
                crate::DragModel::G1 => "G1",
                crate::DragModel::G7 => "G7",
                _ => "G1", // Default to G1 for other models
            };
            return calculate_bc_fallback(
                Some(inputs.bullet_mass),
                Some(inputs.bullet_diameter), 
                Some(bc_type_str)
            );
        }
        return crate::constants::BC_FALLBACK_CONSERVATIVE; // Conservative fallback based on database analysis
    }
    
    if segments.len() == 1 {
        return segments[0].1;
    }
    
    // Sort segments by Mach number to ensure proper interpolation
    let mut sorted_segments = segments.to_vec();
    sorted_segments.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    
    // Handle out-of-range cases first
    if mach <= sorted_segments[0].0 {
        return sorted_segments[0].1;
    }
    if mach >= sorted_segments[sorted_segments.len() - 1].0 {
        return sorted_segments[sorted_segments.len() - 1].1;
    }
    
    // Find the appropriate segment using binary search
    let idx = sorted_segments.partition_point(|(m, _)| *m <= mach);
    if idx == 0 || idx >= sorted_segments.len() {
        // Should not happen given the checks above
        return sorted_segments[0].1;
    }
    
    let (mach1, bc1) = sorted_segments[idx - 1];
    let (mach2, bc2) = sorted_segments[idx];
    
    // Linear interpolation with safe division
    let denominator = mach2 - mach1;
    if denominator.abs() < crate::constants::MIN_DIVISION_THRESHOLD {
        return bc1; // Return first BC value if Mach values are identical
    }
    let t = (mach - mach1) / denominator;
    bc1 + t * (bc2 - bc1)
}

/// Get BC value for current velocity, supporting velocity-based BC segments
fn get_bc_for_velocity(velocity_fps: f64, inputs: &BallisticInputs, bc_used: f64) -> f64 {
    
    // Check if velocity-based BC segments are enabled
    if !inputs.use_bc_segments {
        return bc_used;
    }
    
    // Try direct BC segments data first
    if let Some(ref bc_segments_data) = inputs.bc_segments_data {
        for segment in bc_segments_data {
            if velocity_fps >= segment.velocity_min && velocity_fps <= segment.velocity_max {
                return segment.bc_value;
            }
        }
    }
    
    // Try BC estimation if we have bullet details but no segments
    if inputs.bullet_diameter > 0.0 && inputs.bullet_mass > 0.0 && bc_used > 0.0 {
        // Create a model string from bullet_id or generate generic description
        let model = if let Some(ref bullet_id) = inputs.bullet_id {
            bullet_id.clone()
        } else {
            format!("{}gr bullet", inputs.bullet_mass as i32)
        };
        
        // Estimate segments based on bullet characteristics
        let segments = BCSegmentEstimator::estimate_bc_segments(
            bc_used,
            inputs.bullet_diameter,
            inputs.bullet_mass,
            &model,
            inputs.bc_type_str(),
        );
        
        // Find appropriate segment for current velocity
        for segment in &segments {
            if velocity_fps >= segment.velocity_min && velocity_fps <= segment.velocity_max {
                return segment.bc_value;
            }
        }
    }
    
    // Fallback to constant BC
    bc_used
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DragModel;

    fn create_test_inputs() -> BallisticInputs {
        BallisticInputs {
            bc_value: 0.5,
            bc_type: DragModel::G1,
            bullet_mass: 168.0,
            altitude: 1000.0,
            tipoff_yaw: 0.0,
            tipoff_decay_distance: 20.0,
            ground_threshold: 0.0,
            bc_segments: None,
            muzzle_velocity: 0.0,
            twist_rate: 0.0,
            bullet_length: 0.0,
            bullet_diameter: 0.0,
            target_distance: 0.0,
            muzzle_angle: 0.0,
            use_cluster_bc: false,
            bullet_cluster: None,
            wind_speed: 0.0,
            wind_angle: 0.0,
            temperature: 15.0,
            pressure: 1013.25,
            humidity: 0.0,
            latitude: None,
            enable_advanced_effects: false,
            is_twist_right: true,
            shooting_angle: 0.0,
            use_powder_sensitivity: false,
            powder_temp_sensitivity: 0.0,
            powder_temp: 70.0,
            caliber_inches: 0.0,
            weight_grains: 0.0,
            use_bc_segments: false,
            bullet_id: None,
            bc_segments_data: None,
            use_enhanced_spin_drift: false,
            use_form_factor: false,
            manufacturer: None,
            bullet_model: None,
            enable_wind_shear: false,
            wind_shear_model: "none".to_string(),
        }
    }

    #[test]
    fn test_compute_derivatives_basic() {
        let pos = Vector3::new(0.0, 0.0, 0.0);
        let vel = Vector3::new(800.0, 0.0, 0.0);
        let inputs = create_test_inputs();
        let wind_vector = Vector3::zeros();
        let atmos_params = (288.15, 1013.25, 50.0, 0.0); // Standard conditions
        let bc_used = 0.5;
        
        let result = compute_derivatives(
            pos, vel, &inputs, wind_vector, atmos_params, bc_used, None, 0.0
        ).unwrap();
        
        // Check that we get velocity and acceleration components
        assert_eq!(result.len(), 6);
        
        // Velocity components should match input velocity
        assert!((result[0] - vel[0]).abs() < 1e-10);
        assert!((result[1] - vel[1]).abs() < 1e-10);
        assert!((result[2] - vel[2]).abs() < 1e-10);
        
        // Should have gravitational acceleration
        assert!(result[4] < 0.0); // Negative y acceleration due to gravity
        
        // Should have drag acceleration opposing motion
        assert!(result[3] < 0.0); // Negative x acceleration due to drag
    }

    #[test]
    fn test_compute_derivatives_with_wind() {
        let pos = Vector3::new(0.0, 0.0, 0.0);
        let vel = Vector3::new(800.0, 0.0, 0.0);
        let inputs = create_test_inputs();
        let wind_vector = Vector3::new(10.0, 0.0, 0.0); // Tailwind
        let atmos_params = (288.15, 1013.25, 50.0, 0.0);
        let bc_used = 0.5;
        
        let result = compute_derivatives(
            pos, vel, &inputs, wind_vector, atmos_params, bc_used, None, 0.0
        ).unwrap();
        
        // With tailwind, effective velocity should be lower, thus less drag
        assert!(result[3] > -100.0); // Less negative drag than without wind
    }

    #[test]
    fn test_compute_derivatives_with_coriolis() {
        let pos = Vector3::new(0.0, 0.0, 0.0);
        let vel = Vector3::new(800.0, 0.0, 0.0);
        let inputs = create_test_inputs();
        let wind_vector = Vector3::zeros();
        let atmos_params = (288.15, 1013.25, 50.0, 0.0);
        let bc_used = 0.5;
        let omega = Vector3::new(0.0, 0.0, 7.2921e-5); // Earth's rotation
        
        let result = compute_derivatives(
            pos, vel, &inputs, wind_vector, atmos_params, bc_used, Some(omega), 0.0
        ).unwrap();
        
        // Should have Coriolis effect
        assert!(result[4].abs() > 1e-3); // Should have some y-component from Coriolis
    }

    #[test]
    fn test_interpolated_bc() {
        let segments = vec![
            (0.5, 0.4),
            (1.0, 0.5),
            (1.5, 0.6),
            (2.0, 0.5),
        ];
        
        // Test exact matches
        assert!((interpolated_bc(1.0, &segments, None) - 0.5).abs() < 1e-10);
        
        // Test interpolation
        let bc_075 = interpolated_bc(0.75, &segments, None);
        assert!(bc_075 > 0.4 && bc_075 < 0.5);
        
        // Test out of range
        assert!((interpolated_bc(0.1, &segments, None) - 0.4).abs() < 1e-10);
        assert!((interpolated_bc(3.0, &segments, None) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_interpolated_bc_edge_cases() {
        // Empty segments
        assert!((interpolated_bc(1.0, &[], None) - crate::constants::BC_FALLBACK_CONSERVATIVE).abs() < 1e-10);
        
        // Single segment
        let single = vec![(1.0, 0.7)];
        assert!((interpolated_bc(1.5, &single, None) - 0.7).abs() < 1e-10);
    }

    #[test]
    fn test_magnus_effect() {
        let pos = Vector3::new(0.0, 0.0, 0.0);
        let vel = Vector3::new(822.96, 0.0, 0.0); // 2700 fps
        let mut inputs = create_test_inputs();
        inputs.bullet_diameter = 0.308;  // .308 caliber
        inputs.twist_rate = 10.0;        // 1:10 twist
        inputs.is_twist_right = true;
        
        let wind_vector = Vector3::zeros();
        let atmos_params = (288.15, 1013.25, 50.0, 0.0);
        let bc_used = 0.5;
        
        let result = compute_derivatives(
            pos, vel, &inputs, wind_vector, atmos_params, bc_used, None, 0.0
        ).unwrap();
        
        // Should have Magnus drift in z direction
        assert!(result[5].abs() > 0.01); // Should have some z-acceleration
        assert!(result[5] > 0.0); // Should be positive (right drift) for right-hand twist
        
        // Magnus acceleration should be reasonable (typical: 0.1-1.0 m/s²)
        assert!(result[5] < 2.0); // Less than 2 m/s²
    }

    #[test]
    fn test_magnus_moment_coefficient() {
        // Test at various Mach numbers
        assert!((calculate_magnus_moment_coefficient(0.5) - 0.030).abs() < 1e-6);
        assert!((calculate_magnus_moment_coefficient(0.8) - 0.030).abs() < 1e-6);
        assert!((calculate_magnus_moment_coefficient(1.0) - 0.0225).abs() < 1e-6);
        assert!((calculate_magnus_moment_coefficient(1.2) - 0.015).abs() < 1e-6);
        assert!((calculate_magnus_moment_coefficient(2.0) - 0.0194).abs() < 0.001);
    }
}