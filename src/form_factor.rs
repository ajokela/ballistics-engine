/// Form factor calculations for drag enhancement
use crate::DragModel;

/// Get default form factor based on bullet type
pub fn get_default_form_factor(bullet_name: &str, drag_model: &DragModel) -> f64 {
    let name_upper = bullet_name.to_uppercase();
    
    match drag_model {
        DragModel::G1 => {
            if name_upper.contains("MATCH") || name_upper.contains("SMK") || 
               name_upper.contains("SCENAR") || name_upper.contains("BERGER") {
                0.90  // Match bullets
            } else if name_upper.contains("VLD") || name_upper.contains("HYBRID") || 
                      name_upper.contains("ELD") {
                0.85  // Very low drag bullets
            } else if name_upper.contains("FMJ") {
                1.10  // Full metal jacket
            } else if name_upper.contains("HUNT") || name_upper.contains("SP") {
                1.05  // Hunting bullets
            } else {
                1.00  // Default
            }
        },
        DragModel::G7 => {
            if name_upper.contains("MATCH") || name_upper.contains("SMK") || 
               name_upper.contains("SCENAR") {
                0.95  // Match bullets
            } else if name_upper.contains("VLD") || name_upper.contains("HYBRID") {
                0.90  // Very low drag bullets
            } else if name_upper.contains("FMJ") {
                1.15  // Full metal jacket (less efficient for G7)
            } else {
                1.05  // Default G7
            }
        },
        _ => 1.00  // Default form factor for other models
    }
}

/// Apply form factor to drag coefficient
pub fn apply_form_factor_to_drag(
    cd_base: f64, 
    bullet_name: Option<&str>, 
    drag_model: &DragModel,
    use_form_factor: bool,
) -> f64 {
    if !use_form_factor {
        return cd_base;
    }
    
    // Get form factor based on bullet type
    let form_factor = if let Some(name) = bullet_name {
        get_default_form_factor(name, drag_model)
    } else {
        1.0  // No correction if no bullet name
    };
    
    // Apply form factor to drag coefficient
    cd_base * form_factor
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_g1_form_factors() {
        let drag_model = DragModel::G1;
        
        // Match bullets should have lower form factor
        assert_eq!(get_default_form_factor("168gr SMK", &drag_model), 0.90);
        assert_eq!(get_default_form_factor("Berger Match", &drag_model), 0.90);
        assert_eq!(get_default_form_factor("Lapua Scenar", &drag_model), 0.90);
        
        // VLD bullets should have even lower
        assert_eq!(get_default_form_factor("VLD Target", &drag_model), 0.85);
        assert_eq!(get_default_form_factor("ELD-X", &drag_model), 0.85);  // ELD without "Match"
        assert_eq!(get_default_form_factor("Hybrid OTM", &drag_model), 0.85);  // Hybrid without "Match"
        
        // FMJ should have higher form factor
        assert_eq!(get_default_form_factor("M855 FMJ", &drag_model), 1.10);
        
        // Hunting bullets
        assert_eq!(get_default_form_factor("Hunting SP", &drag_model), 1.05);
        
        // Default
        assert_eq!(get_default_form_factor("Generic Bullet", &drag_model), 1.00);
    }
    
    #[test]
    fn test_g7_form_factors() {
        let drag_model = DragModel::G7;
        
        // Match bullets
        assert_eq!(get_default_form_factor("175gr SMK", &drag_model), 0.95);
        
        // VLD bullets
        assert_eq!(get_default_form_factor("VLD Hunting", &drag_model), 0.90);
        
        // FMJ less efficient with G7
        assert_eq!(get_default_form_factor("147gr FMJ", &drag_model), 1.15);
        
        // Default G7
        assert_eq!(get_default_form_factor("Generic", &drag_model), 1.05);
    }
    
    #[test]
    fn test_apply_form_factor() {
        let cd_base = 0.5;
        let drag_model = DragModel::G1;
        
        // With form factor enabled
        let cd_with_ff = apply_form_factor_to_drag(
            cd_base, 
            Some("168gr SMK"), 
            &drag_model, 
            true
        );
        assert_eq!(cd_with_ff, 0.5 * 0.90); // 0.45
        
        // With form factor disabled
        let cd_without_ff = apply_form_factor_to_drag(
            cd_base,
            Some("168gr SMK"),
            &drag_model,
            false
        );
        assert_eq!(cd_without_ff, 0.5);
        
        // No bullet name
        let cd_no_name = apply_form_factor_to_drag(
            cd_base,
            None,
            &drag_model,
            true
        );
        assert_eq!(cd_no_name, 0.5); // No correction
    }
    
    #[test]
    fn test_case_insensitive() {
        let drag_model = DragModel::G1;
        
        // Should work with various cases
        assert_eq!(get_default_form_factor("smk", &drag_model), 0.90);
        assert_eq!(get_default_form_factor("SMK", &drag_model), 0.90);
        assert_eq!(get_default_form_factor("SmK", &drag_model), 0.90);
    }
}