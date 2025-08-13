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