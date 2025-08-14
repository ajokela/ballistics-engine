/// Drag model enum
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DragModel {
    G1,
    G2,
    G5,
    G6,
    G7,
    G8,
    GI,
    GS,
}

impl DragModel {
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_uppercase().as_str() {
            "G1" => Some(DragModel::G1),
            "G2" => Some(DragModel::G2),
            "G5" => Some(DragModel::G5),
            "G6" => Some(DragModel::G6),
            "G7" => Some(DragModel::G7),
            "G8" => Some(DragModel::G8),
            "GI" => Some(DragModel::GI),
            "GS" => Some(DragModel::GS),
            _ => None,
        }
    }
}

impl std::fmt::Display for DragModel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_drag_model_from_str() {
        // Test valid drag models
        assert_eq!(DragModel::from_str("G1"), Some(DragModel::G1));
        assert_eq!(DragModel::from_str("G2"), Some(DragModel::G2));
        assert_eq!(DragModel::from_str("G5"), Some(DragModel::G5));
        assert_eq!(DragModel::from_str("G6"), Some(DragModel::G6));
        assert_eq!(DragModel::from_str("G7"), Some(DragModel::G7));
        assert_eq!(DragModel::from_str("G8"), Some(DragModel::G8));
        assert_eq!(DragModel::from_str("GI"), Some(DragModel::GI));
        assert_eq!(DragModel::from_str("GS"), Some(DragModel::GS));
    }
    
    #[test]
    fn test_drag_model_from_str_case_insensitive() {
        // Test case insensitivity
        assert_eq!(DragModel::from_str("g1"), Some(DragModel::G1));
        assert_eq!(DragModel::from_str("G1"), Some(DragModel::G1));
        assert_eq!(DragModel::from_str("g7"), Some(DragModel::G7));
        assert_eq!(DragModel::from_str("G7"), Some(DragModel::G7));
        assert_eq!(DragModel::from_str("gi"), Some(DragModel::GI));
        assert_eq!(DragModel::from_str("GI"), Some(DragModel::GI));
        assert_eq!(DragModel::from_str("Gs"), Some(DragModel::GS));
    }
    
    #[test]
    fn test_drag_model_from_str_invalid() {
        // Test invalid inputs
        assert_eq!(DragModel::from_str("G9"), None);
        assert_eq!(DragModel::from_str("G10"), None);
        assert_eq!(DragModel::from_str("X1"), None);
        assert_eq!(DragModel::from_str(""), None);
        assert_eq!(DragModel::from_str("invalid"), None);
        assert_eq!(DragModel::from_str("123"), None);
    }
    
    #[test]
    fn test_drag_model_display() {
        // Test Display implementation
        assert_eq!(format!("{}", DragModel::G1), "G1");
        assert_eq!(format!("{}", DragModel::G2), "G2");
        assert_eq!(format!("{}", DragModel::G5), "G5");
        assert_eq!(format!("{}", DragModel::G6), "G6");
        assert_eq!(format!("{}", DragModel::G7), "G7");
        assert_eq!(format!("{}", DragModel::G8), "G8");
        assert_eq!(format!("{}", DragModel::GI), "GI");
        assert_eq!(format!("{}", DragModel::GS), "GS");
    }
    
    #[test]
    fn test_drag_model_equality() {
        // Test PartialEq implementation
        assert_eq!(DragModel::G1, DragModel::G1);
        assert_eq!(DragModel::G7, DragModel::G7);
        assert_ne!(DragModel::G1, DragModel::G7);
        assert_ne!(DragModel::G5, DragModel::G6);
        
        // Test that from_str produces equal values
        let g1_from_str = DragModel::from_str("G1").unwrap();
        assert_eq!(g1_from_str, DragModel::G1);
    }
    
    #[test]
    fn test_drag_model_clone() {
        // Test Clone implementation
        let original = DragModel::G7;
        let cloned = original.clone();
        assert_eq!(original, cloned);
        
        // Both should be independent
        assert_eq!(format!("{}", original), "G7");
        assert_eq!(format!("{}", cloned), "G7");
    }
    
    #[test]
    fn test_drag_model_copy() {
        // Test Copy implementation
        let original = DragModel::G1;
        let copied = original; // This uses Copy
        let also_copied = original; // Can still use original
        
        assert_eq!(original, copied);
        assert_eq!(original, also_copied);
        assert_eq!(copied, also_copied);
    }
    
    #[test]
    fn test_drag_model_debug() {
        // Test Debug implementation
        assert_eq!(format!("{:?}", DragModel::G1), "G1");
        assert_eq!(format!("{:?}", DragModel::G7), "G7");
        assert_eq!(format!("{:?}", DragModel::GI), "GI");
        assert_eq!(format!("{:?}", DragModel::GS), "GS");
    }
}