//! Turret adjustment-unit conversions and click-value parsing (MBA-1355).
//!
//! Shared by the CLI (`main.rs`) and the WASM terminal (`wasm.rs`). This lives in the
//! library crate — not `main.rs` — so it compiles for `wasm32-unknown-unknown` with no
//! feature gate. To avoid a circular dependency on `main.rs`'s clap `AdjustmentUnit`
//! `ValueEnum`, this module defines its own minimal angular base (`ClickBase`) and factor
//! table; `main.rs` maps its own `AdjustmentUnit` onto `ClickBase`/`adjustment_factor`
//! locally (see `drop_to_adjustment` in `main.rs`).

/// Angular base unit for a turret click graduation (MBA-1355). A click graduation is
/// always expressed in mil, (true) MOA, or SMOA/IPHY — never in whole clicks itself.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ClickBase {
    Mil,
    Moa,
    Smoa,
}

/// The `(drop_yd / range_yd) * factor` conversion factor for a `ClickBase`:
/// - `Mil`: 1000.0 (milliradians)
/// - `Moa`: 3438.0 (true MOA; this crate's locked printed-table constant — see MBA-724 —
///   deliberately not the exact-angle 3437.7467)
/// - `Smoa`: 3600.0 ("shooter's MOA" / inches-per-hundred-yards; exact by definition)
pub fn adjustment_factor(base: ClickBase) -> f64 {
    match base {
        ClickBase::Mil => 1000.0,
        ClickBase::Moa => 3438.0,
        ClickBase::Smoa => 3600.0,
    }
}

/// A turret click graduation, parsed from suffixed CLI/profile syntax
/// like "0.1mil" / "0.25moa" / "0.125smoa" (MBA-1355).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClickValue {
    pub size: f64,
    pub base: ClickBase,
}

/// Parses a suffixed click graduation string. The unit suffix is mandatory — one of
/// `mil`, `moa`, `smoa`, `iphy` (`iphy` is accepted as an alias for `smoa`, the identical
/// unit) — and the magnitude must be a positive, finite number.
pub fn parse_click_value(s: &str) -> Result<ClickValue, String> {
    let t = s.trim().to_ascii_lowercase();
    let (num, base) = if let Some(n) = t.strip_suffix("smoa") {
        (n, ClickBase::Smoa)
    } else if let Some(n) = t.strip_suffix("iphy") {
        (n, ClickBase::Smoa)
    } else if let Some(n) = t.strip_suffix("moa") {
        (n, ClickBase::Moa)
    } else if let Some(n) = t.strip_suffix("mil") {
        (n, ClickBase::Mil)
    } else {
        return Err(format!(
            "click value '{s}' needs a unit suffix: mil, moa, smoa, or iphy (e.g. 0.1mil, 0.25moa)"
        ));
    };
    let size: f64 = num
        .trim()
        .parse()
        .map_err(|_| format!("click value '{s}' has an invalid number '{num}'"))?;
    if !size.is_finite() || size <= 0.0 {
        return Err(format!("click value '{s}' must be a positive, finite graduation"));
    }
    Ok(ClickValue { size, base })
}

/// Whole-click adjustment for a drop at a range: the angular value in the graduation's
/// own base unit divided by the graduation, rounded to the nearest click (ties away from
/// zero). Sign is preserved. Ranges under 1 yard are defined as zero adjustment, matching
/// the CLI's `drop_to_adjustment` short-range guard.
pub fn clicks_for(drop_yd: f64, range_yd: f64, click: &ClickValue) -> i64 {
    if range_yd < 1.0 {
        return 0;
    }
    let angle = (drop_yd / range_yd) * adjustment_factor(click.base);
    (angle / click.size).round() as i64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_click_value_accepts_suffixed_and_rejects_bare() {
        let c = parse_click_value("0.25moa").unwrap();
        assert!((c.size - 0.25).abs() < 1e-12);
        assert!(matches!(c.base, ClickBase::Moa));
        assert!(matches!(parse_click_value("0.1mil").unwrap().base, ClickBase::Mil));
        assert!(matches!(parse_click_value("0.125smoa").unwrap().base, ClickBase::Smoa));
        assert!(matches!(parse_click_value("0.125iphy").unwrap().base, ClickBase::Smoa));
        for bad in ["0.25", "moa", "-0.1mil", "0mil", "0.1mils", ""] {
            assert!(parse_click_value(bad).is_err(), "{bad:?} must be rejected");
        }
        // error message names the accepted suffixes
        let e = parse_click_value("0.25").unwrap_err();
        assert!(e.contains("mil") && e.contains("moa"), "{e}");
    }

    #[test]
    fn clicks_round_to_nearest_whole_click() {
        let quarter_moa = parse_click_value("0.25moa").unwrap();
        // 2.6 TMOA of drop -> 10.4 clicks -> 10
        let drop_yd = 2.6 / 3438.0 * 100.0;
        assert_eq!(clicks_for(drop_yd, 100.0, &quarter_moa), 10);
        // negative (wind left) rounds symmetrically
        assert_eq!(clicks_for(-drop_yd, 100.0, &quarter_moa), -10);
        // 10.5 clicks rounds away from zero
        let drop_yd = 2.625 / 3438.0 * 100.0;
        assert_eq!(clicks_for(drop_yd, 100.0, &quarter_moa), 11);
    }
}
