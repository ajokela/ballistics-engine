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

/// One quantized dial setting: the detent count nearest `angle`, and what remains.
/// `residual` is in the same angular unit as `angle` (the click's base unit) and is
/// DEFINED by the exact reconstruction identity
/// `angle == clicks as f64 * click.size + residual` (bit-for-bit, by construction).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quantized {
    pub clicks: i64,
    pub residual: f64,
}

/// Quantize an angular value (in `click.base` units) onto a click detent.
/// Rounding is nearest, ties away from zero — identical to `clicks_for`, which delegates here.
pub fn quantize_angle(angle: f64, click: &ClickValue) -> Quantized {
    let clicks_f = (angle / click.size).round();
    Quantized { clicks: clicks_f as i64, residual: angle - clicks_f * click.size }
}

/// This click's size in milliradians, the crate's true-angular unit. `adjustment_factor`
/// maps one physical ratio into each display unit, so conversion is a factor ratio —
/// using the LOCKED printed-table MOA constant (3438; MBA-724), deliberately not 3437.7467.
pub fn click_size_mil(click: &ClickValue) -> f64 {
    click.size * adjustment_factor(ClickBase::Mil) / adjustment_factor(click.base)
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
    quantize_angle(angle, click).clicks
}

/// Zero-banner dial values `(MOA, mrad)` for a solved zero angle, corrected by the
/// elevation tracking correction factor (MBA-1358).
///
/// Direction: the stored CF is the published tall-target ratio `actual measured travel
/// / dialed travel` (0.95 = the scope under-tracks by 5%). To obtain a true angular
/// need `N` on a scope whose dialed unit delivers `CF` true units, the number set on
/// the dial must be `N / CF` — so dial-unit OUTPUTS are DIVIDED by the CF.
///
/// Replicates the WASM terminal's historical banner conversions bit-for-bit at
/// `elevation_cf == 1.0`: `deg * 60.0` (true MOA — deliberately NOT this module's locked
/// printed-table `Moa = 3438.0` factor; see the MBA-724 note in the wasm banner sites)
/// and `rad * 1000.0`, each divided by exactly `1.0` (a bit-exact no-op). Lives here,
/// in an ungated module, so the emit rule is host-testable (`wasm.rs` is wasm32-gated)
/// and shared by all three banner sites instead of being edited in parallel.
pub fn zero_banner_dial_values(angle_deg: f64, angle_rad: f64, elevation_cf: f64) -> (f64, f64) {
    (
        angle_deg * 60.0 / elevation_cf,
        angle_rad * 1000.0 / elevation_cf,
    )
}

/// MBA-1358: the ONE accepted band for a scope tracking correction factor, shared by the
/// native CLI and the WASM terminal so the two surfaces cannot drift: strictly between
/// 0.5 and 1.5 (a real scope does not mistrack by 50% — values outside the band mean the
/// tall-target test or a hand edit went wrong) and finite.
pub fn tracking_cf_in_range(value: f64) -> bool {
    value.is_finite() && value > 0.5 && value < 1.5
}

#[cfg(test)]
mod tests {
    use super::*;

    // MBA-1358: the WASM zero-banner emit rule, host-tested here because wasm.rs is
    // wasm32-gated (the drag_coefficient_json_value pattern).
    #[test]
    fn zero_banner_dial_values_divide_by_the_elevation_cf_and_are_exact_at_one() {
        let deg = 0.0717_f64;
        let rad = deg.to_radians();
        // cf = 1.0 is bit-exact against the historical inline expressions.
        let (moa, mrad) = zero_banner_dial_values(deg, rad, 1.0);
        assert_eq!(moa.to_bits(), (deg * 60.0).to_bits());
        assert_eq!(mrad.to_bits(), (rad * 1000.0).to_bits());
        // an under-tracking scope (CF < 1) needs MORE dial: outputs divide by the CF.
        let (moa_cf, mrad_cf) = zero_banner_dial_values(deg, rad, 0.95);
        assert_eq!(moa_cf.to_bits(), (deg * 60.0 / 0.95).to_bits());
        assert_eq!(mrad_cf.to_bits(), (rad * 1000.0 / 0.95).to_bits());
        assert!(moa_cf > moa && mrad_cf > mrad);
    }

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

    #[test]
    fn quantize_reconstruction_identity_is_bit_exact() {
        // residual is DEFINED as angle - clicks*size, so reconstruction is exact by construction.
        let c = ClickValue { size: 0.25, base: ClickBase::Moa };
        for angle in [0.0, 0.24, 0.25, 0.26, -std::f64::consts::PI, 7.4499999, 1234.567] {
            let q = quantize_angle(angle, &c);
            assert_eq!(q.clicks as f64 * c.size + q.residual, angle, "identity broke at {angle}");
            assert!(q.residual.abs() <= c.size / 2.0 + f64::EPSILON, "residual beyond half-click");
        }
    }

    #[test]
    fn quantize_exact_multiples_have_residual_exactly_zero() {
        let c = ClickValue { size: 0.1, base: ClickBase::Mil };
        let q = quantize_angle(0.1 * 27.0, &c);
        assert_eq!(q.clicks, 27);
        assert_eq!(q.residual, 0.0); // bit-exact, not approx
    }

    #[test]
    fn quantize_ties_round_away_from_zero_matching_clicks_for() {
        let c = ClickValue { size: 0.5, base: ClickBase::Mil };
        assert_eq!(quantize_angle(0.25, &c).clicks, 1);   // f64::round: ties away from zero
        assert_eq!(quantize_angle(-0.25, &c).clicks, -1);
    }

    #[test]
    fn clicks_for_is_unchanged_and_delegates() {
        // Historical behavior pinned: sub-1yd guard + rounding. Values chosen off any tie.
        let c = ClickValue { size: 0.25, base: ClickBase::Moa };
        assert_eq!(clicks_for(10.0, 0.5, &c), 0);                       // range_yd < 1.0 guard
        let angle = (7.2 / 300.0) * adjustment_factor(ClickBase::Moa);
        assert_eq!(clicks_for(7.2, 300.0, &c), (angle / 0.25_f64).round() as i64);
    }

    #[test]
    fn click_size_mil_converts_via_the_locked_factor_table() {
        let moa = ClickValue { size: 0.25, base: ClickBase::Moa };
        // 0.25 MOA in mil via the LOCKED 3438 constant (MBA-724), not 3437.7467.
        assert_eq!(click_size_mil(&moa), 0.25 * 1000.0 / 3438.0);
        let mil = ClickValue { size: 0.1, base: ClickBase::Mil };
        assert_eq!(click_size_mil(&mil), 0.1);
    }
}
