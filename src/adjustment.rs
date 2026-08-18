//! Turret adjustment-unit conversions and click-value parsing (MBA-1355).
//!
//! Shared by the CLI (`main.rs`) and the WASM terminal (`wasm.rs`). This lives in the
//! library crate — not `main.rs` — so it compiles for `wasm32-unknown-unknown` with no
//! feature gate. To avoid a circular dependency on `main.rs`'s clap `AdjustmentUnit`
//! `ValueEnum`, this module defines its own minimal angular base (`ClickBase`) and factor
//! table; `main.rs` maps its own `AdjustmentUnit` onto `ClickBase`/`adjustment_factor`
//! locally (see `drop_to_adjustment` in `main.rs`).

use serde::{Deserialize, Deserializer, Serialize, Serializer};

/// Angular base unit for a turret click graduation (MBA-1355). A click graduation is
/// always expressed in mil, (true) MOA, or SMOA/IPHY — never in whole clicks itself.
///
/// Does not itself derive `Serialize`/`Deserialize`: `ClickValue` serializes as ONE
/// suffixed string (`"0.1mil"`, `"0.25moa"`, ...) via its own hand-written impls below,
/// not as a structural `{size, base}` object, so nothing needs this type's wire form on
/// its own (MBA-1348).
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

/// Serializes as the crate's one existing suffixed-string representation — `"0.1mil"`,
/// `"0.25moa"`, `"1smoa"` — the exact shape `parse_click_value` parses and a CLI
/// `--elevation-click` flag already accepts (MBA-1348). Deliberately NOT the structural
/// `{"size": 0.1, "base": "mil"}` a plain derive would produce: the string form is stable
/// against internal field changes and makes a profile file's click entries identical to
/// the CLI flag syntax a shooter already types. `size` uses `f64`'s default `Display`,
/// which is Rust's shortest string that round-trips back to the same value (so `1.0`
/// prints as `"1"`, matching `parse_click_value`'s own accepted syntax).
impl Serialize for ClickValue {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let suffix = match self.base {
            ClickBase::Mil => "mil",
            ClickBase::Moa => "moa",
            ClickBase::Smoa => "smoa",
        };
        let size = self.size;
        serializer.serialize_str(&format!("{size}{suffix}"))
    }
}

/// Deserializes via `parse_click_value` — the ONLY parser for this syntax, not a second,
/// divergent implementation of it — so a malformed string is rejected with the identical
/// message the CLI already gives, and `parse_click_value`'s positive-and-finite rule
/// applies on every deserialize (MBA-1348).
impl<'de> Deserialize<'de> for ClickValue {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let raw = String::deserialize(deserializer)?;
        parse_click_value(&raw).map_err(serde::de::Error::custom)
    }
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

    // MBA-1348: ClickValue's JSON wire form is the same suffixed string
    // `parse_click_value` already parses, not a structural `{size, base}` object.
    #[test]
    fn click_value_json_is_pinned_to_the_suffixed_string() {
        let quarter_moa = ClickValue { size: 0.25, base: ClickBase::Moa };
        assert_eq!(serde_json::to_string(&quarter_moa).unwrap(), "\"0.25moa\"");

        let tenth_mil = ClickValue { size: 0.1, base: ClickBase::Mil };
        assert_eq!(serde_json::to_string(&tenth_mil).unwrap(), "\"0.1mil\"");

        // A whole-number size prints without a trailing ".0" -- `f64`'s shortest Display.
        let one_smoa = ClickValue { size: 1.0, base: ClickBase::Smoa };
        assert_eq!(serde_json::to_string(&one_smoa).unwrap(), "\"1smoa\"");
    }

    #[test]
    fn click_value_round_trips_through_json_including_smoa() {
        let cases = [
            ClickValue { size: 0.1, base: ClickBase::Mil },
            ClickValue { size: 0.25, base: ClickBase::Moa },
            ClickValue { size: 1.0, base: ClickBase::Smoa },
            ClickValue { size: 0.125, base: ClickBase::Smoa },
        ];
        for click in cases {
            let json = serde_json::to_string(&click).unwrap();
            let parsed: ClickValue = serde_json::from_str(&json).unwrap();
            assert_eq!(parsed, click, "round trip broke for {json}");
        }
    }

    #[test]
    fn click_value_deserialize_rejects_what_parse_click_value_rejects() {
        // The Deserialize impl IS parse_click_value -- a non-positive or malformed
        // string must fail to deserialize exactly as it fails to parse.
        for bad in ["\"0mil\"", "\"-0.1mil\"", "\"0.25\"", "\"nonsense\""] {
            assert!(
                serde_json::from_str::<ClickValue>(bad).is_err(),
                "{bad} must be rejected"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Adjustment display layer (moved from main.rs for the bridge card service).
//
// Everything below was previously private to the CLI. The bridge's card
// commands need the exact same unit conversion, zero-set bias, tracking-CF,
// and click-quantization ordering the CLI prints, so the shared boundary now
// lives here and `main.rs` re-imports it. The module doc's historical note
// about avoiding a circular dependency on main.rs's AdjustmentUnit no longer
// applies: the enum itself moved here, and the CLI derives its clap face via
// cfg_attr.
// ---------------------------------------------------------------------------

/// Printed adjustment unit for dial columns (MIL/MOA/SMOA/IPHY or whole clicks).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
pub enum AdjustmentUnit {
    /// Milliradians (1 MIL = 3.6 inches at 100 yards)
    #[default]
    Mil,
    /// Minutes of Angle, true MOA (1 MOA = 1.047 inches at 100 yards)
    Moa,
    /// Shooter's MOA (exactly 1 inch per 100 yards)
    Smoa,
    /// Inches per hundred yards (numerically identical to SMOA)
    Iphy,
    /// Whole turret clicks; requires an elevation click graduation
    /// (--elevation-click-value or the profile's elevation_click)
    Clicks,
}

/// Convert drop to adjustment unit (MIL, MOA, or SMOA/IPHY). `unit` maps onto
/// `ClickBase`/`adjustment_factor` (MBA-1355) so the CLI, the WASM terminal, and the
/// bridge card service share one conversion table.
///
/// `AdjustmentUnit::Clicks` has no fixed factor of its own — clicks depend on a graduation
/// (`ClickValue`) supplied separately — so callers MUST resolve clicks via `clicks_for()`
/// before ever reaching this function with `unit == Clicks`. That arm exists only as a
/// release-safe fallback (falls back to MIL) guarded by a `debug_assert!`.
pub fn drop_to_adjustment(drop_yd: f64, range_yd: f64, unit: AdjustmentUnit) -> f64 {
    if range_yd < 1.0 {
        return 0.0;
    }
    let factor = match unit {
        AdjustmentUnit::Mil => adjustment_factor(ClickBase::Mil),
        AdjustmentUnit::Moa => adjustment_factor(ClickBase::Moa),
        AdjustmentUnit::Smoa | AdjustmentUnit::Iphy => adjustment_factor(ClickBase::Smoa),
        AdjustmentUnit::Clicks => {
            // Callers resolve clicks via clicks_for() BEFORE display; reaching this
            // arm is a scoping bug. Fall back to MIL so release builds stay sane.
            debug_assert!(false, "Clicks must be resolved via clicks_for()");
            adjustment_factor(ClickBase::Mil)
        }
    };
    (drop_yd / range_yd) * factor
}

/// The `adjustment_factor` a non-clicks [`AdjustmentUnit`] renders in, for rescaling a
/// MIL-denominated zero-set bias into the active display unit (MBA-1360). `Clicks`
/// never reaches this (the clicks arm below works on raw drop); the MIL fallback
/// mirrors `drop_to_adjustment`'s own release-safe arm.
pub fn unit_factor_for_bias(unit: AdjustmentUnit) -> f64 {
    match unit {
        AdjustmentUnit::Mil => adjustment_factor(ClickBase::Mil),
        AdjustmentUnit::Moa => adjustment_factor(ClickBase::Moa),
        AdjustmentUnit::Smoa | AdjustmentUnit::Iphy => adjustment_factor(ClickBase::Smoa),
        AdjustmentUnit::Clicks => {
            debug_assert!(false, "Clicks bias is applied on raw drop, not here");
            adjustment_factor(ClickBase::Mil)
        }
    }
}

/// An adjustment value plus the click-quantization detail (Plan B Task 2). `value` is
/// the whole-click count as an `f64` on the clicks arm, the plain angle otherwise.
/// `quantized` carries the sub-click `residual` so a caller that needs it doesn't have
/// to re-derive the click math a second time; `None` on the angular arm, where nothing
/// is quantized.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdjustmentDisplay {
    pub value: f64,
    pub quantized: Option<Quantized>,
}

/// The one boundary where zero-set bias, tracking CF, and click quantization compose —
/// the order is load-bearing (see main.rs's `zero_set_bias_applies_before_the_cf_division`
/// and `clicks_arm_order_of_operations_is_load_bearing` tests) and must not move.
pub fn adjustment_display(
    drop_yd: f64,
    range_yd: f64,
    unit: AdjustmentUnit,
    click: Option<ClickValue>,
    bias_mil: f64,
    cf: f64,
) -> AdjustmentDisplay {
    match click {
        // Clicks convert from RAW drop (not from the angular value), so the CF is
        // applied to the input drop — correcting the angle before click quantization.
        // The zero-set bias joins as its drop-equivalent at this range (1 mil =
        // range/1000), keeping the ordering: (true + bias) first, / CF second, then
        // quantize into whole clicks.
        Some(c) => {
            let biased_drop_yd = if bias_mil != 0.0 {
                drop_yd + bias_mil / 1000.0 * range_yd
            } else {
                drop_yd
            };
            let drop_for_clicks = biased_drop_yd / cf;
            // Inlined from `clicks_for` (guard, then angle formula) so this boundary can
            // quantize once and keep the residual, instead of calling `clicks_for` and
            // throwing the sub-click remainder away.
            if range_yd < 1.0 {
                AdjustmentDisplay {
                    value: 0.0,
                    quantized: Some(Quantized { clicks: 0, residual: 0.0 }),
                }
            } else {
                let angle = (drop_for_clicks / range_yd) * adjustment_factor(c.base);
                let q = quantize_angle(angle, &c);
                AdjustmentDisplay { value: q.clicks as f64, quantized: Some(q) }
            }
        }
        None => {
            let true_need = drop_to_adjustment(drop_yd, range_yd, unit);
            let biased = if bias_mil != 0.0 {
                // Rescale the mil bias into the display unit via the shared factor
                // table (MIL -> unit), so MOA/SMOA/IPHY columns get the same angular
                // correction the MIL column does.
                true_need
                    + bias_mil * unit_factor_for_bias(unit) / adjustment_factor(ClickBase::Mil)
            } else {
                true_need
            };
            AdjustmentDisplay { value: biased / cf, quantized: None }
        }
    }
}

/// Windage-axis adjustment display (Tier 2 review C1 fix): only pass the click
/// graduation through to `adjustment_display` when the windage axis itself is `Clicks`;
/// otherwise this collapses to the angular path. All windage columns (range-table,
/// wind-card, compare, the PDF dope card's wind_adj/lead_adj, and the bridge card
/// service) go through this one function so the guard cannot be dropped again.
pub fn windage_adjustment_display(
    drift_yd: f64,
    range_yd: f64,
    windage_unit: AdjustmentUnit,
    windage_click: Option<ClickValue>,
    windage_bias_mil: f64,
    windage_cf: f64,
) -> AdjustmentDisplay {
    let click = if windage_unit == AdjustmentUnit::Clicks {
        windage_click
    } else {
        None
    };
    adjustment_display(drift_yd, range_yd, windage_unit, click, windage_bias_mil, windage_cf)
}

/// Display label for an [`AdjustmentUnit`] header/column name (MBA-1355/MBA-1410):
/// every surface that prints one shares this, so the "MIL"/"MOA"/"SMOA"/"IPHY"/"CLICKS"
/// spelling cannot drift between commands.
pub fn adjustment_unit_label(unit: AdjustmentUnit) -> String {
    match unit {
        AdjustmentUnit::Mil => "MIL".to_string(),
        AdjustmentUnit::Moa => "MOA".to_string(),
        AdjustmentUnit::Smoa => "SMOA".to_string(),
        AdjustmentUnit::Iphy => "IPHY".to_string(),
        AdjustmentUnit::Clicks => "CLICKS".to_string(),
    }
}
