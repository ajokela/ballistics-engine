//! Power-factor calculator and per-organization rulebook pass/fail (MBA-1372).
//!
//! Power factor (PF) is `bullet weight (grains) * velocity (fps) / 1000` -- a shared,
//! organization-agnostic formula used by essentially every American practical/action
//! shooting sport to gate ammunition into scoring categories. This module computes from
//! bullet weight and velocity only; it does not run a trajectory.
//!
//! ## Truncation convention -- a real rules detail, not rounding noise
//!
//! USPSA's rulebook is explicit: "The final result will ignore all decimal places (e.g.
//! for USPSA purposes, a result of 124.9999 is not 125)" (*USPSA Competition Rules*,
//! March 2026 edition, Sec. 5.6, Rule 36). IDPA's rulebook states the same convention in
//! different words: "...divide by 1000, and ignore numbers to the right of the decimal"
//! (*2024 IDPA Rulebook*, Sec. 8.3.4.7). Both **truncate toward zero (floor)** -- neither
//! rounds.
//!
//! [`scored_power_factor`] reports a single canonical PF as `floor(weight_gr * v_fps /
//! 1000)`, matching those two rulebooks' own scored value exactly. SASS's handbook does
//! not spell out a truncation rule for its own threshold check (its worked examples
//! happen to land on exact values), so applying the same floor to the SASS row too is a
//! deliberate, conservative choice: it never reports a marginal load as passing a SASS
//! minimum that an official chronograph session could truncate below.
//!
//! ## Boundary semantics: `>=` ("meets or exceeds"), never strict `>`
//!
//! Every source rulebook uses "meets or exceeds" / "at least" / "not less than" language,
//! never a strict "greater than":
//!   - USPSA: rule tables read "Minimum power factor for Major: 165" / "...for Minor:
//!     125"; Rule 37: "...until the competitor **meets** the minimum power factor...".
//!   - IDPA, Sec. 8.3.5.1.1: "If two of the three rounds **meet or exceed** the required
//!     power factor, the ammunition is in compliance."
//!   - SASS *Shooter's Handbook*: "not less than a minimum power factor of 60" / "If the
//!     average velocity...**meets or exceeds** the calculated power factor of 60".
//!
//! A load landing exactly ON a threshold (PF == 125, velocity == 400, etc.) therefore
//! **passes**. Every threshold in [`PF_THRESHOLDS`] is evaluated as `pf >= minimum` and
//! every velocity cap as `v <= maximum` / `v >= minimum` (inclusive on both ends) --
//! never a strict inequality. [`evaluate_all`]'s boundary tests exercise both sides of
//! every published number in the table.
//!
//! ## Data sources (published rulebook facts -- no copyright exposure; cited with revision)
//!
//! - USPSA: *USPSA Competition Rules*, March 2026 edition, Sec. 5.6 (per-division rule
//!   tables: "Minimum power factor for Major: 165" / "...for Minor: 125").
//! - IDPA: *2024 IDPA Rulebook*, Sec. 8.3.4 ("Ammunition minimum power factors"):
//!   8.3.4.1 SSP/ESP/CO 125; 8.3.4.2 CDP 165; 8.3.4.3 Stock Revolver/CCP 105;
//!   8.3.4.4 Enhanced Revolver 155; 8.3.4.5 BUG 95; 8.3.4.6 PCC 135.
//! - SASS: *SASS Cowboy Action Shooting Shooter's Handbook*, Version 28, January 2026,
//!   "Power Factors" and Glossary entry "Power factor": minimum PF 60, minimum velocity
//!   400 fps, maximum velocity 1000 fps (revolver/pistol) / 1400 fps (rifle).

/// One organization/class row in the power-factor rulebook table.
#[derive(Debug, Clone, Copy)]
pub struct PfThreshold {
    pub organization: &'static str,
    pub class: &'static str,
    pub min_pf: f64,
    /// Minimum permitted velocity (fps), if the org's rule imposes one (SASS only).
    pub min_velocity_fps: Option<f64>,
    /// Maximum permitted velocity (fps), if the org's rule imposes one (SASS only).
    pub max_velocity_fps: Option<f64>,
    /// Rulebook citation, including edition/revision (kept alongside the numbers so a
    /// future rulebook update is a one-line diff against a named source).
    pub source: &'static str,
}

/// Rulebook power-factor/velocity thresholds. A DATA TABLE by design (per the brief) so
/// updating a rulebook figure is a one-line edit here, not a code change elsewhere. See
/// the module docs for citations.
pub const PF_THRESHOLDS: &[PfThreshold] = &[
    PfThreshold {
        organization: "USPSA",
        class: "Minor",
        min_pf: 125.0,
        min_velocity_fps: None,
        max_velocity_fps: None,
        source: "USPSA Competition Rules, March 2026, Sec. 5.6",
    },
    PfThreshold {
        organization: "USPSA",
        class: "Major",
        min_pf: 165.0,
        min_velocity_fps: None,
        max_velocity_fps: None,
        source: "USPSA Competition Rules, March 2026, Sec. 5.6",
    },
    PfThreshold {
        organization: "IDPA",
        class: "BUG",
        min_pf: 95.0,
        min_velocity_fps: None,
        max_velocity_fps: None,
        source: "2024 IDPA Rulebook, Sec. 8.3.4.5",
    },
    PfThreshold {
        organization: "IDPA",
        class: "Stock Revolver / CCP",
        min_pf: 105.0,
        min_velocity_fps: None,
        max_velocity_fps: None,
        source: "2024 IDPA Rulebook, Sec. 8.3.4.3",
    },
    PfThreshold {
        organization: "IDPA",
        class: "SSP / ESP / CO",
        min_pf: 125.0,
        min_velocity_fps: None,
        max_velocity_fps: None,
        source: "2024 IDPA Rulebook, Sec. 8.3.4.1",
    },
    PfThreshold {
        organization: "IDPA",
        class: "PCC",
        min_pf: 135.0,
        min_velocity_fps: None,
        max_velocity_fps: None,
        source: "2024 IDPA Rulebook, Sec. 8.3.4.6",
    },
    PfThreshold {
        organization: "IDPA",
        class: "Enhanced Revolver",
        min_pf: 155.0,
        min_velocity_fps: None,
        max_velocity_fps: None,
        source: "2024 IDPA Rulebook, Sec. 8.3.4.4",
    },
    PfThreshold {
        organization: "IDPA",
        class: "CDP",
        min_pf: 165.0,
        min_velocity_fps: None,
        max_velocity_fps: None,
        source: "2024 IDPA Rulebook, Sec. 8.3.4.2",
    },
    PfThreshold {
        organization: "SASS",
        class: "Smokeless - Pistol-Revolver",
        min_pf: 60.0,
        min_velocity_fps: Some(400.0),
        max_velocity_fps: Some(1000.0),
        source: "SASS Cowboy Action Shooting Shooter's Handbook, Version 28, Jan 2026",
    },
    PfThreshold {
        organization: "SASS",
        class: "Smokeless - Rifle",
        min_pf: 60.0,
        min_velocity_fps: Some(400.0),
        max_velocity_fps: Some(1400.0),
        source: "SASS Cowboy Action Shooting Shooter's Handbook, Version 28, Jan 2026",
    },
];

/// Raw power factor: `weight_gr * velocity_fps / 1000`, no rounding/truncation.
pub fn power_factor(weight_gr: f64, velocity_fps: f64) -> f64 {
    weight_gr * velocity_fps / 1000.0
}

/// The truncated (floored) PF used for rulebook scoring/threshold comparisons. See the
/// module docs' "Truncation convention" section for why floor (not round) is correct,
/// and why it is applied uniformly across all three organizations.
pub fn scored_power_factor(weight_gr: f64, velocity_fps: f64) -> f64 {
    power_factor(weight_gr, velocity_fps).floor()
}

/// One rulebook row's evaluation against a computed load.
#[derive(Debug, Clone, Copy)]
pub struct PfEvaluation {
    pub organization: &'static str,
    pub class: &'static str,
    pub min_pf: f64,
    pub pf_pass: bool,
    pub min_velocity_fps: Option<f64>,
    pub max_velocity_fps: Option<f64>,
    /// `None` when the row has no velocity constraint at all (only SASS rows do).
    pub velocity_pass: Option<bool>,
    /// Overall pass: `pf_pass && velocity_pass.unwrap_or(true)`.
    pub pass: bool,
}

/// Evaluate a load (bullet weight in grains, velocity in fps) against every row in
/// [`PF_THRESHOLDS`]. Boundary semantics: `>=`/`<=` throughout (see module docs).
pub fn evaluate_all(weight_gr: f64, velocity_fps: f64) -> Vec<PfEvaluation> {
    let pf = scored_power_factor(weight_gr, velocity_fps);
    PF_THRESHOLDS
        .iter()
        .map(|t| {
            let pf_pass = pf >= t.min_pf;
            let min_ok = t.min_velocity_fps.is_none_or(|m| velocity_fps >= m);
            let max_ok = t.max_velocity_fps.is_none_or(|m| velocity_fps <= m);
            let velocity_pass = if t.min_velocity_fps.is_some() || t.max_velocity_fps.is_some() {
                Some(min_ok && max_ok)
            } else {
                None
            };
            PfEvaluation {
                organization: t.organization,
                class: t.class,
                min_pf: t.min_pf,
                pf_pass,
                min_velocity_fps: t.min_velocity_fps,
                max_velocity_fps: t.max_velocity_fps,
                velocity_pass,
                pass: pf_pass && velocity_pass.unwrap_or(true),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// SASS's own worked examples (Shooter's Handbook, "Power Factors" section).
    #[test]
    fn sass_worked_examples() {
        assert!((power_factor(100.0, 600.0) - 60.0).abs() < 1e-9);
        assert!((power_factor(77.0, 800.0) - 61.6).abs() < 1e-9);
        assert!((power_factor(200.0, 400.0) - 80.0).abs() < 1e-9);
    }

    /// USPSA's own worked truncation example (Sec. 5.6, Rule 36): a raw PF of 124.9999
    /// truncates to 124, i.e. it is NOT 125.
    #[test]
    fn uspsa_truncation_example() {
        // weight*velocity/1000 = 124.9999 exactly.
        let raw = 124.9999_f64;
        let weight_gr = 1000.0;
        let velocity_fps = raw; // weight=1000 makes velocity numerically equal to the raw PF
        assert!((power_factor(weight_gr, velocity_fps) - 124.9999).abs() < 1e-9);
        assert_eq!(scored_power_factor(weight_gr, velocity_fps), 124.0);
    }

    /// IDPA's own worked truncation example (Sec. 8.3.4.7): 230.1gr @ 794.7fps ->
    /// 182.86047, reported as 182 power factor (decimal digits ignored/truncated).
    #[test]
    fn idpa_truncation_example() {
        let raw = power_factor(230.1, 794.7);
        assert!((raw - 182.86047).abs() < 1e-3);
        assert_eq!(scored_power_factor(230.1, 794.7), 182.0);
    }

    /// Every published threshold: exactly AT the boundary must PASS (>=, not >), and one
    /// PF point below must FAIL. Exercises both sides of every number in the table.
    #[test]
    fn every_pf_threshold_boundary_both_sides() {
        for t in PF_THRESHOLDS {
            // A bullet weight of 1000gr makes velocity_fps numerically equal to the raw
            // (and, since it's already a whole number, scored) power factor.
            let at_threshold = evaluate_all(1000.0, t.min_pf);
            let row_at = at_threshold
                .iter()
                .find(|e| e.organization == t.organization && e.class == t.class)
                .expect("row present");
            assert!(
                row_at.pf_pass,
                "{} {} at exactly its min_pf {} must PASS (>=, not >)",
                t.organization, t.class, t.min_pf
            );

            let below_threshold = evaluate_all(1000.0, t.min_pf - 1.0);
            let row_below = below_threshold
                .iter()
                .find(|e| e.organization == t.organization && e.class == t.class)
                .expect("row present");
            assert!(
                !row_below.pf_pass,
                "{} {} one PF point below min_pf {} must FAIL",
                t.organization, t.class, t.min_pf
            );
        }
    }

    /// SASS velocity caps: exactly at the min/max boundary must PASS; one fps outside
    /// must FAIL. A PF of 60 (100gr @ 600fps => 60.0 exactly) keeps the PF check passing
    /// throughout so only the velocity boundary is under test.
    #[test]
    fn sass_velocity_boundaries_both_sides() {
        // Minimum velocity (400 fps): exactly at it passes; one fps under fails.
        let weight_gr = 200.0; // PF = weight*v/1000; chosen so PF stays >= 60 even at 399 fps.
        let at_min = evaluate_all(weight_gr, 400.0);
        let below_min = evaluate_all(weight_gr, 399.0);
        for class in ["Smokeless - Pistol-Revolver", "Smokeless - Rifle"] {
            let row = at_min
                .iter()
                .find(|e| e.organization == "SASS" && e.class == class)
                .unwrap();
            assert!(row.velocity_pass == Some(true), "{class} at 400 fps must pass");
            assert!(row.pass, "{class} at 400 fps must pass overall");

            let row = below_min
                .iter()
                .find(|e| e.organization == "SASS" && e.class == class)
                .unwrap();
            assert_eq!(row.velocity_pass, Some(false), "{class} at 399 fps must fail");
            assert!(!row.pass);
        }

        // Maximum velocity: pistol/revolver cap 1000 fps, rifle cap 1400 fps.
        let pistol_weight = 65.0; // 65gr @ 1000fps -> PF 65.0 (>=60); @1001fps -> 65.065 (still >=60)
        let at_pistol_max = evaluate_all(pistol_weight, 1000.0);
        let over_pistol_max = evaluate_all(pistol_weight, 1001.0);
        let row = at_pistol_max
            .iter()
            .find(|e| e.organization == "SASS" && e.class == "Smokeless - Pistol-Revolver")
            .unwrap();
        assert_eq!(row.velocity_pass, Some(true));
        let row = over_pistol_max
            .iter()
            .find(|e| e.organization == "SASS" && e.class == "Smokeless - Pistol-Revolver")
            .unwrap();
        assert_eq!(row.velocity_pass, Some(false));

        let rifle_weight = 65.0;
        let at_rifle_max = evaluate_all(rifle_weight, 1400.0);
        let over_rifle_max = evaluate_all(rifle_weight, 1401.0);
        let row = at_rifle_max
            .iter()
            .find(|e| e.organization == "SASS" && e.class == "Smokeless - Rifle")
            .unwrap();
        assert_eq!(row.velocity_pass, Some(true));
        let row = over_rifle_max
            .iter()
            .find(|e| e.organization == "SASS" && e.class == "Smokeless - Rifle")
            .unwrap();
        assert_eq!(row.velocity_pass, Some(false));

        // The pistol cap must not leak onto the rifle row and vice versa: 1200 fps is
        // over the pistol/revolver max (1000) but under the rifle max (1400).
        let mixed = evaluate_all(rifle_weight, 1200.0);
        let pistol_row = mixed
            .iter()
            .find(|e| e.organization == "SASS" && e.class == "Smokeless - Pistol-Revolver")
            .unwrap();
        let rifle_row = mixed
            .iter()
            .find(|e| e.organization == "SASS" && e.class == "Smokeless - Rifle")
            .unwrap();
        assert_eq!(pistol_row.velocity_pass, Some(false));
        assert_eq!(rifle_row.velocity_pass, Some(true));
    }

    /// Non-SASS rows carry no velocity constraint at all: `velocity_pass` is `None` and
    /// never gates the overall pass/fail.
    #[test]
    fn non_sass_rows_have_no_velocity_gate() {
        let rows = evaluate_all(1000.0, 200.0);
        for row in rows.iter().filter(|r| r.organization != "SASS") {
            assert_eq!(row.velocity_pass, None);
            assert_eq!(row.pass, row.pf_pass);
        }
    }

    /// A 9mm 147gr @ 900fps load: PF = 132.3 -> scored 132. Passes USPSA Minor (125),
    /// IDPA SSP/ESP/CO (125) and PCC (135)? No -- 132 < 135, so PCC fails while
    /// USPSA-Minor and SSP/ESP/CO pass. Cross-checks evaluate_all end to end.
    #[test]
    fn realistic_9mm_load_cross_check() {
        let rows = evaluate_all(147.0, 900.0);
        let uspsa_minor = rows
            .iter()
            .find(|e| e.organization == "USPSA" && e.class == "Minor")
            .unwrap();
        assert!(uspsa_minor.pass);
        let uspsa_major = rows
            .iter()
            .find(|e| e.organization == "USPSA" && e.class == "Major")
            .unwrap();
        assert!(!uspsa_major.pass);
        let idpa_ssp = rows
            .iter()
            .find(|e| e.organization == "IDPA" && e.class == "SSP / ESP / CO")
            .unwrap();
        assert!(idpa_ssp.pass);
        let idpa_pcc = rows
            .iter()
            .find(|e| e.organization == "IDPA" && e.class == "PCC")
            .unwrap();
        assert!(!idpa_pcc.pass);
    }

    #[test]
    fn table_has_expected_row_count_and_orgs() {
        assert_eq!(PF_THRESHOLDS.len(), 10);
        assert!(PF_THRESHOLDS.iter().any(|t| t.organization == "USPSA"));
        assert!(PF_THRESHOLDS.iter().any(|t| t.organization == "IDPA"));
        assert!(PF_THRESHOLDS.iter().any(|t| t.organization == "SASS"));
    }
}
