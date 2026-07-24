use crate::BCSegmentData;

/// Resolve a velocity-keyed BC table without assuming segment order.
///
/// Bands are half-open (`velocity_min <= v < velocity_max`), so a shared boundary belongs to
/// the upper band. Below/above global coverage, clamp to the lowest/highest-velocity band;
/// an interior coverage gap or empty table uses the caller's projectile-specific scalar BC.
///
/// MBA-1404: away from a transition this is byte-identical to the original step/clamp
/// lookup (see [`raw_velocity_segment_bc`], preserved verbatim). Near every discontinuity
/// class the raw lookup can produce -- an interior band-to-band boundary, either edge of an
/// interior band-to-fallback gap, or a global coverage entry/exit edge (below the lowest
/// band's `velocity_min`, above the highest band's `velocity_max`) -- this instead blends
/// across a small velocity margin centered on the boundary with a Hermite smoothstep
/// (`3t^2 - 2t^3`), so an integrator sampling BC every step does not see an instantaneous
/// drag-coefficient jolt exactly at a band edge. Tables with fewer than two segments have no
/// adjacent band to blend against and return the raw value unchanged.
pub(crate) fn velocity_segment_bc(
    velocity_fps: f64,
    segments: &[BCSegmentData],
    fallback_bc: f64,
) -> f64 {
    if segments.len() < 2 {
        // No adjacent band to blend against: single-band and empty tables stay exactly the
        // pre-MBA-1404 behavior.
        return raw_velocity_segment_bc(velocity_fps, segments, fallback_bc);
    }

    for (boundary, margin) in transition_boundaries(segments) {
        if margin <= 0.0 {
            // Degenerate (zero/negative-width) adjacent band: no blend, keep the hard step.
            continue;
        }
        let half = margin / 2.0;
        let lo_v = boundary - half;
        let hi_v = boundary + half;
        if velocity_fps < lo_v || velocity_fps > hi_v {
            continue;
        }

        let lo = raw_velocity_segment_bc(lo_v, segments, fallback_bc);
        let hi = raw_velocity_segment_bc(hi_v, segments, fallback_bc);
        if lo == hi {
            // Nothing actually jumps here (e.g. a continuous coverage-clamp edge): stay flat
            // rather than run smoothstep algebra that would return the same value anyway.
            return lo;
        }

        let t = ((velocity_fps - lo_v) / margin).clamp(0.0, 1.0);
        return lo + smoothstep(t) * (hi - lo);
    }

    raw_velocity_segment_bc(velocity_fps, segments, fallback_bc)
}

/// The pre-MBA-1404 step/clamp lookup, unchanged. Also used by [`velocity_segment_bc`] to
/// sample the value on each side of a boundary when building a blend.
fn raw_velocity_segment_bc(velocity_fps: f64, segments: &[BCSegmentData], fallback_bc: f64) -> f64 {
    if let Some(segment) = segments.iter().find(|segment| {
        velocity_fps >= segment.velocity_min && velocity_fps < segment.velocity_max
    }) {
        return segment.bc_value;
    }

    let lowest = segments
        .iter()
        .min_by(|a, b| a.velocity_min.total_cmp(&b.velocity_min));
    if let Some(segment) = lowest {
        if velocity_fps < segment.velocity_min {
            return segment.bc_value;
        }
    }

    let highest = segments
        .iter()
        .max_by(|a, b| a.velocity_max.total_cmp(&b.velocity_max));
    if let Some(segment) = highest {
        if velocity_fps >= segment.velocity_max {
            return segment.bc_value;
        }
    }

    fallback_bc
}

/// Every velocity at which [`raw_velocity_segment_bc`] can jump, paired with the margin to
/// blend across (`0.0` means "no blend": a degenerate/zero-width adjacent band). Order of
/// `segments` does not matter -- the boundaries are derived from a `velocity_min`-sorted
/// copy, so ascending- and descending-stored tables produce identical boundaries.
///
/// Discontinuity classes (see `velocity_segment_bc`'s doc comment): interior band-to-band
/// boundaries, the two edges of every interior band-to-fallback gap, and the two global
/// coverage entry/exit edges.
fn transition_boundaries(segments: &[BCSegmentData]) -> Vec<(f64, f64)> {
    let width = |s: &BCSegmentData| (s.velocity_max - s.velocity_min).max(0.0);

    let mut sorted: Vec<&BCSegmentData> = segments.iter().collect();
    sorted.sort_by(|a, b| a.velocity_min.total_cmp(&b.velocity_min));

    let mut boundaries = Vec::with_capacity(sorted.len() * 2 + 2);

    // Coverage entry: below the lowest-velocity_min band. The clamp region above has no
    // finite width, so the margin is driven entirely by the lowest band's own width.
    if let Some(first) = sorted.first() {
        boundaries.push((first.velocity_min, smoothing_margin(width(first))));
    }
    // Coverage exit: above the highest-velocity_max band. Mirrors raw()'s own "highest"
    // selection (max by `velocity_max`), which need not be `sorted.last()` for a
    // malformed/overlapping table.
    if let Some(highest) = segments
        .iter()
        .max_by(|a, b| a.velocity_max.total_cmp(&b.velocity_max))
    {
        boundaries.push((highest.velocity_max, smoothing_margin(width(highest))));
    }

    for pair in sorted.windows(2) {
        let (left, right) = (pair[0], pair[1]);
        if right.velocity_min > left.velocity_max {
            // Interior band-to-fallback gap: two edges, each capped by the gap's own width
            // in addition to its bordering band's width, so a narrow gap never lets the
            // blend eat into the band on the far side of it.
            let gap_width = right.velocity_min - left.velocity_max;
            boundaries.push((
                left.velocity_max,
                smoothing_margin(width(left).min(gap_width)),
            ));
            boundaries.push((
                right.velocity_min,
                smoothing_margin(gap_width.min(width(right))),
            ));
        } else if right.velocity_min == left.velocity_max {
            // Contiguous band-to-band boundary.
            boundaries.push((
                left.velocity_max,
                smoothing_margin(width(left).min(width(right))),
            ));
        }
        // Overlapping bands (right.velocity_min < left.velocity_max) are already an
        // ambiguous shape for raw()'s first-match lookup and aren't produced by any current
        // caller; MBA-1404 does not add blending for that out-of-scope configuration.
    }

    boundaries
}

/// `min(50.0 fps, 0.25 * narrower_adjacent_width)`; `0.0` (no blend) for a zero/negative
/// width, which is how a degenerate adjacent band disables blending at its boundaries.
fn smoothing_margin(narrower_adjacent_width: f64) -> f64 {
    if narrower_adjacent_width <= 0.0 {
        return 0.0;
    }
    (0.25 * narrower_adjacent_width).min(50.0)
}

/// Hermite smoothstep, `3t^2 - 2t^3`, clamped to `[0, 1]`.
fn smoothstep(t: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

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
                drop: 0.075, // 7.5% total drop for match boat tail
                transition_curve: 0.3,
            },
            BulletType::MatchFlatBase => BulletTypeFactors {
                drop: 0.10, // 10% for match flat base
                transition_curve: 0.35,
            },
            BulletType::HuntingBoatTail => BulletTypeFactors {
                drop: 0.15, // 15% for hunting boat tail
                transition_curve: 0.45,
            },
            BulletType::HuntingFlatBase => BulletTypeFactors {
                drop: 0.20, // 20% for hunting flat base
                transition_curve: 0.5,
            },
            BulletType::VldHighBc => BulletTypeFactors {
                drop: 0.05, // 5% for VLD (very low drag)
                transition_curve: 0.25,
            },
            BulletType::Hybrid => BulletTypeFactors {
                drop: 0.06, // 6% for hybrid designs
                transition_curve: 0.28,
            },
            BulletType::FMJ => BulletTypeFactors {
                drop: 0.12, // 12% for military ball
                transition_curve: 0.4,
            },
            BulletType::RoundNose => BulletTypeFactors {
                drop: 0.35, // 35% for round nose
                transition_curve: 0.7,
            },
            BulletType::Unknown => BulletTypeFactors {
                drop: 0.15, // Conservative 15%
                transition_curve: 0.5,
            },
        }
    }
}

/// BC segment estimator based on physics and known patterns
pub struct BCSegmentEstimator;

impl BCSegmentEstimator {
    /// Identify bullet type from model name and characteristics.
    ///
    /// This compatibility entry point interprets `bc_value` as a G1 BC. Call
    /// [`Self::identify_bullet_type_for_drag_model`] when the reference drag
    /// model is known.
    pub fn identify_bullet_type(
        model: &str,
        weight: f64,
        caliber: f64,
        bc_value: Option<f64>,
    ) -> BulletType {
        Self::identify_bullet_type_for_drag_model(model, weight, caliber, bc_value, "G1")
    }

    /// Identify bullet type while interpreting the BC in its reference-model space.
    pub fn identify_bullet_type_for_drag_model(
        model: &str,
        weight: f64,
        caliber: f64,
        bc_value: Option<f64>,
        drag_model: &str,
    ) -> BulletType {
        let model_lower = model.to_lowercase();

        // VLD/High BC bullets
        if model_lower.contains("vld")
            || model_lower.contains("berger")
            || model_lower.contains("hybrid")
            || model_lower.contains("elite")
        {
            if model_lower.contains("hybrid") {
                return BulletType::Hybrid;
            }
            return BulletType::VldHighBc;
        }

        // Match bullets (competition/target)
        if model_lower.contains("smk")
            || model_lower.contains("matchking")
            || model_lower.contains("match")
            || model_lower.contains("bthp")
            || model_lower.contains("competition")
            || model_lower.contains("target")
            || model_lower.contains("a-max")
            || model_lower.contains("eld-m")
            || model_lower.contains("scenar")
            || model_lower.contains("x-ring")
        {
            // Check for boat tail
            if model_lower.contains("bt") || model_lower.contains("boat") {
                return BulletType::MatchBoatTail;
            }
            // Check if high BC indicates boat tail (guard sd>0: calculate_sectional_density
            // returns 0 for non-positive caliber, which would make bc/sd == +Inf).
            if let Some(bc) = bc_value {
                let sd = Self::calculate_sectional_density(weight, caliber);
                if sd > 0.0 && Self::classification_bc_sd_ratio(bc, sd, drag_model) > 1.6 {
                    return BulletType::MatchBoatTail;
                }
            }
            return BulletType::MatchFlatBase;
        }

        // Hunting bullets (expanding)
        if model_lower.contains("gameking")
            || model_lower.contains("hunting")
            || model_lower.contains("sst")
            || model_lower.contains("eld-x")
            || model_lower.contains("partition")
            || model_lower.contains("accubond")
            || model_lower.contains("core-lokt")
            || model_lower.contains("ballistic tip")
            || model_lower.contains("v-max")
            || model_lower.contains("hornady sp")
            || model_lower.contains("interlock")
            || model_lower.contains("tsx")
        {
            // Check for boat tail
            if model_lower.contains("bt")
                || model_lower.contains("boat")
                || model_lower.contains("sst")
                || model_lower.contains("accubond")
            {
                return BulletType::HuntingBoatTail;
            }
            return BulletType::HuntingFlatBase;
        }

        // FMJ/Military
        if model_lower.contains("fmj")
            || model_lower.contains("ball")
            || model_lower.contains("m80")
            || model_lower.contains("m855")
            || model_lower.contains("tracer")
        {
            return BulletType::FMJ;
        }

        // Round nose
        if model_lower.contains("rn")
            || model_lower.contains("round nose")
            || model_lower.contains("rnsp")
        {
            return BulletType::RoundNose;
        }

        // Use BC value as hint if available. Guard sd>0 (zero for non-positive caliber)
        // so a degenerate input falls through to Unknown instead of dividing by zero
        // (bc/0 == +Inf, which would silently classify as VldHighBc).
        if let Some(bc) = bc_value {
            let sd = Self::calculate_sectional_density(weight, caliber);
            if sd > 0.0 {
                let bc_to_sd_ratio = Self::classification_bc_sd_ratio(bc, sd, drag_model);

                if bc_to_sd_ratio > 1.8 {
                    return BulletType::VldHighBc;
                } else if bc_to_sd_ratio > 1.5 {
                    return BulletType::MatchBoatTail;
                } else if bc_to_sd_ratio < 1.2 {
                    return BulletType::HuntingFlatBase;
                }
            }
        }

        BulletType::Unknown
    }

    /// Convert the reference-model-dependent BC/SD ratio into the G1 space used
    /// by the legacy coarse classification thresholds above. Typical boat-tail
    /// G1 BCs are approximately twice their G7 BCs; this normalization prevents
    /// ordinary G7 match bullets from looking like low-BC G1 flat-base bullets.
    fn classification_bc_sd_ratio(bc: f64, sd: f64, drag_model: &str) -> f64 {
        let g1_equivalent_bc = if drag_model.eq_ignore_ascii_case("G7") {
            bc * 2.0
        } else {
            bc
        };
        g1_equivalent_bc / sd
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
    #[allow(clippy::manual_clamp)] // max/min intentionally maps a NaN SD to the lower fallback.
    pub fn estimate_bc_segments(
        base_bc: f64,
        caliber: f64,
        weight: f64,
        model: &str,
        drag_model: &str,
    ) -> Vec<BCSegmentData> {
        // Identify bullet type
        let bullet_type = Self::identify_bullet_type_for_drag_model(
            model,
            weight,
            caliber,
            Some(base_bc),
            drag_model,
        );
        let type_factors = bullet_type.get_factors();

        // Calculate sectional density
        let sd = Self::calculate_sectional_density(weight, caliber);

        // Adjust BC drop based on sectional density
        // Higher SD = more stable BC
        let sd_factor = (sd / 0.25).max(0.7).min(1.3);
        let nominal_drop = type_factors.drop;

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
            }
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
            }
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
            }
            _ => {
                // Default degradation profile
                segments.push(BCSegmentData {
                    velocity_min: 2800.0,
                    velocity_max: 5000.0,
                    bc_value: base_bc,
                });

                let transonic_bc = base_bc * (1.0 - nominal_drop * 0.3);
                segments.push(BCSegmentData {
                    velocity_min: 1800.0,
                    velocity_max: 2800.0,
                    bc_value: transonic_bc,
                });

                let subsonic_bc = base_bc * (1.0 - nominal_drop);
                segments.push(BCSegmentData {
                    velocity_min: 0.0,
                    velocity_max: 1800.0,
                    bc_value: subsonic_bc,
                });
            }
        }

        // G7 reference drag follows modern boat-tail projectiles more closely,
        // so their banded BC varies less than the G1-shaped ladders above. Scale
        // each loss from nominal rather than the BC itself: this leaves the muzzle
        // band unchanged by this adjustment and makes the model effective for every
        // named and default profile. Do not run the identity algebra for G1, so
        // its established floating-point outputs remain bit-for-bit unchanged.
        if drag_model.eq_ignore_ascii_case("G7") {
            const G7_DROP_SCALE: f64 = 0.8;
            for segment in &mut segments {
                let drop_from_nominal = base_bc - segment.bc_value;
                segment.bc_value = base_bc - drop_from_nominal * G7_DROP_SCALE;
            }
        }

        // Sectional density shapes only the degradation depth. Scaling the BC
        // itself would alter the user's published muzzle value and would apply SD
        // twice in the default profile. Skip identity algebra so SD=0.25 keeps
        // established output bits unchanged.
        if sd_factor != 1.0 {
            for segment in &mut segments {
                let drop_from_nominal = base_bc - segment.bc_value;
                segment.bc_value = base_bc - drop_from_nominal / sd_factor;
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
        assert_eq!(
            BCSegmentEstimator::identify_bullet_type("168gr SMK", 168.0, 0.308, None),
            BulletType::MatchFlatBase
        );
        assert_eq!(
            BCSegmentEstimator::identify_bullet_type("168gr SMK BT", 168.0, 0.308, None),
            BulletType::MatchBoatTail
        );
        assert_eq!(
            BCSegmentEstimator::identify_bullet_type("150gr SST", 150.0, 0.308, None),
            BulletType::HuntingBoatTail
        );
        assert_eq!(
            BCSegmentEstimator::identify_bullet_type("147gr FMJ", 147.0, 0.308, None),
            BulletType::FMJ
        );
        assert_eq!(
            BCSegmentEstimator::identify_bullet_type("180gr RN", 180.0, 0.308, None),
            BulletType::RoundNose
        );
        assert_eq!(
            BCSegmentEstimator::identify_bullet_type("168gr VLD", 168.0, 0.308, None),
            BulletType::VldHighBc
        );
        assert_eq!(
            BCSegmentEstimator::identify_bullet_type("Some bullet", 150.0, 0.308, None),
            BulletType::Unknown
        );
    }

    #[test]
    fn test_sectional_density() {
        let sd = BCSegmentEstimator::calculate_sectional_density(168.0, 0.308);
        assert!((sd - 0.253).abs() < 0.001);
    }

    #[test]
    fn test_bc_estimation() {
        let segments =
            BCSegmentEstimator::estimate_bc_segments(0.450, 0.308, 168.0, "168gr SMK", "G1");

        // Match rifles typically have 4 segments
        assert!(segments.len() >= 3);
        // First segment should be close to base BC
        assert!((segments[0].bc_value - 0.450).abs() < 0.05);
        // BC should degrade at lower velocities
        assert!(segments[segments.len() - 1].bc_value < segments[0].bc_value);
    }

    #[test]
    fn g7_transition_adjustment_softens_each_band_drop() {
        let base_bc = 0.5;
        // SD = 0.25 exactly, so the independent sectional-density adjustment is neutral.
        let caliber = 1.0;
        let weight = 1750.0;

        for model in ["SMK BT", "FMJ"] {
            let g1 =
                BCSegmentEstimator::estimate_bc_segments(base_bc, caliber, weight, model, "G1");
            let g7 =
                BCSegmentEstimator::estimate_bc_segments(base_bc, caliber, weight, model, "G7");
            let lowercase_g7 =
                BCSegmentEstimator::estimate_bc_segments(base_bc, caliber, weight, model, "g7");
            assert_eq!(g7.len(), g1.len());
            assert_eq!(lowercase_g7.len(), g7.len());

            for ((g1_band, g7_band), lowercase_band) in
                g1.iter().zip(&g7).zip(&lowercase_g7)
            {
                assert_eq!(g7_band.velocity_min, g1_band.velocity_min);
                assert_eq!(g7_band.velocity_max, g1_band.velocity_max);
                assert_eq!(lowercase_band.velocity_min, g7_band.velocity_min);
                assert_eq!(lowercase_band.velocity_max, g7_band.velocity_max);
                assert_eq!(lowercase_band.bc_value.to_bits(), g7_band.bc_value.to_bits());
                let expected_g7 = base_bc - (base_bc - g1_band.bc_value) * 0.8;
                assert!(
                    (g7_band.bc_value - expected_g7).abs() < 1e-12,
                    "{model} band {}-{} did not soften the G1 loss: G1={}, G7={}, expected={expected_g7}",
                    g1_band.velocity_min,
                    g1_band.velocity_max,
                    g1_band.bc_value,
                    g7_band.bc_value
                );
            }
        }
    }

    #[test]
    fn sectional_density_shapes_only_band_degradation() {
        let base_bc = 0.5;
        let caliber = 0.224;
        let weight = 77.0;
        let sd = BCSegmentEstimator::calculate_sectional_density(weight, caliber);
        let sd_factor = (sd / 0.25).clamp(0.7, 1.3);

        for drag_model in ["G1", "G7"] {
            let model_drop_scale = if drag_model == "G7" { 0.8 } else { 1.0 };
            for (model, raw_retentions) in [
                ("SMK BT", &[1.0, 0.985, 0.965, 0.945, 0.925][..]),
                ("FMJ", &[1.0, 0.964, 0.88][..]),
            ] {
                let segments = BCSegmentEstimator::estimate_bc_segments(
                    base_bc,
                    caliber,
                    weight,
                    model,
                    drag_model,
                );
                assert_eq!(segments.len(), raw_retentions.len());

                for (segment, raw_retention) in segments.iter().zip(raw_retentions) {
                    let raw_drop = base_bc * (1.0 - raw_retention) * model_drop_scale;
                    let expected = base_bc - raw_drop / sd_factor;
                    assert!(
                        (segment.bc_value - expected).abs() < 1e-12,
                        "{drag_model} {model} band {}-{} misapplied SD: got {}, expected {expected}",
                        segment.velocity_min,
                        segment.velocity_max,
                        segment.bc_value
                    );
                }
                assert_eq!(
                    segments[0].bc_value.to_bits(),
                    base_bc.to_bits(),
                    "published muzzle BC must remain exact for {drag_model} {model}"
                );
            }
        }

        // High SD used to multiply every band above nominal and then cap them
        // all to base_bc, erasing the degradation ladder.
        let high_sd_base_bc = 0.3;
        let high_sd_segments = BCSegmentEstimator::estimate_bc_segments(
            high_sd_base_bc,
            0.308,
            220.0,
            "SMK BT",
            "G7",
        );
        assert_eq!(high_sd_segments[0].bc_value.to_bits(), high_sd_base_bc.to_bits());
        assert!(high_sd_segments.last().unwrap().bc_value < high_sd_base_bc);
        assert!((high_sd_segments.last().unwrap().bc_value - 0.28615384615384615).abs() < 1e-12);
    }

    #[test]
    fn generic_g7_bc_uses_g7_classification_space() {
        // A representative 175 gr .308 match bullet. Its G7 BC is ordinary for a
        // boat-tail projectile, but the same numeric value looks like a blunt
        // flat-base bullet when interpreted with the G1 BC/SD thresholds.
        let base_bc = 0.243;
        let segments = BCSegmentEstimator::estimate_bc_segments(base_bc, 0.308, 175.0, "", "G7");

        assert!(
            segments.len() >= 4,
            "G7 match bullet should use a near-flat match/VLD ladder: {segments:?}"
        );
        let subsonic_bc = segments.last().unwrap().bc_value;
        assert!(
            subsonic_bc >= base_bc * 0.92,
            "G7 match bullet was over-degraded: {base_bc} -> {subsonic_bc}"
        );

        // The normalization is deliberately equivalent to classifying the
        // approximate G1 BC, while the G7 transition adjustment then softens
        // the loss within that same ladder and the legacy entry point remains
        // G1-compatible.
        let g1_segments =
            BCSegmentEstimator::estimate_bc_segments(base_bc * 2.0, 0.308, 175.0, "", "G1");
        assert_eq!(segments.len(), g1_segments.len());
        let mut saw_g7_softening = false;
        for (g7, g1) in segments.iter().zip(&g1_segments) {
            assert_eq!(g7.velocity_min.to_bits(), g1.velocity_min.to_bits());
            assert_eq!(g7.velocity_max.to_bits(), g1.velocity_max.to_bits());
            let g7_retention = g7.bc_value / base_bc;
            let g1_retention = g1.bc_value / (base_bc * 2.0);
            assert!(g7_retention + 1e-12 >= g1_retention);
            saw_g7_softening |= g7_retention > g1_retention + 1e-12;
        }
        assert!(saw_g7_softening);

        let legacy_g1 = BCSegmentEstimator::identify_bullet_type("", 175.0, 0.308, Some(base_bc));
        assert_eq!(legacy_g1, BulletType::HuntingFlatBase);
        assert_eq!(
            legacy_g1,
            BCSegmentEstimator::identify_bullet_type_for_drag_model(
                "",
                175.0,
                0.308,
                Some(base_bc),
                "G1",
            )
        );
        assert_eq!(
            BCSegmentEstimator::identify_bullet_type_for_drag_model(
                "175gr SMK",
                175.0,
                0.308,
                Some(base_bc),
                "g7",
            ),
            BulletType::MatchBoatTail
        );
        assert_eq!(
            BCSegmentEstimator::identify_bullet_type_for_drag_model(
                "",
                175.0,
                0.0,
                Some(base_bc),
                "G7",
            ),
            BulletType::Unknown
        );
    }

    // MBA-1404: smoothstep continuity battery for `velocity_segment_bc`. Bands used here
    // pick bc_values that are exact binary fractions (0.25/0.5/0.75/etc.) wherever a test
    // hand-derives the blended midpoint, so assertions can use `assert_eq!` rather than an
    // epsilon tolerance.

    #[test]
    fn margin_caps_at_50_fps_for_wide_adjacent_bands() {
        // 0.25 * 1000 = 250, well above the 50 fps cap.
        assert_eq!(smoothing_margin(1000.0), 50.0);
        assert_eq!(smoothing_margin(200.0), 50.0); // 0.25*200 == 50, right at the cap edge
    }

    #[test]
    fn margin_uses_25_percent_rule_for_narrow_adjacent_bands() {
        // 0.25 * 40 = 10, under the 50 fps cap, so the 25% rule governs.
        assert_eq!(smoothing_margin(40.0), 10.0);
        assert_eq!(smoothing_margin(4.0), 1.0);
    }

    #[test]
    fn margin_is_zero_for_degenerate_widths() {
        assert_eq!(smoothing_margin(0.0), 0.0);
        assert_eq!(smoothing_margin(-5.0), 0.0); // defensive: malformed negative width
    }

    #[test]
    fn smoothstep_is_centered_and_matches_hermite_formula() {
        assert_eq!(smoothstep(0.0), 0.0);
        assert_eq!(smoothstep(1.0), 1.0);
        assert_eq!(smoothstep(0.5), 0.5); // 3*0.25 - 2*0.125 = 0.75 - 0.25 = 0.5
        // Out-of-range t is clamped rather than extrapolated.
        assert_eq!(smoothstep(-1.0), 0.0);
        assert_eq!(smoothstep(2.0), 1.0);
    }

    #[test]
    fn ascending_band_boundary_blends_symmetrically_around_the_boundary() {
        // Two contiguous 1000 fps-wide bands: margin = min(50, 0.25*1000) = 50, half = 25.
        let segments = vec![
            BCSegmentData {
                velocity_min: 0.0,
                velocity_max: 1000.0,
                bc_value: 0.25,
            },
            BCSegmentData {
                velocity_min: 1000.0,
                velocity_max: 2000.0,
                bc_value: 0.75,
            },
        ];

        // Deep mid-band: exactly flat, bit-identical to the plain band value.
        assert_eq!(velocity_segment_bc(200.0, &segments, 0.9).to_bits(), 0.25f64.to_bits());
        assert_eq!(velocity_segment_bc(1800.0, &segments, 0.9).to_bits(), 0.75f64.to_bits());

        // Exactly at the boundary (t = 0.5): the midpoint of the two band values.
        assert_eq!(velocity_segment_bc(1000.0, &segments, 0.9), 0.5);

        // At the low edge of the margin window (t = 0): matches the lower band exactly.
        assert_eq!(
            velocity_segment_bc(975.0, &segments, 0.9).to_bits(),
            0.25f64.to_bits()
        );
        // At the high edge of the margin window (t = 1): matches the upper band exactly.
        assert_eq!(
            velocity_segment_bc(1025.0, &segments, 0.9).to_bits(),
            0.75f64.to_bits()
        );

        // Strictly inside the window: strictly between the two band values (monotonic).
        let just_below = velocity_segment_bc(999.0, &segments, 0.9);
        let just_above = velocity_segment_bc(1001.0, &segments, 0.9);
        assert!(just_below > 0.25 && just_below < 0.5);
        assert!(just_above > 0.5 && just_above < 0.75);
    }

    #[test]
    fn descending_stored_order_matches_ascending_boundary_blend() {
        // Same table as above, stored high-to-low: the helper's semantics are documented as
        // order-independent, and MBA-1404 must not break that for the new blend either.
        let ascending = vec![
            BCSegmentData {
                velocity_min: 0.0,
                velocity_max: 1000.0,
                bc_value: 0.25,
            },
            BCSegmentData {
                velocity_min: 1000.0,
                velocity_max: 2000.0,
                bc_value: 0.75,
            },
        ];
        let mut descending = ascending.clone();
        descending.reverse();

        for v in [200.0, 975.0, 999.0, 1000.0, 1001.0, 1025.0, 1800.0] {
            assert_eq!(
                velocity_segment_bc(v, &ascending, 0.9),
                velocity_segment_bc(v, &descending, 0.9),
                "order must not affect the blended result at v={v}"
            );
        }
    }

    #[test]
    fn gapped_table_blends_both_edges_of_the_fallback_gap_and_stays_flat_mid_gap() {
        // Band widths 1000 each; gap from 1000 to 1200 (width 200). Exit-edge margin =
        // min(50, 0.25*min(1000,200)) = 50; entry-edge margin = min(50, 0.25*min(200,1000)) = 50.
        let segments = vec![
            BCSegmentData {
                velocity_min: 0.0,
                velocity_max: 1000.0,
                bc_value: 0.25,
            },
            BCSegmentData {
                velocity_min: 1200.0,
                velocity_max: 2200.0,
                bc_value: 0.75,
            },
        ];
        let fallback = 0.5; // exact binary fraction: keeps the blend math exact here too.

        // Deep mid-gap: exactly the fallback, untouched by either edge's margin.
        assert_eq!(
            velocity_segment_bc(1100.0, &segments, fallback).to_bits(),
            fallback.to_bits()
        );

        // Exit edge (band -> gap) at v=1000: t=0.5 blend between the band value and fallback.
        assert_eq!(velocity_segment_bc(1000.0, &segments, fallback), 0.375); // (0.25+0.5)/2

        // Entry edge (gap -> band) at v=1200: t=0.5 blend between fallback and the band value.
        assert_eq!(velocity_segment_bc(1200.0, &segments, fallback), 0.625); // (0.5+0.75)/2

        // Just outside each margin window: exactly flat again.
        assert_eq!(
            velocity_segment_bc(974.0, &segments, fallback).to_bits(),
            0.25f64.to_bits()
        );
        assert_eq!(
            velocity_segment_bc(1226.0, &segments, fallback).to_bits(),
            0.75f64.to_bits()
        );
    }

    #[test]
    fn single_segment_table_never_blends_and_is_byte_identical_to_the_raw_lookup() {
        let single = vec![BCSegmentData {
            velocity_min: 1000.0,
            velocity_max: 2000.0,
            bc_value: 0.5,
        }];

        for v in [-1e6, 500.0, 1000.0, 1500.0, 1999.999, 2000.0, 1e6] {
            assert_eq!(
                velocity_segment_bc(v, &single, 0.9).to_bits(),
                raw_velocity_segment_bc(v, &single, 0.9).to_bits(),
                "single-band tables must never blend (v={v})"
            );
        }
    }

    #[test]
    fn empty_table_never_blends_and_always_returns_the_fallback_exactly() {
        let empty: Vec<BCSegmentData> = vec![];
        for v in [-1e6, 0.0, 500.0, 1e6] {
            assert_eq!(velocity_segment_bc(v, &empty, 0.73).to_bits(), 0.73f64.to_bits());
        }
    }

    #[test]
    fn zero_width_adjacent_band_disables_blending_at_its_boundaries() {
        // The middle "band" is degenerate (velocity_min == velocity_max), so it can never be
        // matched directly, but it still touches both of its neighbors' boundaries -- both of
        // those boundaries must fall back to the hard step (margin 0), matching the pre-MBA-1404
        // raw lookup exactly, rather than average toward the unmatchable degenerate value.
        let segments = vec![
            BCSegmentData {
                velocity_min: 0.0,
                velocity_max: 1000.0,
                bc_value: 0.5,
            },
            BCSegmentData {
                velocity_min: 1000.0,
                velocity_max: 1000.0,
                bc_value: 0.9,
            },
            BCSegmentData {
                velocity_min: 1000.0,
                velocity_max: 2000.0,
                bc_value: 0.6,
            },
        ];

        for v in [999.0, 999.99, 1000.0, 1000.01, 1001.0] {
            assert_eq!(
                velocity_segment_bc(v, &segments, 0.99).to_bits(),
                raw_velocity_segment_bc(v, &segments, 0.99).to_bits(),
                "a degenerate adjacent band must disable blending, not average toward it (v={v})"
            );
        }
    }

    #[test]
    fn coverage_entry_and_exit_edges_stay_flat_since_the_clamp_matches_the_bordering_band() {
        // Below the lowest band and above the highest band, the raw lookup clamps to that
        // same band's own value, so the "coverage entry/exit" boundary is a no-op blend
        // (lo == hi): it must be exactly flat, not just close, all the way up to (and past)
        // the boundary's own margin window.
        let segments = vec![
            BCSegmentData {
                velocity_min: 1000.0,
                velocity_max: 2000.0,
                bc_value: 0.25,
            },
            BCSegmentData {
                velocity_min: 2000.0,
                velocity_max: 3000.0,
                bc_value: 0.75,
            },
        ];

        for v in [-1e6, -100.0, 999.0, 1000.0, 1001.0] {
            assert_eq!(
                velocity_segment_bc(v, &segments, 0.5).to_bits(),
                0.25f64.to_bits(),
                "below-coverage clamp must stay exactly flat (v={v})"
            );
        }
        for v in [2999.0, 3000.0, 3001.0, 1e6] {
            assert_eq!(
                velocity_segment_bc(v, &segments, 0.5).to_bits(),
                0.75f64.to_bits(),
                "above-coverage clamp must stay exactly flat (v={v})"
            );
        }
    }
}
