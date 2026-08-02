//! The single input taxonomy for the decision-support train.
//!
//! Two levels: `InputAxis` leaves (what MBA-1347 propagates individually) grouped
//! into `InputGroup` buckets (what MBA-1345 attributes by). Defined once so the two
//! features cannot drift apart.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum InputGroup {
    ProjectileDrag, MuzzleVelocity, ZeroSightGeometry,
    Atmosphere, Wind, ShotGeometry, Effects,
}

impl InputGroup {
    pub const ALL: &'static [InputGroup] = &[
        InputGroup::ProjectileDrag, InputGroup::MuzzleVelocity, InputGroup::ZeroSightGeometry,
        InputGroup::Atmosphere, InputGroup::Wind, InputGroup::ShotGeometry, InputGroup::Effects,
    ];
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum InputAxis {
    Mass, Diameter, Length, BallisticCoefficient,
    MuzzleVelocityMps,
    SightHeight, ZeroDistance, ZeroPoiUp, ZeroPoiRight, SightOffsetLateral, MuzzleAngle,
    Altitude, Temperature, Pressure, RelativeHumidity, Latitude,
    WindSpeed, WindDirection, WindVertical,
    TargetDistance, ShootingAngle, Cant, ShotAzimuth, AimAzimuth,
    MagnusEnabled, CoriolisEnabled, EnhancedSpinDriftEnabled,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AxisKind {
    /// A real-valued input. `default_rel_step`/`min_abs_step` follow the existing
    /// convention h = (|x| * rel).max(min_abs) from src/truing.rs:1266.
    Continuous { unit: &'static str, default_rel_step: f64, min_abs_step: f64 },
    /// A boolean or enumerated input. Participates in counterfactuals; never differentiated.
    Categorical,
}

#[derive(Debug, Clone, Copy)]
pub struct AxisMeta {
    pub group: InputGroup,
    pub kind: AxisKind,
    /// True when changing this axis invalidates the zero, so the perturbed solve must
    /// re-zero before it is comparable. This flag is the cost control for the kernel.
    pub requires_rezero: bool,
}

const fn cont(unit: &'static str, rel: f64, min_abs: f64) -> AxisKind {
    AxisKind::Continuous { unit, default_rel_step: rel, min_abs_step: min_abs }
}

pub fn axis_meta(axis: InputAxis) -> AxisMeta {
    use InputAxis::*;
    use InputGroup::*;
    let (group, kind, requires_rezero) = match axis {
        Mass                 => (ProjectileDrag, cont("kg", 1e-3, 1e-7), true),
        Diameter             => (ProjectileDrag, cont("m", 1e-3, 1e-7), true),
        Length               => (ProjectileDrag, cont("m", 1e-3, 1e-6), false),
        BallisticCoefficient => (ProjectileDrag, cont("", 1e-3, 1e-4), true),
        MuzzleVelocityMps    => (MuzzleVelocity, cont("m/s", 1e-3, 0.15), true),
        SightHeight          => (ZeroSightGeometry, cont("m", 1e-3, 1e-5), true),
        ZeroDistance         => (ZeroSightGeometry, cont("m", 1e-3, 0.1), true),
        ZeroPoiUp            => (ZeroSightGeometry, cont("m", 1e-3, 1e-5), true),
        ZeroPoiRight         => (ZeroSightGeometry, cont("m", 1e-3, 1e-5), true),
        SightOffsetLateral   => (ZeroSightGeometry, cont("m", 1e-3, 1e-5), true),
        MuzzleAngle          => (ZeroSightGeometry, cont("rad", 1e-3, 1e-7), false),
        Altitude             => (Atmosphere, cont("m", 1e-3, 1.0), false),
        Temperature          => (Atmosphere, cont("K", 1e-3, 0.05), false),
        Pressure             => (Atmosphere, cont("Pa", 1e-3, 10.0), false),
        RelativeHumidity     => (Atmosphere, cont("", 1e-3, 1e-3), false),
        Latitude             => (Atmosphere, cont("rad", 1e-3, 1e-5), false),
        WindSpeed            => (Wind, cont("m/s", 1e-3, 0.05), false),
        WindDirection        => (Wind, cont("rad", 1e-3, 1e-4), false),
        WindVertical         => (Wind, cont("m/s", 1e-3, 0.05), false),
        TargetDistance       => (ShotGeometry, cont("m", 1e-3, 0.5), false),
        ShootingAngle        => (ShotGeometry, cont("rad", 1e-3, 1e-4), false),
        Cant                 => (ShotGeometry, cont("rad", 1e-3, 1e-4), false),
        ShotAzimuth          => (ShotGeometry, cont("rad", 1e-3, 1e-4), false),
        AimAzimuth           => (ShotGeometry, cont("rad", 1e-3, 1e-4), false),
        MagnusEnabled            => (Effects, AxisKind::Categorical, false),
        CoriolisEnabled          => (Effects, AxisKind::Categorical, false),
        EnhancedSpinDriftEnabled => (Effects, AxisKind::Categorical, false),
    };
    AxisMeta { group, kind, requires_rezero }
}

pub fn axes_in_group(group: InputGroup) -> &'static [InputAxis] {
    use InputAxis::*;
    match group {
        InputGroup::ProjectileDrag => &[Mass, Diameter, Length, BallisticCoefficient],
        InputGroup::MuzzleVelocity => &[MuzzleVelocityMps],
        InputGroup::ZeroSightGeometry =>
            &[SightHeight, ZeroDistance, ZeroPoiUp, ZeroPoiRight, SightOffsetLateral, MuzzleAngle],
        InputGroup::Atmosphere => &[Altitude, Temperature, Pressure, RelativeHumidity, Latitude],
        InputGroup::Wind => &[WindSpeed, WindDirection, WindVertical],
        InputGroup::ShotGeometry =>
            &[TargetDistance, ShootingAngle, Cant, ShotAzimuth, AimAzimuth],
        InputGroup::Effects =>
            &[MagnusEnabled, CoriolisEnabled, EnhancedSpinDriftEnabled],
    }
}

// KNOWN LIMITATIONS: Two axis/mode combinations are physically wrong for counterfactuals and
// need guards in a later task (not this one — this task is data only):
//
// (a) Altitude when the original request used QNH pressure: the rebuilt request carries an
//     absolute station pressure, so perturbing altitude changes density-by-altitude but NOT
//     station pressure, which is the opposite of what a QNH-entering user means.
//
// (b) ShotAzimuth when the original request used compass-referenced wind: the rebuilt request
//     carries shooter-relative wind, so the wind rotates WITH the rifle instead of staying
//     earth-fixed — the counterfactual is physically inverted.

#[cfg(test)]
mod tests {
    use super::*;

    /// Every axis belongs to exactly one group, and every group lists it.
    #[test]
    fn taxonomy_is_a_partition() {
        let mut seen = Vec::new();
        for g in InputGroup::ALL {
            for a in axes_in_group(*g) {
                assert_eq!(axis_meta(*a).group, *g, "{a:?} listed under the wrong group");
                assert!(!seen.contains(a), "{a:?} appears in more than one group");
                seen.push(*a);
            }
        }
        assert!(seen.contains(&InputAxis::MuzzleVelocityMps));
        assert!(seen.contains(&InputAxis::Cant));
    }

    /// Effects are boolean toggles and must never be differentiated (spec D7).
    #[test]
    fn effects_are_categorical() {
        for a in axes_in_group(InputGroup::Effects) {
            assert!(matches!(axis_meta(*a).kind, AxisKind::Categorical),
                    "{a:?} must be Categorical");
        }
    }

    /// Cost control: only zero-affecting axes force a re-zero.
    #[test]
    fn only_zero_affecting_axes_require_rezero() {
        assert!(axis_meta(InputAxis::SightHeight).requires_rezero);
        assert!(axis_meta(InputAxis::ZeroDistance).requires_rezero);
        assert!(!axis_meta(InputAxis::WindSpeed).requires_rezero);
        assert!(!axis_meta(InputAxis::TargetDistance).requires_rezero);
    }
}
