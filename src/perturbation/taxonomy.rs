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
    Mass, Diameter, Length, BallisticCoefficient, TwistRate, TwistDirection, DragModel,
    MuzzleVelocityMps,
    SightHeight, ZeroDistance, ZeroPoiUp, ZeroPoiRight, SightOffsetLateral, MuzzleHeight, MuzzleAngle,
    Altitude, Temperature, Pressure, RelativeHumidity, Latitude,
    WindSpeed, WindDirection, WindVertical,
    TargetDistance, ShootingAngle, Cant, ShotAzimuth, AimAzimuth, TargetHeight,
    MagnusEnabled, CoriolisEnabled, EnhancedSpinDriftEnabled,
}

impl InputAxis {
    pub const ALL: &'static [InputAxis] = &[
        InputAxis::Mass, InputAxis::Diameter, InputAxis::Length, InputAxis::BallisticCoefficient,
        InputAxis::TwistRate, InputAxis::TwistDirection, InputAxis::DragModel,
        InputAxis::MuzzleVelocityMps,
        InputAxis::SightHeight, InputAxis::ZeroDistance, InputAxis::ZeroPoiUp, InputAxis::ZeroPoiRight,
        InputAxis::SightOffsetLateral, InputAxis::MuzzleHeight, InputAxis::MuzzleAngle,
        InputAxis::Altitude, InputAxis::Temperature, InputAxis::Pressure, InputAxis::RelativeHumidity,
        InputAxis::Latitude,
        InputAxis::WindSpeed, InputAxis::WindDirection, InputAxis::WindVertical,
        InputAxis::TargetDistance, InputAxis::ShootingAngle, InputAxis::Cant, InputAxis::ShotAzimuth,
        InputAxis::AimAzimuth, InputAxis::TargetHeight,
        InputAxis::MagnusEnabled, InputAxis::CoriolisEnabled, InputAxis::EnhancedSpinDriftEnabled,
    ];
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
        TwistRate            => (ProjectileDrag, cont("m/turn", 1e-3, 1e-5), false),
        TwistDirection       => (ProjectileDrag, AxisKind::Categorical, false),
        DragModel            => (ProjectileDrag, AxisKind::Categorical, true),
        MuzzleVelocityMps    => (MuzzleVelocity, cont("m/s", 1e-3, 0.15), true),
        SightHeight          => (ZeroSightGeometry, cont("m", 1e-3, 1e-5), true),
        ZeroDistance         => (ZeroSightGeometry, cont("m", 1e-3, 0.1), true),
        ZeroPoiUp            => (ZeroSightGeometry, cont("m", 1e-3, 1e-5), true),
        ZeroPoiRight         => (ZeroSightGeometry, cont("m", 1e-3, 1e-5), true),
        SightOffsetLateral   => (ZeroSightGeometry, cont("m", 1e-3, 1e-5), true),
        MuzzleHeight         => (ZeroSightGeometry, cont("m", 1e-3, 1e-5), true),
        MuzzleAngle          => (ZeroSightGeometry, cont("rad", 1e-3, 1e-4), false),
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
        TargetHeight         => (ShotGeometry, cont("m", 1e-3, 1e-3), false),
        MagnusEnabled            => (Effects, AxisKind::Categorical, false),
        CoriolisEnabled          => (Effects, AxisKind::Categorical, false),
        EnhancedSpinDriftEnabled => (Effects, AxisKind::Categorical, false),
    };
    AxisMeta { group, kind, requires_rezero }
}

pub fn axes_in_group(group: InputGroup) -> &'static [InputAxis] {
    use InputAxis::*;
    match group {
        InputGroup::ProjectileDrag => &[Mass, Diameter, Length, BallisticCoefficient, TwistRate, TwistDirection, DragModel],
        InputGroup::MuzzleVelocity => &[MuzzleVelocityMps],
        InputGroup::ZeroSightGeometry =>
            &[SightHeight, ZeroDistance, ZeroPoiUp, ZeroPoiRight, SightOffsetLateral, MuzzleHeight, MuzzleAngle],
        InputGroup::Atmosphere => &[Altitude, Temperature, Pressure, RelativeHumidity, Latitude],
        InputGroup::Wind => &[WindSpeed, WindDirection, WindVertical],
        InputGroup::ShotGeometry =>
            &[TargetDistance, ShootingAngle, Cant, ShotAzimuth, AimAzimuth, TargetHeight],
        InputGroup::Effects =>
            &[MagnusEnabled, CoriolisEnabled, EnhancedSpinDriftEnabled],
    }
}

// KNOWN LIMITATIONS: Several axis/mode combinations require guards in a later task (not this one —
// this task is data only):
//
// (a) Altitude when the original request used QNH pressure: the rebuilt request carries an
//     absolute station pressure, so perturbing altitude changes density-by-altitude but NOT
//     station pressure, which is the opposite of what a QNH-entering user means.
//
// (b) ShotAzimuth when the original request used compass-referenced wind: the rebuilt request
//     carries shooter-relative wind, so the wind rotates WITH the rifle instead of staying
//     earth-fixed — the counterfactual is physically inverted.
//
// (c) WindSpeed/WindDirection/WindVertical when the original request used segmented wind: there
//     is no single scalar to read or perturb. These axes return None from read_axis under
//     segmented wind and the kernel treats them as absent.
//
// (d) Magnus + EnhancedSpinDrift together: validate_effects (src/solve_json.rs:1544-1552)
//     rejects magnus: true + enhanced_spin_drift: true as ConflictingFields. A future with_axis
//     that flips one to true while the other is already true produces a request that fails
//     validation. This cross-category constraint cannot be expressed purely in the taxonomy.

#[cfg(test)]
mod tests {
    use super::*;

    /// Every axis belongs to exactly one group, and every group lists it.
    /// Set equality: InputAxis::ALL == union of axes_in_group over all groups.
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
        // Set equality: every axis in ALL is listed exactly once, and every listed axis is in ALL.
        assert_eq!(seen.len(), InputAxis::ALL.len(),
                   "axes_in_group() covers {} axes; InputAxis::ALL has {}",
                   seen.len(), InputAxis::ALL.len());
        for axis in InputAxis::ALL {
            assert!(seen.contains(axis), "{axis:?} is in InputAxis::ALL but not in any group");
        }
        // Also verify axis_meta is callable for every member of ALL.
        for axis in InputAxis::ALL {
            let _ = axis_meta(*axis);
        }
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
