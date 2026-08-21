//! Transport-free tall-target service.
//!
//! Replicates the CLI's `tall-target` arithmetic exactly — same imperial/metric divisor
//! branch, same four validation guards, same tracking-CF band check via
//! [`crate::profile::validate_tracking_cf`] — as a versioned request/result pair the bridge
//! can expose to embedded consumers. The CLI's `Commands::TallTarget` arm keeps its
//! `println!` formatting but calls [`tall_target_v1`] for the arithmetic, so the two callers
//! cannot drift apart.
//!
//! No feature gate: must compile for wasm32-unknown-unknown. Depends only on
//! `crate::adjustment` and `crate::profile`, both unconditional and filesystem-free.

use serde::{Deserialize, Serialize};

use crate::adjustment::{drop_to_adjustment, AdjustmentUnit};

/// Request for a tall-target correction-factor computation.
///
/// `metric` selects the measured-travel unit exactly like the CLI's `--units` flag:
/// `false` measures `measured` in inches at yards (divisor 36.0), `true` measures it in
/// centimeters at meters (divisor 100.0).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TallTargetRequestV1 {
    pub dialed: f64,
    pub measured: f64,
    pub range: f64,
    pub unit: AdjustmentUnit,
    pub metric: bool,
}

/// Result of a tall-target correction-factor computation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct TallTargetResultV1 {
    pub dialed: f64,
    pub actual: f64,
    pub correction_factor: f64,
    pub within_accepted_band: bool,
}

/// Errors returned by [`tall_target_v1`].
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum TallTargetErrorV1 {
    #[error("invalid tall-target request: {0}")]
    InvalidInput(String),
}

/// Computes the tall-target correction factor: same arithmetic and validation guards as the
/// CLI's `tall-target` command, minus the printing.
pub fn tall_target_v1(req: &TallTargetRequestV1) -> Result<TallTargetResultV1, TallTargetErrorV1> {
    if req.unit == AdjustmentUnit::Clicks {
        return Err(TallTargetErrorV1::InvalidInput(
            "clicks is not an angular unit; enter the dialed travel in mil, moa, smoa, or iphy"
                .to_string(),
        ));
    }
    if !req.dialed.is_finite() || req.dialed <= 0.0 {
        return Err(TallTargetErrorV1::InvalidInput(
            "dialed must be a positive angular travel".to_string(),
        ));
    }
    if !req.measured.is_finite() || req.measured <= 0.0 {
        return Err(TallTargetErrorV1::InvalidInput(
            "measured must be a positive measured travel".to_string(),
        ));
    }
    if !req.range.is_finite() || req.range < 1.0 {
        return Err(TallTargetErrorV1::InvalidInput(
            "range must be at least 1 yard/meter".to_string(),
        ));
    }
    let drop_len = if req.metric {
        req.measured / 100.0
    } else {
        req.measured / 36.0
    };
    let actual = drop_to_adjustment(drop_len, req.range, req.unit);
    let correction_factor = actual / req.dialed;
    Ok(TallTargetResultV1 {
        dialed: req.dialed,
        actual,
        correction_factor,
        within_accepted_band: crate::profile::validate_tracking_cf(
            correction_factor,
            "the computed factor",
        )
        .is_ok(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tall_target_matches_the_cli_arithmetic_and_rejects_clicks() {
        // 30 inches measured at 100 yd against 10 mil dialed.
        let req = TallTargetRequestV1 {
            dialed: 10.0,
            measured: 30.0,
            range: 100.0,
            unit: AdjustmentUnit::Mil,
            metric: false,
        };
        let r = tall_target_v1(&req).expect("computes");
        let expected_actual =
            crate::adjustment::drop_to_adjustment(30.0 / 36.0, 100.0, AdjustmentUnit::Mil);
        assert!((r.actual - expected_actual).abs() < 1e-12);
        assert!((r.correction_factor - expected_actual / 10.0).abs() < 1e-12);

        // Clicks is not an angular unit.
        let bad = TallTargetRequestV1 {
            unit: AdjustmentUnit::Clicks,
            ..req
        };
        assert!(tall_target_v1(&bad).is_err());

        // Non-positive and sub-1 range are rejected, matching the CLI guards.
        assert!(tall_target_v1(&TallTargetRequestV1 {
            dialed: 0.0,
            ..req
        })
        .is_err());
        assert!(tall_target_v1(&TallTargetRequestV1 {
            measured: -1.0,
            ..req
        })
        .is_err());
        assert!(tall_target_v1(&TallTargetRequestV1 {
            range: 0.5,
            ..req
        })
        .is_err());
    }

    #[test]
    fn service_reproduces_the_cli_imperial_and_metric_ratios() {
        // Imperial: inches at yards, 36 in/yd. Metric: cm at metres, 100 cm/m.
        for (metric, measured, divisor) in [(false, 30.0, 36.0), (true, 76.2, 100.0)] {
            let r = tall_target_v1(&TallTargetRequestV1 {
                dialed: 10.0,
                measured,
                range: 100.0,
                unit: AdjustmentUnit::Mil,
                metric,
            })
            .expect("computes");
            let expected = crate::adjustment::drop_to_adjustment(
                measured / divisor,
                100.0,
                AdjustmentUnit::Mil,
            );
            assert!((r.actual - expected).abs() < 1e-12, "metric={metric}");
        }
    }
}
