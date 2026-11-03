//! The two 0.49.0 names 0.50.0 keeps for one release as `#[deprecated]` shims instead of
//! removing them: the `BallisticInputs::humidity_percent()` accessor, which returns the percent
//! field, and the 24-argument `wez::compute_wez`, which flies the 50 % default every 0.49.0
//! sweep flew.

#![allow(deprecated)]

use ballistics_engine::wez::{compute_wez, compute_wez_with_humidity, TargetSizeMetric};
use ballistics_engine::{BallisticInputs, DragModel, RelativeHumidity};

fn with_humidity(percent: f64) -> BallisticInputs {
    BallisticInputs {
        humidity_percent: percent,
        ..Default::default()
    }
}

/// The accessor returns the percent field as stored: 0.50.0 clamps a humidity nowhere, so a
/// value beyond the sensor band comes back as it is, for the solve to refuse.
#[test]
fn the_deprecated_humidity_accessor_returns_the_percent_field() {
    for percent in [0.0, 37.5, 100.0, 102.0, -3.0, 150.0, -0.0] {
        let read = with_humidity(percent).humidity_percent();
        assert_eq!(read.to_bits(), percent.to_bits(), "{percent} %");
    }
    assert!(with_humidity(f64::NAN).humidity_percent().is_nan());
    assert_eq!(BallisticInputs::default().humidity_percent(), 50.0);
}

/// From 0 to 100 it reads what 0.49.0's accessor read for the same air (the 0-1 fraction field
/// as a percent), bit for bit.
#[cfg(not(target_arch = "wasm32"))]
#[test]
fn the_deprecated_humidity_accessor_reads_what_0_49_0_read_from_0_to_100() {
    for fraction in [0.0, 0.25, 0.375, 0.5, 0.75, 1.0] {
        let then = ballistics_engine_049::BallisticInputs {
            humidity: fraction,
            ..Default::default()
        };
        assert_eq!(
            with_humidity(fraction * 100.0).humidity_percent().to_bits(),
            then.humidity_percent().to_bits(),
            "a fraction of {fraction}"
        );
    }
}

/// A short WEZ sweep through `run` (the deprecated form, or the new one at a humidity), as its
/// JSON: every float written at its shortest round-trip, so equal text is equal bits.
fn sweep(
    run: impl Fn() -> Result<ballistics_engine::wez::WezResult, Box<dyn std::error::Error>>,
) -> String {
    serde_json::to_string(&run().expect("the sweep runs")).unwrap()
}

macro_rules! wez_arguments {
    ($function:ident $(, $humidity:expr)?) => {
        $function(
            823.0,
            0.0,
            0.243,
            0.0113,
            0.00782,
            20,
            5.0,
            0.0001,
            0.005,
            1.0,
            0.05,
            3.0,
            90.0,
            0.0,
            1.5,
            TargetSizeMetric::Rect {
                width_m: 0.5,
                height_m: 0.75,
            },
            300.0,
            400.0,
            100.0,
            DragModel::G7,
            None,
            1.0,
            0.0,
            0.0,
            $($humidity)?
        )
    };
}

/// The 24-argument form is the 25-argument one at the default humidity, to the bit, and the
/// humidity argument is not ignored.
#[test]
fn the_deprecated_compute_wez_flies_the_default_humidity() {
    let deprecated = sweep(|| wez_arguments!(compute_wez));
    let at_default = sweep(|| wez_arguments!(compute_wez_with_humidity, RelativeHumidity::DEFAULT));
    assert_eq!(deprecated, at_default);
    let saturated = sweep(|| {
        wez_arguments!(
            compute_wez_with_humidity,
            RelativeHumidity::try_percent(100.0).unwrap()
        )
    });
    assert_ne!(saturated, at_default, "the humidity reaches the sweep");
}
