//! What 0.50.0 makes of every legacy kernel 4-tuple 0.49.0 flew.
//!
//! 0.49.0's kernels took the air as a bare `(f64, f64, f64, f64)`: `(density, sound speed, 0,
//! 0)` in direct mode, `(altitude, temperature, pressure, density ratio)` otherwise. 0.50.0
//! replaces the tuple with [`KernelAtmosphere`], and the documented one-line migration is
//! [`KernelAtmosphere::from_legacy_params`]. Some tuples 0.49.0 flew are now refused, by design;
//! this file lists every such class, so a later narrowing of the migration is a deliberate edit
//! here and in the CHANGELOG, never a silent one.
//!
//! The "0.49.0" column was measured, not assumed: each tuple below was run through 0.49.0's own
//! `compute_derivatives` (a scratch crate on a `git archive` of 261c1eb, 800 m/s level flight,
//! `BallisticInputs::default()`, BC 0.442), and every tuple in [`FLOWN_BY_0_49_0`] gave a finite
//! acceleration there, every one in [`NAN_IN_0_49_0`] a NaN.
//!
//! Exact vacuum is the one class 0.50.0 keeps that a strict density bound would have dropped:
//! `uniform(0, c)` is direct mode with no air at all, the 0.49.0 way to switch drag off.

use ballistics_engine::atmosphere::AtmosphereErrorReason;
use ballistics_engine::derivatives::compute_derivatives;
use ballistics_engine::trajectory_integration::{try_integrate_trajectory, TrajectoryParams};
use ballistics_engine::{BallisticInputs, DragModel, KernelAtmosphere, RelativeHumidity};
use nalgebra::Vector3;

/// What `from_legacy_params` does with a tuple in 0.50.0.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Now {
    /// Converts, and flies a finite acceleration.
    Flies,
    /// Refused, naming this parameter, for this reason.
    Refused(&'static str, AtmosphereErrorReason),
}

/// A legacy kernel atmosphere 4-tuple.
type Tuple = (f64, f64, f64, f64);

use AtmosphereErrorReason as R;
use Now::{Flies, Refused};

/// `(what it is, tuple, humidity percent, 0.50.0 outcome)` for every tuple class 0.49.0's
/// `compute_derivatives` flew to a finite acceleration.
const FLOWN_BY_0_49_0: &[(&str, Tuple, f64, Now)] = &[
    // Direct mode, (density, sound speed, 0, 0).
    ("direct air", (1.0, 340.29, 0.0, 0.0), 0.0, Flies),
    (
        "direct vacuum, density 0",
        (0.0, 340.29, 0.0, 0.0),
        0.0,
        Flies,
    ),
    (
        "direct vacuum, density -0",
        (-0.0, 340.29, 0.0, 0.0),
        0.0,
        Flies,
    ),
    // 0.49.0 flew a negative density as a FORWARD drag of +446.9 m/s^2 at 800 m/s.
    (
        "direct negative density",
        (-1.0, 340.29, 0.0, 0.0),
        0.0,
        Refused("density_kg_m3", R::UniformDensityOutOfRange),
    ),
    // 0.49.0 flew an infinite sound speed as Mach 0.
    (
        "direct infinite sound speed",
        (1.0, f64::INFINITY, 0.0, 0.0),
        0.0,
        Refused("speed_of_sound_mps", R::NotFinite),
    ),
    // Standard mode, (altitude, temperature, pressure, ratio).
    ("standard air", (0.0, 15.0, 1013.25, 1.0), 50.0, Flies),
    // A finite ratio <= 0 means "not supplied" and flies 1.0, as in 0.49.0 (MBA-1157).
    ("standard ratio 0", (0.0, 15.0, 1013.25, 0.0), 50.0, Flies),
    ("standard ratio -1", (0.0, 15.0, 1013.25, -1.0), 50.0, Flies),
    (
        "standard ratio 1.999",
        (0.0, 15.0, 1013.25, 1.999),
        0.0,
        Flies,
    ),
    // 0.49.0 flew a NaN or -inf ratio as 1.0.
    (
        "standard ratio NaN",
        (0.0, 15.0, 1013.25, f64::NAN),
        50.0,
        Refused("density_ratio", R::NotFinite),
    ),
    (
        "standard ratio -inf",
        (0.0, 15.0, 1013.25, f64::NEG_INFINITY),
        50.0,
        Refused("density_ratio", R::NotFinite),
    ),
    // 0.49.0's compute_derivatives flew any ratio above 0 (its fast_integrate refused >= 2).
    (
        "standard ratio 2",
        (0.0, 15.0, 1013.25, 2.0),
        0.0,
        Refused("density_ratio", R::DensityRatioOutOfRange),
    ),
    (
        "standard ratio 2.5",
        (0.0, 15.0, 1013.25, 2.5),
        0.0,
        Refused("density_ratio", R::DensityRatioOutOfRange),
    ),
    // 0.49.0 flew a negative station pressure on its density ratio alone.
    (
        "standard negative pressure",
        (0.0, 15.0, -5.0, 1.0),
        0.0,
        Refused("pressure_hpa", R::NonPositivePressure),
    ),
    (
        "standard temperature below absolute zero",
        (0.0, -300.0, 1013.25, 1.0),
        0.0,
        Refused("temperature_c", R::AtOrBelowAbsoluteZero),
    ),
    // 0.49.0 read the NaN saturation pressure above 647.096 K as pure steam, at any humidity.
    (
        "standard temperature above water's critical point",
        (0.0, 400.0, 1013.25, 1.0),
        0.0,
        Refused("temperature_c", R::ImpossibleAir),
    ),
    // 0.49.0's kernel flew every station dry, so it never saw the humidity; 0.50.0 flies it.
    (
        "standard hot air, dry",
        (0.0, 120.0, 1013.25, 1.0),
        0.0,
        Flies,
    ),
    (
        "standard hot air, saturated: air that cannot exist",
        (0.0, 120.0, 1013.25, 1.0),
        100.0,
        Refused("humidity_percent", R::ImpossibleAir),
    ),
    ("standard cold air", (0.0, -150.0, 1013.25, 1.0), 0.0, Flies),
    (
        "standard high station",
        (5000.0, 15.0, 1013.25, 1.0),
        50.0,
        Flies,
    ),
    (
        "standard huge pressure",
        (0.0, 15.0, 1.0e6, 1.0),
        50.0,
        Flies,
    ),
];

/// Tuples 0.49.0 accepted but flew to a NaN acceleration, and the refusal 0.50.0 gives instead.
const NAN_IN_0_49_0: &[(&str, Tuple, Now)] = &[
    (
        "direct NaN density",
        (f64::NAN, 340.29, 0.0, 0.0),
        Refused("density_kg_m3", R::NotFinite),
    ),
    (
        "direct density 2",
        (2.0, 340.29, 0.0, 0.0),
        Refused("density_kg_m3", R::UniformDensityOutOfRange),
    ),
    (
        "direct density 2.5",
        (2.5, 340.29, 0.0, 0.0),
        Refused("density_kg_m3", R::UniformDensityOutOfRange),
    ),
    (
        "direct sound speed 150",
        (1.0, 150.0, 0.0, 0.0),
        Refused("speed_of_sound_mps", R::UniformSpeedOfSoundOutOfRange),
    ),
    (
        "direct sound speed 200",
        (1.0, 200.0, 0.0, 0.0),
        Refused("speed_of_sound_mps", R::UniformSpeedOfSoundOutOfRange),
    ),
    (
        "standard ratio +inf",
        (0.0, 15.0, 1013.25, f64::INFINITY),
        Refused("density_ratio", R::NotFinite),
    ),
    (
        "standard pressure 0",
        (0.0, 15.0, 0.0, 1.0),
        Refused("pressure_hpa", R::NonPositivePressure),
    ),
    (
        "standard pressure NaN",
        (0.0, 15.0, f64::NAN, 1.0),
        Refused("pressure_hpa", R::NotFinite),
    ),
    (
        "standard pressure +inf",
        (0.0, 15.0, f64::INFINITY, 1.0),
        Refused("pressure_hpa", R::NotFinite),
    ),
    (
        "standard temperature at absolute zero",
        (0.0, -273.15, 1013.25, 1.0),
        Refused("temperature_c", R::AtOrBelowAbsoluteZero),
    ),
    (
        "standard temperature NaN",
        (0.0, f64::NAN, 1013.25, 1.0),
        Refused("temperature_c", R::NotFinite),
    ),
    (
        "standard altitude NaN",
        (f64::NAN, 15.0, 1013.25, 1.0),
        Refused("altitude_m", R::NotFinite),
    ),
    (
        "standard altitude +inf",
        (f64::INFINITY, 15.0, 1013.25, 1.0),
        Refused("altitude_m", R::NotFinite),
    ),
];

fn humidity(percent: f64) -> RelativeHumidity {
    RelativeHumidity::try_percent(percent).expect("a humidity in range")
}

/// 0.49.0's probe: level flight at 800 m/s from the origin.
fn acceleration(atmosphere: &KernelAtmosphere) -> [f64; 3] {
    let d = compute_derivatives(
        Vector3::zeros(),
        Vector3::new(800.0, 0.0, 0.0),
        &BallisticInputs::default(),
        Vector3::zeros(),
        atmosphere,
        0.442,
        None,
        0.0,
        None,
    );
    [d[3], d[4], d[5]]
}

fn outcome(tuple: Tuple, humidity_percent: f64) -> Now {
    match KernelAtmosphere::from_legacy_params(tuple, humidity(humidity_percent)) {
        Ok(atmosphere) => {
            let a = acceleration(&atmosphere);
            assert!(
                a.iter().all(|x| x.is_finite()),
                "{tuple:?} converts but flies a non-finite acceleration {a:?}"
            );
            Flies
        }
        Err(error) => Refused(
            match error.field() {
                "density_kg_m3" => "density_kg_m3",
                "speed_of_sound_mps" => "speed_of_sound_mps",
                "density_ratio" => "density_ratio",
                "pressure_hpa" => "pressure_hpa",
                "temperature_c" => "temperature_c",
                "altitude_m" => "altitude_m",
                "humidity_percent" => "humidity_percent",
                other => panic!("{tuple:?} refused on an unlisted field {other}: {error}"),
            },
            error.reason(),
        ),
    }
}

#[test]
fn every_tuple_class_0_49_0_flew_converts_or_is_refused_as_listed() {
    for &(what, tuple, humidity_percent, expected) in FLOWN_BY_0_49_0 {
        assert_eq!(
            outcome(tuple, humidity_percent),
            expected,
            "{what}: {tuple:?} at {humidity_percent} %"
        );
    }
}

#[test]
fn every_tuple_0_49_0_flew_to_nan_is_now_a_named_refusal() {
    for &(what, tuple, expected) in NAN_IN_0_49_0 {
        assert_eq!(outcome(tuple, 50.0), expected, "{what}: {tuple:?}");
    }
}

#[test]
fn exact_vacuum_is_direct_mode_with_no_drag_at_all() {
    for density in [0.0, -0.0] {
        let vacuum = KernelAtmosphere::uniform(density, 340.29)
            .unwrap_or_else(|e| panic!("uniform({density}, 340.29) is a vacuum: {e}"));
        assert_eq!(
            KernelAtmosphere::from_legacy_params(
                (density, 340.29, 0.0, 0.0),
                RelativeHumidity::DRY
            )
            .expect("the legacy vacuum tuple converts"),
            vacuum
        );
        // 0.49.0 gave exactly [0, -9.80665, 0]: gravity and nothing else.
        let [ax, ay, az] = acceleration(&vacuum);
        assert_eq!(ax, 0.0, "no drag along the flight");
        assert_eq!(az, 0.0, "no drag across it");
        assert_eq!(ay, -9.80665, "gravity alone");
    }
}

#[test]
fn a_vacuum_flight_is_the_closed_form() {
    let params = TrajectoryParams {
        mass_kg: 0.011_34,
        bc: 0.442,
        drag_model: DragModel::G7,
        wind_segments: vec![],
        atmosphere: KernelAtmosphere::uniform(0.0, 340.29).expect("a vacuum"),
        omega_vector: None,
        enable_spin_drift: false,
        enable_magnus: false,
        enable_coriolis: false,
        target_distance_m: 5_000.0,
        enable_wind_shear: false,
        wind_shear_model: "none".into(),
        shooter_altitude_m: 0.0,
        is_twist_right: true,
        shooting_angle: 0.0,
        bullet_diameter: 0.007_823_2,
        bullet_length: 0.031_496,
        twist_rate: 10.0,
        custom_drag_table: None,
        cd_scale: 1.0,
        custom_drag_reference: false,
        bc_segments: None,
        use_bc_segments: false,
        bc_segments_data: None,
        ground_threshold: -1.0e9,
        atmo_sock: None,
    };
    let flight = try_integrate_trajectory(
        [0.0, 0.0, 0.0, 800.0, 0.0, 0.0],
        (0.0, 1.0),
        params,
        "RK4",
        1e-8,
        0.001,
    )
    .expect("a vacuum flies");
    let (t, state) = flight.last().expect("a flight has points");
    assert_eq!(*t, 1.0);
    // x = v t and vx = v exactly; y = -g t^2 / 2 to round-off.
    assert_eq!(state[3], 800.0, "no drag ever slows it");
    assert!((state[0] - 800.0).abs() < 1e-9, "x = {}", state[0]);
    assert!((state[1] + 9.80665 / 2.0).abs() < 1e-9, "y = {}", state[1]);
}

#[test]
fn a_negative_or_nan_density_is_refused_by_name() {
    for density in [-1e-300, -1.0, f64::NEG_INFINITY] {
        let error = KernelAtmosphere::uniform(density, 340.29).unwrap_err();
        assert_eq!(error.field(), "density_kg_m3");
        assert!(
            matches!(error.reason(), R::UniformDensityOutOfRange | R::NotFinite),
            "{density}: {error}"
        );
    }
    let nan = KernelAtmosphere::uniform(f64::NAN, 340.29).unwrap_err();
    assert_eq!((nan.field(), nan.reason()), ("density_kg_m3", R::NotFinite));
    // The message states the closed lower edge: 0 is a vacuum, not a refusal.
    let negative = KernelAtmosphere::uniform(-1.0, 340.29)
        .unwrap_err()
        .to_string();
    assert!(
        negative.contains("0 (a vacuum)"),
        "the refusal states the range that flies: {negative}"
    );
}
