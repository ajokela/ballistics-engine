//! 0.49.0: the presence-aware twins and fields of the library entry points.
//!
//! A station temperature or pressure the caller supplied is never read as omitted. Every
//! entry point that takes a presence-blind `f64` pair keeps the legacy default-sentinel rule
//! under `StationAtmospherePresence::NEITHER`, and that is checked against values captured at
//! ab4eb0c, never against the legacy name alone (which now delegates to the twin, so the two
//! would move together). With `BOTH`, an in-band reading must fly as typed: the result has to
//! lie strictly between its out-of-band neighbours, which every rule flies as typed.

use ballistics_engine::hold_curve::{
    run_sampled_flight, run_sampled_flight_with_station_presence, HoldCurve, HoldCurveLoad,
    SampledFlight,
};
use ballistics_engine::moving_target::calculate_lead_with_station_presence;
use ballistics_engine::wind::parse_wind_segment_str;
use ballistics_engine::wind_scenarios::{
    solve_robust_hold, CorridorLoad, NamedWindScenario, RobustHoldRequest, WindScenarioSetV1,
};
use ballistics_engine::{
    calculate_lead, AtmosphericConditions, BallisticInputs, DragModel, StationAtmospherePresence,
    WindConditions,
};

const NEITHER: StationAtmospherePresence = StationAtmospherePresence::NEITHER;
const BOTH: StationAtmospherePresence = StationAtmospherePresence::BOTH;

fn assert_strictly_between(low: f64, value: f64, high: f64, what: &str) {
    let (lo, hi) = if low < high { (low, high) } else { (high, low) };
    assert!(
        lo < value && value < hi,
        "{what}: {value:?} is not strictly between {low:?} and {high:?}"
    );
}

// ---------------------------------------------------------------------------------------
// hold_curve::run_sampled_flight
// ---------------------------------------------------------------------------------------

fn flight(
    temperature: f64,
    pressure: f64,
    altitude: f64,
    presence: StationAtmospherePresence,
) -> SampledFlight {
    run_sampled_flight_with_station_presence(
        823.0,
        0.475,
        0.0109,
        0.00782,
        DragModel::G1,
        0.0381,
        temperature,
        pressure,
        50.0,
        altitude,
        0.0,
        0.0,
        600.0,
        10.0,
        0.002,
        None,
        None,
        None,
        0.0,
        0.0,
        0.0,
        None,
        presence,
    )
    .unwrap()
}

#[test]
fn sampled_flight_twin_reproduces_ab4eb0c_with_neither() {
    let legacy = run_sampled_flight(
        823.0,
        0.475,
        0.0109,
        0.00782,
        DragModel::G1,
        0.0381,
        20.0,
        1013.0,
        50.0,
        150.0,
        0.0,
        0.0,
        600.0,
        10.0,
        0.002,
        None,
        None,
        None,
        0.0,
        0.0,
        0.0,
        None,
    )
    .unwrap();
    let neither = flight(20.0, 1013.0, 150.0, NEITHER);
    for f in [&legacy, &neither] {
        // captured at ab4eb0c
        assert_eq!(f.samples[50].velocity_mps.to_bits(), 0x4080fcf76d85503f);
        assert_eq!(f.samples[50].drop_m.to_bits(), 0x3ff70fe529243c70);
        assert_eq!(f.reach_m.to_bits(), 0x4082c00000000000);
    }
}

#[test]
fn sampled_flight_twin_flies_a_supplied_in_band_reading() {
    let v = |t, p, alt, presence| flight(t, p, alt, presence).samples[50].velocity_mps;
    let both = v(20.0, 1013.0, 150.0, BOTH);
    assert_strictly_between(
        v(20.0, 1012.0, 150.0, NEITHER),
        both,
        v(20.0, 1014.0, 150.0, NEITHER),
        "pressure",
    );
    assert_ne!(both.to_bits(), v(20.0, 1013.0, 150.0, NEITHER).to_bits());

    let both = v(15.0, 850.0, 1500.0, BOTH);
    assert_strictly_between(
        v(14.8, 850.0, 1500.0, NEITHER),
        both,
        v(15.2, 850.0, 1500.0, NEITHER),
        "temperature",
    );
    assert_ne!(both.to_bits(), v(15.0, 850.0, 1500.0, NEITHER).to_bits());
}

// ---------------------------------------------------------------------------------------
// hold_curve::HoldCurve (HoldCurveLoad::station_presence)
// ---------------------------------------------------------------------------------------

fn hold_load(
    temperature_c: f64,
    pressure_hpa: f64,
    altitude_m: f64,
    station_presence: StationAtmospherePresence,
) -> HoldCurveLoad {
    HoldCurveLoad {
        velocity_mps: 823.0,
        bc: 0.475,
        mass_kg: 0.0109,
        diameter_m: 0.00782,
        drag_model: DragModel::G1,
        sight_height_m: 0.0381,
        zero_distance_m: 100.0,
        temperature_c,
        pressure_hpa,
        humidity: 50.0,
        altitude_m,
        wind_speed_mps: 0.0,
        wind_direction_deg: 0.0,
        station_presence,
    }
}

fn hold_drop_mil(t: f64, p: f64, alt: f64, presence: StationAtmospherePresence) -> f64 {
    HoldCurve::solve(&hold_load(t, p, alt, presence), 700.0)
        .unwrap()
        .at_range(500.0)
        .unwrap()
        .drop_mil
}

#[test]
fn hold_curve_presence_field() {
    // captured at ab4eb0c (the load had no presence then: the legacy rule)
    assert_eq!(
        hold_drop_mil(20.0, 1013.0, 150.0, NEITHER).to_bits(),
        0x400de828b955092f
    );
    let both = hold_drop_mil(20.0, 1013.0, 150.0, BOTH);
    assert_strictly_between(
        hold_drop_mil(20.0, 1012.0, 150.0, NEITHER),
        both,
        hold_drop_mil(20.0, 1014.0, 150.0, NEITHER),
        "pressure",
    );
    let both = hold_drop_mil(15.0, 850.0, 1500.0, BOTH);
    assert_strictly_between(
        hold_drop_mil(14.8, 850.0, 1500.0, NEITHER),
        both,
        hold_drop_mil(15.2, 850.0, 1500.0, NEITHER),
        "temperature",
    );
}

// ---------------------------------------------------------------------------------------
// moving_target::calculate_lead
// ---------------------------------------------------------------------------------------

fn rifle() -> BallisticInputs {
    BallisticInputs {
        muzzle_velocity: 823.0,
        bc_value: 0.475,
        bc_type: DragModel::G1,
        bullet_mass: 0.0109,
        bullet_diameter: 0.00782,
        sight_height: 0.0381,
        ..Default::default()
    }
}

fn air(temperature: f64, pressure: f64, altitude: f64) -> AtmosphericConditions {
    AtmosphericConditions {
        temperature,
        pressure,
        humidity: 50.0,
        altitude,
    }
}

fn lead_m(t: f64, p: f64, alt: f64, presence: StationAtmospherePresence) -> f64 {
    calculate_lead_with_station_presence(
        rifle(),
        WindConditions::default(),
        air(t, p, alt),
        4.0,
        90.0,
        500.0,
        presence,
    )
    .unwrap()
    .lead_m
}

#[test]
fn lead_twin() {
    let legacy = calculate_lead(
        rifle(),
        WindConditions::default(),
        air(20.0, 1013.0, 150.0),
        4.0,
        90.0,
        500.0,
    )
    .unwrap();
    // captured at ab4eb0c
    assert_eq!(legacy.lead_m.to_bits(), 0x4007f2c005a107e0);
    assert_eq!(legacy.time_of_flight_s.to_bits(), 0x3fe7f2c005a107e0);
    assert_eq!(
        lead_m(20.0, 1013.0, 150.0, NEITHER).to_bits(),
        0x4007f2c005a107e0
    );

    let both = lead_m(20.0, 1013.0, 150.0, BOTH);
    assert_strictly_between(
        lead_m(20.0, 1012.0, 150.0, NEITHER),
        both,
        lead_m(20.0, 1014.0, 150.0, NEITHER),
        "pressure",
    );
    let both = lead_m(15.0, 850.0, 1500.0, BOTH);
    assert_strictly_between(
        lead_m(14.8, 850.0, 1500.0, NEITHER),
        both,
        lead_m(15.2, 850.0, 1500.0, NEITHER),
        "temperature",
    );
}

// ---------------------------------------------------------------------------------------
// wind_scenarios::solve_robust_hold (CorridorLoad::station_presence)
// ---------------------------------------------------------------------------------------

fn corridor_elevation_mil(t: f64, p: f64, alt: f64, presence: StationAtmospherePresence) -> f64 {
    let scenario = |name: &str, token: &str| NamedWindScenario {
        name: name.to_string(),
        segments: vec![parse_wind_segment_str(token, true).unwrap()],
    };
    let request = RobustHoldRequest {
        scenarios: WindScenarioSetV1 {
            version: 1,
            scenarios: vec![scenario("low", "4:90:1000"), scenario("high", "14:90:1000")],
            nominal: None,
        },
        ranges_m: vec![182.88, 365.76, 548.64],
        target: None,
        load: CorridorLoad {
            muzzle_velocity_mps: 823.0,
            bc: 0.475,
            drag_model: DragModel::G1,
            mass_kg: 0.010_886,
            diameter_m: 0.007_823,
            bullet_length_m: 0.031,
            zero_distance_m: 91.44,
            sight_height_m: 0.0508,
            temperature_c: t,
            pressure_hpa: p,
            humidity_pct: 50.0,
            altitude_m: alt,
            station_presence: presence,
        },
    };
    solve_robust_hold(&request).unwrap().rows[2].elevation_min_mil
}

#[test]
fn corridor_presence_field() {
    // captured at ab4eb0c
    assert_eq!(
        corridor_elevation_mil(20.0, 1013.0, 150.0, NEITHER).to_bits(),
        0x40112220dde83d2a
    );
    let both = corridor_elevation_mil(20.0, 1013.0, 150.0, BOTH);
    assert_strictly_between(
        corridor_elevation_mil(20.0, 1012.0, 150.0, NEITHER),
        both,
        corridor_elevation_mil(20.0, 1014.0, 150.0, NEITHER),
        "pressure",
    );
    let both = corridor_elevation_mil(15.0, 850.0, 1500.0, BOTH);
    assert_strictly_between(
        corridor_elevation_mil(14.8, 850.0, 1500.0, NEITHER),
        both,
        corridor_elevation_mil(15.2, 850.0, 1500.0, NEITHER),
        "temperature",
    );
}
