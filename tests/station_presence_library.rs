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

// ---------------------------------------------------------------------------------------
// truing: correct_chrono_velocity_fps, calculate_true_velocity_local,
// run_multi_observation_truing_core, solve_for_dsf, solve_wind_truing
// ---------------------------------------------------------------------------------------

use ballistics_engine::truing::{
    calculate_true_velocity_local, calculate_true_velocity_local_with_station_presence,
    correct_chrono_velocity_fps, correct_chrono_velocity_fps_with_station_presence,
    run_multi_observation_truing_core, run_multi_observation_truing_core_with_station_presence,
    DragModelArg, DropUnit, TruingModelInputsV1, TruingObservation, TruingTwist,
};
use ballistics_engine::truing_dsf::{solve_for_dsf, DsfSolveInputs};
use ballistics_engine::truing_wind::{
    modeled_miss_right_m, solve_wind_truing, WindObservation, WindTruingRequest,
};
use ballistics_engine::StationReading;

/// 150 m and 1500 m in the truing core's feet.
const ALT150_FT: f64 = 150.0 / 0.3048;
const ALT1500_FT: f64 = 1500.0 / 0.3048;
const HPA_PER_INHG: f64 = 33.8639;
/// 15.0 °C and its out-of-band neighbours 14.8 / 15.2 °C, in °F.
const F15: f64 = 59.0;
const F14_8: f64 = 14.8 * 9.0 / 5.0 + 32.0;
const F15_2: f64 = 15.2 * 9.0 / 5.0 + 32.0;

fn inhg(hpa: f64) -> f64 {
    hpa / HPA_PER_INHG
}

fn chrono(t: f64, p_inhg: f64, alt_ft: f64, presence: StationAtmospherePresence) -> f64 {
    correct_chrono_velocity_fps_with_station_presence(
        2650.0,
        3.0,
        0.475,
        DragModelArg::G1,
        168.0,
        0.308,
        t,
        p_inhg,
        50.0,
        alt_ft,
        &None,
        presence,
    )
    .unwrap()
    .muzzle_velocity_fps
}

#[test]
fn chrono_correction_twin() {
    let legacy = correct_chrono_velocity_fps(
        2650.0,
        3.0,
        0.475,
        DragModelArg::G1,
        168.0,
        0.308,
        68.0,
        inhg(1013.0),
        50.0,
        ALT150_FT,
        &None,
    )
    .unwrap();
    // captured at ab4eb0c
    assert_eq!(legacy.muzzle_velocity_fps.to_bits(), 0x40a4c02ed353e514);
    assert_eq!(
        chrono(68.0, inhg(1013.0), ALT150_FT, NEITHER).to_bits(),
        0x40a4c02ed353e514
    );
    assert_strictly_between(
        chrono(68.0, inhg(1012.0), ALT150_FT, NEITHER),
        chrono(68.0, inhg(1013.0), ALT150_FT, BOTH),
        chrono(68.0, inhg(1014.0), ALT150_FT, NEITHER),
        "pressure",
    );
    assert_strictly_between(
        chrono(F14_8, inhg(850.0), ALT1500_FT, NEITHER),
        chrono(F15, inhg(850.0), ALT1500_FT, BOTH),
        chrono(F15_2, inhg(850.0), ALT1500_FT, NEITHER),
        "temperature",
    );
}

fn true_velocity(
    t: f64,
    p_inhg: f64,
    alt_ft: f64,
    presence: StationAtmospherePresence,
) -> ballistics_engine::truing::TrueVelocityLocalResult {
    calculate_true_velocity_local_with_station_presence(
        3.6,
        500.0,
        0.475,
        DragModelArg::G1,
        168.0,
        0.308,
        100.0,
        1.5,
        t,
        p_inhg,
        50.0,
        alt_ft,
        &None,
        presence,
    )
    .unwrap()
}

#[test]
fn true_velocity_twin() {
    let legacy = calculate_true_velocity_local(
        3.6,
        500.0,
        0.475,
        DragModelArg::G1,
        168.0,
        0.308,
        100.0,
        1.5,
        68.0,
        inhg(1013.0),
        50.0,
        ALT150_FT,
        &None,
    )
    .unwrap();
    let neither = true_velocity(68.0, inhg(1013.0), ALT150_FT, NEITHER);
    for r in [&legacy, &neither] {
        // captured at ab4eb0c
        assert_eq!(r.effective_velocity_fps.to_bits(), 0x40a43efbc0000000);
        assert_eq!(r.calculated_drop_mil.to_bits(), 0x400ca2bc89427843);
    }
    // The bisection stops at the first midpoint within 0.01 mil, so the fitted velocity is
    // quantized; the supplied air must still move it off the ICAO-air answer.
    let both = true_velocity(68.0, inhg(1013.0), ALT150_FT, BOTH);
    assert_ne!(
        both.effective_velocity_fps.to_bits(),
        neither.effective_velocity_fps.to_bits()
    );
    // 1013 hPa is denser than the ICAO 995 hPa it used to be swapped for: the same measured
    // drop needs a faster bullet.
    assert!(both.effective_velocity_fps > neither.effective_velocity_fps);
}

fn multi(t: f64, p_inhg: f64, alt_ft: f64, presence: StationAtmospherePresence) -> f64 {
    let observations = [
        TruingObservation {
            range_yd: 500.0,
            drop: 3.6,
        },
        TruingObservation {
            range_yd: 700.0,
            drop: 6.6,
        },
    ];
    run_multi_observation_truing_core_with_station_presence(
        &observations,
        DropUnit::Mil,
        0.475,
        DragModelArg::G1,
        168.0,
        0.308,
        100.0,
        1.5,
        t,
        p_inhg,
        50.0,
        alt_ft,
        &None,
        presence,
    )
    .unwrap()
    .fitted_mv_fps
}

#[test]
fn multi_observation_truing_twin() {
    let observations = [
        TruingObservation {
            range_yd: 500.0,
            drop: 3.6,
        },
        TruingObservation {
            range_yd: 700.0,
            drop: 6.6,
        },
    ];
    let legacy = run_multi_observation_truing_core(
        &observations,
        DropUnit::Mil,
        0.475,
        DragModelArg::G1,
        168.0,
        0.308,
        100.0,
        1.5,
        68.0,
        inhg(1013.0),
        50.0,
        ALT150_FT,
        &None,
    )
    .unwrap();
    // captured at ab4eb0c
    assert_eq!(legacy.fitted_mv_fps.to_bits(), 0x40a5268cde4e141c);
    assert_eq!(legacy.rms.to_bits(), 0x3f1d732fcd3ff3eb);
    assert_eq!(
        multi(68.0, inhg(1013.0), ALT150_FT, NEITHER).to_bits(),
        0x40a5268cde4e141c
    );
    assert_strictly_between(
        multi(68.0, inhg(1012.0), ALT150_FT, NEITHER),
        multi(68.0, inhg(1013.0), ALT150_FT, BOTH),
        multi(68.0, inhg(1014.0), ALT150_FT, NEITHER),
        "pressure",
    );
    assert_strictly_between(
        multi(F14_8, inhg(850.0), ALT1500_FT, NEITHER),
        multi(F15, inhg(850.0), ALT1500_FT, BOTH),
        multi(F15_2, inhg(850.0), ALT1500_FT, NEITHER),
        "temperature",
    );
}

fn truing_model(
    temperature_f: StationReading,
    pressure_inhg: StationReading,
    altitude_ft: f64,
) -> TruingModelInputsV1 {
    TruingModelInputsV1 {
        muzzle_velocity_fps: 2700.0,
        ballistic_coefficient: 0.243,
        drag_model: DragModelArg::G7,
        mass_gr: 168.0,
        diameter_in: 0.308,
        zero_distance_yd: 100.0,
        sight_height_in: 1.5,
        temperature_f,
        pressure_inhg,
        humidity_pct: 50.0,
        altitude_ft,
    }
}

fn dsf_impact(t: StationReading, p: StationReading, alt_ft: f64) -> f64 {
    solve_for_dsf(&DsfSolveInputs::from(&truing_model(t, p, alt_ft)), 800.0)
        .unwrap()
        .impact_velocity
}

#[test]
fn dsf_solve_reads_presence_from_the_readings() {
    use StationReading as R;
    // captured at ab4eb0c, when these were bare f64s (the legacy rule)
    assert_eq!(
        dsf_impact(R::omitted(68.0), R::omitted(inhg(1013.0)), ALT150_FT).to_bits(),
        0x40819eba8b0a92de
    );
    assert_strictly_between(
        dsf_impact(R::supplied(68.0), R::supplied(inhg(1012.0)), ALT150_FT),
        dsf_impact(R::supplied(68.0), R::supplied(inhg(1013.0)), ALT150_FT),
        dsf_impact(R::supplied(68.0), R::supplied(inhg(1014.0)), ALT150_FT),
        "pressure",
    );
    assert_strictly_between(
        dsf_impact(R::supplied(F14_8), R::supplied(inhg(850.0)), ALT1500_FT),
        dsf_impact(R::supplied(F15), R::supplied(inhg(850.0)), ALT1500_FT),
        dsf_impact(R::supplied(F15_2), R::supplied(inhg(850.0)), ALT1500_FT),
        "temperature",
    );
}

fn wind_request(t: StationReading, p: StationReading, alt_ft: f64) -> WindTruingRequest {
    WindTruingRequest {
        observations: vec![WindObservation {
            range_m: 457.2,
            miss_right_m: 0.15,
            sigma_m: None,
        }],
        muzzle_velocity_fps: 2700.0,
        bc: 0.475,
        drag_model: DragModelArg::G7,
        mass_gr: 168.0,
        diameter_in: 0.308,
        zero_distance_yd: 100.0,
        sight_height_in: 2.0,
        temperature_f: t,
        pressure_inhg: p,
        humidity_pct: 50.0,
        altitude_ft: alt_ft,
        twist: TruingTwist {
            rate_in: 11.0,
            right_hand: true,
        },
        earth: None,
        called_crosswind_mph: None,
    }
}

#[test]
fn wind_truing_reads_presence_from_the_readings() {
    use StationReading as R;
    let mean = |t, p, alt| {
        solve_wind_truing(&wind_request(t, p, alt))
            .unwrap()
            .mean_crosswind_mph
    };
    let miss = |t, p, alt| modeled_miss_right_m(&wind_request(t, p, alt), 10.0, 457.2).unwrap();
    // captured at ab4eb0c, when these were bare f64s (the legacy rule)
    assert_eq!(
        mean(R::omitted(68.0), R::omitted(inhg(1013.0)), ALT150_FT).to_bits(),
        0x4011af8332a75fb5
    );
    assert_eq!(
        miss(R::omitted(68.0), R::omitted(inhg(1013.0)), ALT150_FT).to_bits(),
        0x3fd20d53e5fa7542
    );
    assert_strictly_between(
        miss(R::supplied(68.0), R::supplied(inhg(1012.0)), ALT150_FT),
        miss(R::supplied(68.0), R::supplied(inhg(1013.0)), ALT150_FT),
        miss(R::supplied(68.0), R::supplied(inhg(1014.0)), ALT150_FT),
        "pressure",
    );
    assert_strictly_between(
        miss(R::supplied(F14_8), R::supplied(inhg(850.0)), ALT1500_FT),
        miss(R::supplied(F15), R::supplied(inhg(850.0)), ALT1500_FT),
        miss(R::supplied(F15_2), R::supplied(inhg(850.0)), ALT1500_FT),
        "temperature",
    );
    assert_ne!(
        mean(R::supplied(68.0), R::supplied(inhg(1013.0)), ALT150_FT).to_bits(),
        0x4011af8332a75fb5
    );
}

#[test]
fn truing_dtos_accept_omitted_readings_on_the_wire() {
    let mut json = serde_json::to_value(truing_model(
        StationReading::supplied(59.0),
        StationReading::supplied(29.92),
        5000.0,
    ))
    .unwrap();
    assert_eq!(json["temperature_f"], serde_json::json!(59.0));
    assert_eq!(json["pressure_inhg"], serde_json::json!(29.92));
    let map = json.as_object_mut().unwrap();
    map.remove("temperature_f");
    map.remove("pressure_inhg");
    let omitted: TruingModelInputsV1 = serde_json::from_value(json.clone()).unwrap();
    assert_eq!(omitted.temperature_f, StationReading::omitted(59.0));
    assert_eq!(omitted.pressure_inhg, StationReading::omitted(29.92));
    // An omitted reading is not written back out.
    let round_trip = serde_json::to_value(omitted).unwrap();
    assert!(round_trip.get("temperature_f").is_none());
    assert!(round_trip.get("pressure_inhg").is_none());
    // deny_unknown_fields still holds. A `null` reading is omitted, as on `card.*`; a
    // string is still refused.
    let mut unknown = json.clone();
    unknown["pressure_hpa"] = serde_json::json!(1013.25);
    assert!(serde_json::from_value::<TruingModelInputsV1>(unknown).is_err());
    let mut null = json.clone();
    null["temperature_f"] = serde_json::Value::Null;
    null["pressure_inhg"] = serde_json::Value::Null;
    let from_null: TruingModelInputsV1 = serde_json::from_value(null).unwrap();
    assert_eq!(from_null.temperature_f, StationReading::omitted(59.0));
    assert_eq!(from_null.pressure_inhg, StationReading::omitted(29.92));
    let mut string = json;
    string["pressure_inhg"] = serde_json::json!("29.92");
    assert!(serde_json::from_value::<TruingModelInputsV1>(string).is_err());

    let mut wind = serde_json::to_value(wind_request(
        StationReading::supplied(59.0),
        StationReading::supplied(29.92),
        5000.0,
    ))
    .unwrap();
    let map = wind.as_object_mut().unwrap();
    map.remove("temperature_f");
    map.remove("pressure_inhg");
    let omitted: WindTruingRequest = serde_json::from_value(wind.clone()).unwrap();
    assert_eq!(omitted.temperature_f, StationReading::omitted(59.0));
    assert_eq!(omitted.pressure_inhg, StationReading::omitted(29.92));
    let mut null = wind.clone();
    null["temperature_f"] = serde_json::Value::Null;
    null["pressure_inhg"] = serde_json::Value::Null;
    let from_null: WindTruingRequest = serde_json::from_value(null).unwrap();
    assert_eq!(from_null.temperature_f, StationReading::omitted(59.0));
    assert_eq!(from_null.pressure_inhg, StationReading::omitted(29.92));
    wind["bogus"] = serde_json::json!(1);
    assert!(serde_json::from_value::<WindTruingRequest>(wind).is_err());
}
