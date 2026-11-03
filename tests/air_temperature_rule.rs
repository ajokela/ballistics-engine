//! One air-temperature rule on every surface.
//!
//! - **The rule.** A station temperature flies when it is finite and above absolute zero and the
//!   air it describes can exist (the impossible-air rule). Nothing narrower: no surface refuses a
//!   temperature the solve flies, and `atmosphere.density_altitude` answers for exactly the air
//!   a solve with the same `atmosphere` object flies (outside 173.15-373.15 K the humidity terms
//!   are extrapolated, never refused).
//! - **Dry air has no vapour at any temperature.** At 0 % the water-vapour mole fraction is
//!   exactly 0, so dry air above water's critical point (647.096 K) flies; only humid air there
//!   cannot exist.
//! - **Air that cannot exist aloft.** A station whose lapse takes the air to absolute zero
//!   within a flight is refused at the step that reaches it, naming the temperature, on every
//!   kernel: never a NaN trajectory, never an unrelated error. A flight that stays below that
//!   height flies as it always did.

use ballistics_engine::atmosphere::calculate_atmosphere;
use ballistics_engine::fast_trajectory::{
    fast_integrate, fast_integrate_with_segments, FastIntegrationParams,
};
use ballistics_engine::monte_carlo::solve_trajectory_for_monte_carlo_with_station_presence;
use ballistics_engine::solve_json::decode_solve_request_v1;
use ballistics_engine::trajectory_integration::{
    try_integrate_trajectory, KernelInputError, TrajectoryParams,
};
use ballistics_engine::wind::WindSock;
use ballistics_engine::{
    AtmosphereErrorReason, AtmosphericConditions, BallisticInputs, DragModel, KernelAtmosphere,
    RelativeHumidity, StationAir, StationAtmospherePresence, TrajectorySolver, WindConditions,
};
use serde_json::{json, Value};

/// The words every refusal of air that cannot exist aloft carries.
const ALOFT: &str = "the air aloft would fall below absolute zero";

fn solve_request(atmosphere: Value, muzzle_angle_rad: Option<f64>, method: &str) -> Value {
    let shot = match muzzle_angle_rad {
        Some(angle) => json!({"max_range_m": 20000.0, "muzzle_angle_rad": angle}),
        None => json!({"max_range_m": 300.0, "zero_distance_m": 100.0}),
    };
    json!({
        "schema_version": 1,
        "projectile": {"mass_kg": 0.0486, "diameter_m": 0.0127, "drag_model": "G1",
                       "ballistic_coefficient": 1.05},
        "rifle": {"muzzle_velocity_mps": 823.0},
        "shot": shot,
        "atmosphere": atmosphere,
        "wind": {"speed_mps": 0.0, "direction_from_rad": 0.0},
        "solver": {"method": method, "time_step_s": 0.001},
        "effects": {"magnus": false, "coriolis": false, "enhanced_spin_drift": false},
        "sampling": {"interval_m": 500.0}
    })
}

/// A solve of `request`, as the wire reports it: `Ok(summary)` or `Err(error)` (code, path and
/// message).
fn solve(request: &Value) -> Result<Value, Value> {
    let decoded = decode_solve_request_v1(&request.to_string())
        .map_err(|envelope| serde_json::to_value(envelope).unwrap()["error"].clone())?;
    ballistics_engine::solve_v1(decoded)
        .map(|success| serde_json::to_value(success).unwrap()["summary"].clone())
        .map_err(|envelope| serde_json::to_value(envelope).unwrap()["error"].clone())
}

/// `atmosphere.density_altitude` of `atmosphere`, as the wire reports it.
fn density_altitude(atmosphere: &Value) -> Result<Value, Value> {
    use ballistics_engine::atmosphere_service::{
        decode_density_altitude_request_v1, density_altitude_v1,
    };
    decode_density_altitude_request_v1(&json!({ "atmosphere": atmosphere }))
        .and_then(|request| density_altitude_v1(&request))
        .map(|response| serde_json::to_value(response).unwrap())
        .map_err(|envelope| serde_json::to_value(envelope).unwrap()["error"].clone())
}

// --- dry air above water's critical point --------------------------------------------------------

/// 400 °C is above water's critical point (647.096 K), where no saturation pressure is defined.
/// Dry air holds no vapour there either, so it is ordinary hot air: CIPM-2007 within 0.2 % of
/// the ideal gas, and the dry sound speed to the bit. Humid air at that temperature cannot exist.
#[test]
fn dry_air_above_waters_critical_point_flies_and_humid_air_there_is_refused() {
    let temperature_c = 400.0;
    let temperature_k = temperature_c + 273.15;
    let (density, sound) = calculate_atmosphere(0.0, Some(temperature_c), Some(1013.25), 0.0);
    let ideal = 101_325.0 / (287.0531 * temperature_k);
    assert!(
        (density / ideal - 1.0).abs() < 2e-3,
        "dry air at 400 °C: {density} kg/m³ against the ideal gas's {ideal}"
    );
    assert_eq!(
        sound.to_bits(),
        (1.4_f64 * 287.0531 * temperature_k).sqrt().to_bits(),
        "dry sound speed at 400 °C: {sound}"
    );

    let station = StationAir::new(0.0, temperature_c, 1013.25, RelativeHumidity::DRY)
        .expect("dry air at 400 °C exists");
    assert_eq!(station.density_kg_m3().to_bits(), density.to_bits());
    assert_eq!(station.speed_of_sound_mps().to_bits(), sound.to_bits());

    let humid = RelativeHumidity::try_percent(1.0).unwrap();
    let error = StationAir::new(0.0, temperature_c, 1013.25, humid)
        .expect_err("humid air above water's critical point cannot exist");
    assert_eq!(error.reason(), AtmosphereErrorReason::ImpossibleAir);
    assert_eq!(error.field(), "temperature_c");

    // The wire: dry air at 700 K flies; 1 % there is refused at the temperature.
    let dry = solve(&solve_request(
        json!({"temperature_k": 700.0, "pressure_pa": 101325.0, "relative_humidity": 0.0}),
        None,
        "rk45",
    ));
    assert!(dry.is_ok(), "dry air at 700 K flies: {dry:?}");
    let humid = solve(&solve_request(
        json!({"temperature_k": 700.0, "pressure_pa": 101325.0, "relative_humidity": 0.01}),
        None,
        "rk45",
    ))
    .expect_err("1 % at 700 K cannot exist");
    assert_eq!(humid["code"], "invalid_value", "{humid}");
    assert_eq!(humid["path"], "$.atmosphere.temperature_k", "{humid}");
}

// --- one temperature rule on the wire ------------------------------------------------------------

/// `atmosphere.density_altitude` accepts every temperature a solve flies: 110 °C (an app's
/// meter or a typed 230 °F) and 160 K under 500 hPa used to be refused as "outside 173.15-373.15
/// K"; a solve flew both. A temperature either command refuses, the other refuses at the same
/// path with the same words.
#[test]
fn density_altitude_answers_for_every_temperature_a_solve_flies() {
    for atmosphere in [
        json!({"temperature_k": 383.15, "pressure_pa": 101325.0, "relative_humidity": 0.5}),
        json!({"temperature_k": 160.0, "pressure_pa": 50000.0, "relative_humidity": 0.5}),
        json!({"temperature_k": 500.0, "pressure_pa": 101325.0, "relative_humidity": 0.0}),
    ] {
        let flown = solve(&solve_request(atmosphere.clone(), None, "rk45"));
        assert!(flown.is_ok(), "{atmosphere}: the solve flies it: {flown:?}");
        let answered = density_altitude(&atmosphere);
        let answered = answered.unwrap_or_else(|error| {
            panic!("{atmosphere}: a solve flies this air, so the command answers: {error}")
        });
        assert!(answered["air_density_kg_m3"].as_f64().unwrap() > 0.0);
    }
    for atmosphere in [
        json!({"temperature_k": 1e-300, "pressure_pa": 101325.0, "relative_humidity": 0.5}),
        json!({"temperature_k": 400.0, "pressure_pa": 101325.0, "relative_humidity": 0.9}),
        json!({"temperature_k": 700.0, "pressure_pa": 101325.0, "relative_humidity": 0.5}),
    ] {
        let flown = solve(&solve_request(atmosphere.clone(), None, "rk45"))
            .expect_err("the solve refuses this air");
        let answered = density_altitude(&atmosphere).expect_err("the command refuses this air");
        assert_eq!(flown["path"], "$.atmosphere.temperature_k", "{flown}");
        assert_eq!(
            (&answered["code"], &answered["path"], &answered["message"]),
            (&flown["code"], &flown["path"], &flown["message"]),
            "{atmosphere}"
        );
    }
}

// --- air that cannot exist aloft -----------------------------------------------------------------

/// A 10 K station under 10 hPa: the standard lapse takes its air to 0 K about 1.5 km up.
const COLD_ATMOSPHERE_K: f64 = 10.0;
const COLD_PRESSURE_PA: f64 = 1000.0;

fn cold_atmosphere() -> Value {
    json!({"temperature_k": COLD_ATMOSPHERE_K, "pressure_pa": COLD_PRESSURE_PA,
           "relative_humidity": 0.0})
}

/// A lob from that station climbs past the height where its air would fall below absolute zero:
/// every solver method refuses it at the temperature, in the wire's words.
#[test]
fn a_lob_into_air_below_absolute_zero_is_refused_at_the_temperature() {
    for method in ["rk45", "rk4", "euler"] {
        let error = solve(&solve_request(cold_atmosphere(), Some(0.7), method))
            .expect_err("the lob reaches air below absolute zero");
        assert_eq!(error["code"], "invalid_value", "{method}: {error}");
        assert_eq!(
            error["path"], "$.atmosphere.temperature_k",
            "{method}: {error}"
        );
        let message = error["message"].as_str().unwrap_or_default();
        assert!(
            message.starts_with("temperature_k = 10 K: ") && message.contains(ALOFT),
            "{method}: {message}"
        );
    }
}

/// The same station flies a flat shot, which never climbs that high, exactly as before.
#[test]
fn a_flight_that_stays_below_that_height_flies() {
    for method in ["rk45", "rk4", "euler"] {
        let summary = solve(&solve_request(cold_atmosphere(), Some(0.02), method))
            .unwrap_or_else(|error| panic!("{method}: the flat shot flies: {error}"));
        assert!(
            summary["maximum_height_m"].as_f64().unwrap() < 1000.0,
            "{summary}"
        );
    }
}

fn cold_station() -> StationAir {
    StationAir::new(
        0.0,
        COLD_ATMOSPHERE_K - 273.15,
        COLD_PRESSURE_PA / 100.0,
        RelativeHumidity::DRY,
    )
    .expect("a 10 K station flies")
}

fn lob_inputs() -> BallisticInputs {
    BallisticInputs {
        bc_value: 1.05,
        bc_type: DragModel::G1,
        bullet_mass: 0.0486,
        bullet_diameter: 0.0127,
        muzzle_velocity: 823.0,
        muzzle_angle: 0.7,
        target_distance: 20_000.0,
        temperature: COLD_ATMOSPHERE_K - 273.15,
        pressure: COLD_PRESSURE_PA / 100.0,
        humidity_percent: 0.0,
        ..BallisticInputs::default()
    }
}

/// `TrajectorySolver` names the temperature (it used to blame "an absurd bore/muzzle height").
#[test]
fn the_solver_names_the_temperature() {
    let inputs = lob_inputs();
    let mut solver = TrajectorySolver::new_with_station_presence(
        inputs.clone(),
        WindConditions::default(),
        AtmosphericConditions {
            temperature: inputs.temperature,
            pressure: inputs.pressure,
            humidity: 0.0,
            altitude: 0.0,
        },
        StationAtmospherePresence::BOTH,
    );
    solver.set_max_range(20_000.0);
    let error = solver.solve().expect_err("the lob is refused");
    let message = error.to_string();
    assert!(message.contains(ALOFT), "{message}");
    assert!(!message.contains("bore/muzzle height"), "{message}");
}

/// `fast_integrate` stops at the step that reaches it: a failed solution, never a NaN state.
#[test]
fn fast_integrate_fails_rather_than_flying_nan() {
    let v = 823.0_f64;
    let params = FastIntegrationParams {
        horiz: 20_000.0,
        vert: 0.0,
        initial_state: [0.0, 0.0, 0.0, v * 0.7_f64.cos(), v * 0.7_f64.sin(), 0.0],
        t_span: (0.0, 60.0),
        atmosphere: KernelAtmosphere::Station(cold_station()),
        atmo_sock: None,
    };
    let solution = fast_integrate(&lob_inputs(), &WindSock::new(vec![]), params);
    assert!(
        !solution.success,
        "a flight into air below absolute zero is no solution"
    );
    let non_finite = solution
        .y
        .iter()
        .flatten()
        .filter(|x| !x.is_finite())
        .count();
    assert_eq!(non_finite, 0, "no NaN state is returned");
}

fn lob_params() -> TrajectoryParams {
    TrajectoryParams {
        mass_kg: 0.0486,
        bc: 1.05,
        drag_model: DragModel::G1,
        wind_segments: vec![],
        atmosphere: KernelAtmosphere::Station(cold_station()),
        omega_vector: None,
        enable_spin_drift: false,
        enable_magnus: false,
        enable_coriolis: false,
        target_distance_m: 20_000.0,
        enable_wind_shear: false,
        wind_shear_model: "none".into(),
        shooter_altitude_m: 0.0,
        is_twist_right: true,
        shooting_angle: 0.0,
        bullet_diameter: 0.0127,
        bullet_length: 0.0,
        twist_rate: 15.0,
        custom_drag_table: None,
        cd_scale: 1.0,
        custom_drag_reference: false,
        bc_segments: None,
        use_bc_segments: false,
        bc_segments_data: None,
        ground_threshold: -1000.0,
        atmo_sock: None,
    }
}

/// The checked kernel refuses with a typed atmosphere error: RK45 used to stop silently and
/// report the climb as a trajectory, RK4 to return 60 s of NaN.
#[test]
fn the_checked_kernel_refuses_with_a_typed_error() {
    let v = 823.0_f64;
    let initial_state = [0.0, 0.0, 0.0, v * 0.7_f64.cos(), v * 0.7_f64.sin(), 0.0];
    for method in ["RK45", "RK4"] {
        match try_integrate_trajectory(initial_state, (0.0, 60.0), lob_params(), method, 1e-6, 0.01)
        {
            Err(KernelInputError::Atmosphere(error)) => {
                assert_eq!(error.field(), "temperature_c", "{method}: {error}");
                assert!(error.to_string().contains(ALOFT), "{method}: {error}");
            }
            other => panic!("{method}: refused as air that cannot exist aloft, got {other:?}"),
        }
    }
}

/// The lean Monte Carlo kernel refuses the run instead of reporting a NaN impact.
#[test]
fn the_monte_carlo_kernel_refuses_rather_than_reporting_nan() {
    let error = solve_trajectory_for_monte_carlo_with_station_presence(
        &lob_inputs(),
        StationAtmospherePresence::BOTH,
    )
    .expect_err("the lob is refused");
    assert!(error.contains(ALOFT), "{error}");
}

/// A zero search whose steeper trials climb into air that cannot exist names the temperature. It
/// used to try them, discard each refusal and report "target beyond effective range even at
/// maximum angle". A target 1000 m up at 2000 m from the 10 K station is one: every bracket from
/// 0.4 rad up climbs past 1539 m.
#[test]
fn a_zero_search_that_climbs_into_air_below_absolute_zero_names_the_temperature() {
    let inputs = BallisticInputs {
        muzzle_angle: 0.0,
        target_distance: 2000.0,
        ..lob_inputs()
    };
    let error = ballistics_engine::calculate_zero_angle_with_resolved_conditions(
        inputs,
        2000.0,
        1000.0,
        WindConditions::default(),
        AtmosphericConditions {
            temperature: COLD_ATMOSPHERE_K - 273.15,
            pressure: COLD_PRESSURE_PA / 100.0,
            humidity: 0.0,
            altitude: 0.0,
        },
    )
    .expect_err("no zero exists in air that exists");
    let refusal = error
        .atmosphere_error()
        .unwrap_or_else(|| panic!("the refusal of the air, not `{error}`"));
    assert_eq!(
        refusal.reason(),
        AtmosphereErrorReason::AirAloftBelowAbsoluteZero
    );
    assert!(error.to_string().contains(ALOFT), "{error}");
}

/// Any step that fails to a non-finite state ends a fast solution as a failure, on both fast
/// kernels and for any cause: 1e200 m/s overflows the drag at standard air.
/// `fast_integrate_with_segments` used to report it as `success = true` whenever the cause was
/// not air below absolute zero; `fast_integrate` already failed it.
#[test]
fn both_fast_kernels_fail_a_flight_whose_step_goes_non_finite() {
    let station = StationAir::new(0.0, 15.0, 1013.25, RelativeHumidity::DRY).unwrap();
    let v = 1e200;
    let inputs = BallisticInputs {
        bc_value: 0.5,
        bc_type: DragModel::G1,
        bullet_mass: 0.01,
        bullet_diameter: 0.00782,
        muzzle_velocity: v,
        target_distance: 1000.0,
        ..BallisticInputs::default()
    };
    let params = || FastIntegrationParams {
        horiz: 1000.0,
        vert: 0.0,
        initial_state: [0.0, 0.0, 0.0, v, 0.0, 0.0],
        t_span: (0.0, 10.0),
        atmosphere: KernelAtmosphere::Station(station),
        atmo_sock: None,
    };
    for (name, solution) in [
        (
            "fast_integrate",
            fast_integrate(&inputs, &WindSock::new(vec![]), params()),
        ),
        (
            "fast_integrate_with_segments",
            fast_integrate_with_segments(&inputs, vec![], params()),
        ),
    ] {
        assert!(
            !solution.success,
            "{name}: a non-finite step is no solution"
        );
        assert!(
            solution.y.iter().flatten().all(|x| x.is_finite()),
            "{name}: no NaN state is returned"
        );
    }
}

// --- one rule on every surface, at the temperatures a unit mistake or a meter gives --------------

/// The air every surface must treat alike: °C sent as kelvin and meter readings far outside
/// -100..100 °C (100 K, 150 K, the humid-air formulas' 173.15 K edge, 383.15 K), and air above
/// water's critical point (647 K, 700 K), at sea level and under 500 hPa, dry and at 50 %:
/// `(temperature_k, pressure_pa, relative_humidity)`.
fn named_atmospheres() -> Vec<(f64, f64, f64)> {
    let mut atmospheres = Vec::new();
    for pressure_pa in [101_325.0, 50_000.0] {
        for relative_humidity in [0.0, 0.5] {
            for temperature_k in [100.0, 150.0, 173.15, 383.15, 647.0, 700.0] {
                atmospheres.push((temperature_k, pressure_pa, relative_humidity));
            }
        }
    }
    atmospheres
}

/// How many of [`named_atmospheres`] can exist, and so fly: all six temperatures dry; at 50 %,
/// none above the critical point, and at 500 hPa not 383.15 K either (its saturation pressure
/// is above half the air's).
const NAMED_ATMOSPHERES_THAT_FLY: usize = 19;

/// `card.range_table` of the same air through the bridge (an imperial card, its temperature in
/// °C): `Ok` when it built.
#[cfg(feature = "bridge")]
fn card(temperature_k: f64, pressure_pa: f64, relative_humidity: f64) -> Result<(), Value> {
    let request = json!({"api_version": 1, "command": "card.range_table", "request": {
        "muzzle_velocity": 2600.0, "ballistic_coefficient": 0.243, "mass": 175.0,
        "diameter": 0.308, "zero_distance": 100.0, "start": 100.0, "end": 300.0, "step": 100.0,
        "temperature": temperature_k - 273.15, "temperature_unit": "celsius",
        "pressure": pressure_pa / 100.0 / 33.8639, "humidity": relative_humidity * 100.0
    }});
    let out: Value = serde_json::from_str(&ballistics_engine::bridge::bridge_call(
        &request.to_string(),
    ))
    .unwrap();
    if out["ok"] == true {
        Ok(())
    } else {
        Err(out["error"].clone())
    }
}

/// Each named atmosphere is flown, or refused, alike by the library's station, `solve`,
/// `atmosphere.density_altitude` and `card.range_table`: no surface refuses air another flies.
/// `atmosphere.density_altitude` answers with finite figures and the solver's own density. It
/// used to refuse 100 K, 150 K and 173.15 K at sea level, and 647 K and 700 K dry at both
/// pressures, all of which `solve` flew, because their density altitudes lay outside the
/// troposphere. What one surface refuses, `solve` and `atmosphere.density_altitude` refuse with
/// the same code, path and words.
#[test]
fn every_surface_flies_or_refuses_each_named_atmosphere_alike() {
    let mut flown = 0;
    for (temperature_k, pressure_pa, relative_humidity) in named_atmospheres() {
        let label = format!("{temperature_k} K, {pressure_pa} Pa, RH {relative_humidity}");
        let atmosphere = json!({"temperature_k": temperature_k, "pressure_pa": pressure_pa,
                                "relative_humidity": relative_humidity});
        let station = StationAir::new(
            0.0,
            temperature_k - 273.15,
            pressure_pa / 100.0,
            RelativeHumidity::try_fraction(relative_humidity).unwrap(),
        );
        let solved = solve(&solve_request(atmosphere.clone(), None, "rk45"));
        let answered = density_altitude(&atmosphere);
        assert_eq!(
            station.is_ok(),
            solved.is_ok(),
            "{label}: the library's station and solve disagree: {solved:?}"
        );
        match (&solved, &answered) {
            (Ok(_), Ok(answer)) => {
                flown += 1;
                let station = station.unwrap();
                for figure in [
                    &answer["pressure_altitude"],
                    &answer["density_altitude"]["faa_rule"],
                    &answer["density_altitude"]["density_matched"],
                ] {
                    for unit in ["m", "ft"] {
                        assert!(
                            figure[unit].as_f64().is_some_and(f64::is_finite),
                            "{label}: {answer}"
                        );
                    }
                }
                assert_eq!(
                    answer["air_density_kg_m3"].as_f64().map(f64::to_bits),
                    Some(station.density_kg_m3().to_bits()),
                    "{label}: {answer}"
                );
            }
            (Err(refused), Err(declined)) => assert_eq!(
                (&declined["code"], &declined["path"], &declined["message"]),
                (&refused["code"], &refused["path"], &refused["message"]),
                "{label}"
            ),
            _ => panic!("{label}: solve {solved:?}, but atmosphere.density_altitude {answered:?}"),
        }
        #[cfg(feature = "bridge")]
        {
            let carded = card(temperature_k, pressure_pa, relative_humidity);
            assert_eq!(
                carded.is_ok(),
                solved.is_ok(),
                "{label}: card.range_table disagrees with solve: {carded:?}"
            );
        }
    }
    assert_eq!(flown, NAMED_ATMOSPHERES_THAT_FLY);
}

// --- density altitude beyond the troposphere -----------------------------------------------------

/// The ICAO standard atmosphere at a geopotential height, m, computed here independently of the
/// engine from its published layers (base height m, lapse K/m), the pressure chained from
/// 101325 Pa and 288.15 K with g0 = 9.80665 m/s² and R = 287.0531 J/(kg·K): `(T K, P Pa)`.
fn standard_air(height_m: f64) -> (f64, f64) {
    const LAYERS: [(f64, f64); 7] = [
        (0.0, -0.0065),
        (11_000.0, 0.0),
        (20_000.0, 0.001),
        (32_000.0, 0.0028),
        (47_000.0, 0.0),
        (51_000.0, -0.0028),
        (71_000.0, -0.002),
    ];
    let (g, r) = (9.80665, 287.0531);
    let (mut temperature, mut pressure) = (288.15, 101_325.0);
    for (index, &(base, lapse)) in LAYERS.iter().enumerate() {
        let top = LAYERS.get(index + 1).map_or(f64::INFINITY, |next| next.0);
        let rise = height_m.min(top) - base;
        let next = temperature + lapse * rise;
        pressure *= if lapse == 0.0 {
            (-g * rise / (r * temperature)).exp()
        } else {
            (next / temperature).powf(-g / (lapse * r))
        };
        temperature = next;
        if height_m <= top {
            break;
        }
    }
    (temperature, pressure)
}

/// `density_matched` is the standard atmosphere's altitude of equal density, so standard air
/// reads its own altitude in every layer up to 84 km, within the few metres the CIPM-2007
/// compressibility puts between its density and the ideal gas's. Above the tropopause too:
/// the troposphere's closed form it used to continue there read standard air at 20 km as 20451
/// m. And it falls steadily as the air thins, across every layer boundary.
#[test]
fn density_matched_reads_standard_air_at_its_own_altitude_in_every_layer() {
    use ballistics_engine::atmosphere::density_matched_altitude_m;
    for height_m in [
        -4000.0, 0.0, 5000.0, 10_999.0, 11_000.0, 11_001.0, 15_000.0, 20_000.0, 25_000.0, 32_000.0,
        40_000.0, 47_000.0, 49_000.0, 51_000.0, 60_000.0, 71_000.0, 80_000.0, 84_000.0,
    ] {
        let (temperature_k, pressure_pa) = standard_air(height_m);
        let altitude = density_matched_altitude_m(temperature_k - 273.15, pressure_pa / 100.0, 0.0);
        assert!(
            (altitude - height_m).abs() < 5.0,
            "standard air at {height_m} m reads {altitude} m"
        );
    }
    let mut last = f64::NEG_INFINITY;
    for step in 0..=600 {
        let pressure_hpa = 2000.0 * 10f64.powf(-f64::from(step) / 100.0);
        let altitude = density_matched_altitude_m(-56.5, pressure_hpa, 0.0);
        assert!(
            altitude.is_finite() && altitude > last,
            "{pressure_hpa} hPa reads {altitude} m, not above {last} m"
        );
        last = altitude;
    }
}

/// Up to 11 km `density_matched` is the troposphere's closed form exactly as 0.49.0 computed it,
/// to the bit: every figure `atmosphere.density_altitude` gave before is unchanged.
#[test]
fn density_matched_is_the_troposphere_closed_form_to_the_bit_up_to_11_km() {
    use ballistics_engine::atmosphere::density_matched_altitude_m;
    let (reference, _) = calculate_atmosphere(0.0, Some(15.0), Some(1013.25), 0.0);
    let exponent = 9.80665 / (287.0531 * 0.0065) - 1.0;
    let mut compared = 0;
    for temperature_c in (-100..=100).step_by(5).map(f64::from) {
        for pressure_hpa in (230..=1100).step_by(10).map(f64::from) {
            for humidity in [0.0, 50.0, 100.0] {
                let (density, _) =
                    calculate_atmosphere(0.0, Some(temperature_c), Some(pressure_hpa), humidity);
                let closed_form =
                    (288.15 / 0.0065) * (1.0 - (density / reference).powf(1.0 / exponent));
                if closed_form.is_nan() || closed_form > 11_000.0 {
                    continue;
                }
                compared += 1;
                assert_eq!(
                    density_matched_altitude_m(temperature_c, pressure_hpa, humidity).to_bits(),
                    closed_form.to_bits(),
                    "{temperature_c} C, {pressure_hpa} hPa, {humidity} %"
                );
            }
        }
    }
    assert!(compared > 3000, "{compared}");
}

// --- the figures the docs quote, from the engine -------------------------------------------------

fn read(relative: &str) -> String {
    std::fs::read_to_string(std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(relative))
        .unwrap_or_else(|error| panic!("{relative}: {error}"))
}

/// `text` with every run of whitespace one space, so a quote may wrap.
fn normalized(text: &str) -> String {
    text.split_whitespace().collect::<Vec<_>>().join(" ")
}

fn assert_quotes(document: &str, quotes: &[String]) {
    let text = normalized(&read(document));
    for quote in quotes {
        assert!(
            text.contains(&normalized(quote)),
            "{document} does not quote the engine: `{quote}`"
        );
    }
}

/// Dry air at 400 °C as the CHANGELOG quotes it, now and in the released 0.49.0 (the 0.49.0
/// dev-dependency is native only).
#[cfg(not(target_arch = "wasm32"))]
#[test]
fn the_changelog_quotes_dry_hot_air_from_the_engine_and_from_0_49_0() {
    let (density, sound) = calculate_atmosphere(0.0, Some(400.0), Some(1013.25), 0.0);
    let (steam_density, steam_sound) = ballistics_engine_049::atmosphere::calculate_atmosphere(
        0.0,
        Some(400.0),
        Some(1013.25),
        0.0,
    );
    assert_quotes(
        "CHANGELOG.md",
        &[
            format!("at 400 °C and 1013.25 hPa, {density:.5} kg/m³"),
            format!("and {sound:.3} m/s, the dry sound speed"),
            format!(
                "pure steam, {steam_density:.5} kg/m³ and {steam_sound:.3} m/s, at any humidity"
            ),
        ],
    );
}

/// Every other figure the CHANGELOG, SOLVE_JSON_V1.md and ATMOSPHERE_DENSITY_ALTITUDE.md quote
/// for the one temperature rule, recomputed: the density altitudes `atmosphere.density_altitude`
/// gives the air it now answers for, wherever they land, and the refusal of air that cannot exist
/// aloft.
#[test]
fn the_docs_quote_the_engine() {
    let answer = |atmosphere: &Value| {
        density_altitude(atmosphere).unwrap_or_else(|error| panic!("{atmosphere}: {error}"))
    };
    let figure = |atmosphere: &Value, name: &str| {
        answer(atmosphere)["density_altitude"][name]["m"]
            .as_f64()
            .unwrap()
    };
    // A figure as the prose quotes it: whole metres, a negative one with a typographic minus.
    let prose = |metres: f64| format!("{metres:.0}").replace('-', "\u{2212}");
    let hot = json!({"temperature_k": 383.15, "pressure_pa": 101325.0, "relative_humidity": 0.5});
    let cold = json!({"temperature_k": 160.0, "pressure_pa": 50000.0, "relative_humidity": 0.5});
    let kelvin_as_celsius = json!({"temperature_k": 15.0});
    let dry_100 =
        json!({"temperature_k": 100.0, "pressure_pa": 101325.0, "relative_humidity": 0.0});
    let dry_700 =
        json!({"temperature_k": 700.0, "pressure_pa": 101325.0, "relative_humidity": 0.0});
    let dry_700_low =
        json!({"temperature_k": 700.0, "pressure_pa": 50000.0, "relative_humidity": 0.0});
    let dry_700_figures = format!(
        "{:.0} m and {:.0} m at sea level and {:.0} m and {:.0} m under 500 hPa",
        figure(&dry_700, "faa_rule"),
        figure(&dry_700, "density_matched"),
        figure(&dry_700_low, "faa_rule"),
        figure(&dry_700_low, "density_matched")
    );
    let (standard_20_km, _) = standard_air(20_000.0);
    let standard_20_km = ballistics_engine::atmosphere::density_matched_altitude_m(
        standard_20_km - 273.15,
        standard_air(20_000.0).1 / 100.0,
        0.0,
    );
    assert_quotes(
        "docs/ATMOSPHERE_DENSITY_ALTITUDE.md",
        &[
            format!(
                "383.15 K at 50 % (an app's 230 °F) reads `faa_rule` {:.0} m and `density_matched` {:.0} m",
                figure(&hot, "faa_rule"),
                figure(&hot, "density_matched")
            ),
            format!(
                "160 K at 50 % under 500 hPa reads {:.0} m and {:.0} m",
                figure(&cold, "faa_rule"),
                figure(&cold, "density_matched")
            ),
            format!(
                "°C sent as `temperature_k` (15 K) reads {} m and {} m, 100 K dry at sea level {} m \
                 and {} m, and dry air at 700 K, which only dry air can be, {dry_700_figures}",
                prose(figure(&kelvin_as_celsius, "faa_rule")),
                prose(figure(&kelvin_as_celsius, "density_matched")),
                prose(figure(&dry_100, "faa_rule")),
                prose(figure(&dry_100, "density_matched")),
            ),
            format!(
                "`faa_rule` {} m and `density_matched` {} m, at an `air_density_kg_m3` of {:.2}",
                prose(figure(&kelvin_as_celsius, "faa_rule")),
                prose(figure(&kelvin_as_celsius, "density_matched")),
                answer(&kelvin_as_celsius)["air_density_kg_m3"]
                    .as_f64()
                    .unwrap()
            ),
            format!("it now reads {standard_20_km:.0} m"),
            // hPa sent as pascals, dry, is answered as the solve flies it (0.50.0): the
            // command no longer refuses a station pressure for its pressure altitude.
            format!(
                "answers for it with a `pressure_altitude` of {:.0} m",
                answer(&json!({"pressure_pa": 1013.25, "relative_humidity": 0.0}))
                    ["pressure_altitude"]["m"]
                    .as_f64()
                    .unwrap()
            ),
        ],
    );
    assert_quotes(
        "CHANGELOG.md",
        &[
            format!(
                "(15 K) reads `faa_rule` {:.0} m and `density_matched` {:.0} m, dry air at 700 K \
                 {:.0} m and {:.0} m",
                figure(&kelvin_as_celsius, "faa_rule"),
                figure(&kelvin_as_celsius, "density_matched"),
                figure(&dry_700, "faa_rule"),
                figure(&dry_700, "density_matched")
            ),
            format!("where it now reads {standard_20_km:.0} m"),
        ],
    );

    let aloft = solve(&solve_request(cold_atmosphere(), Some(0.7), "rk45")).expect_err("refused")
        ["message"]
        .as_str()
        .unwrap()
        .to_string();
    assert_quotes("docs/SOLVE_JSON_V1.md", std::slice::from_ref(&aloft));
    assert_quotes(
        "CHANGELOG.md",
        &[aloft
            .split(": from this temperature")
            .next()
            .unwrap()
            .to_string()],
    );
}

/// What 0.49.0 gave, as the docs quote it, from the released crate (the 0.49.0 dev-dependency is
/// native only): its closed-form troposphere read standard air at 20 km some 450 m high;
/// `fast_integrate_with_segments` reported a flight that overflowed at its first step as a
/// success; and the zero search for a target 1000 m up at 2000 m from the 10 K station blamed
/// the range.
#[cfg(not(target_arch = "wasm32"))]
#[test]
fn the_docs_quote_0_49_0_for_what_changed() {
    use ballistics_engine_049 as released;
    let (temperature_k, pressure_pa) = standard_air(20_000.0);
    let then = released::atmosphere::density_matched_altitude_m(
        temperature_k - 273.15,
        pressure_pa / 100.0,
        0.0,
    );
    let quote = format!("read standard air at 20 km as {then:.0} m");
    assert_quotes(
        "docs/ATMOSPHERE_DENSITY_ALTITUDE.md",
        std::slice::from_ref(&quote),
    );
    assert_quotes("CHANGELOG.md", std::slice::from_ref(&quote));

    let v = 1e200;
    let inputs = released::BallisticInputs {
        bc_value: 0.5,
        bc_type: released::DragModel::G1,
        bullet_mass: 0.01,
        bullet_diameter: 0.00782,
        muzzle_velocity: v,
        target_distance: 1000.0,
        ..released::BallisticInputs::default()
    };
    let overflowed = released::fast_trajectory::fast_integrate_with_segments(
        &inputs,
        vec![],
        released::fast_trajectory::FastIntegrationParams {
            horiz: 1000.0,
            vert: 0.0,
            initial_state: [0.0, 0.0, 0.0, v, 0.0, 0.0],
            t_span: (0.0, 10.0),
            atmo_params: (0.0, 15.0, 1013.25, 0.0),
            atmo_sock: None,
        },
    );
    assert!(overflowed.success && overflowed.t.len() == 1);
    assert_quotes(
        "CHANGELOG.md",
        &["still report `success = true` (a muzzle velocity of 1e200 m/s".to_string()],
    );

    let blamed = released::calculate_zero_angle_with_station_presence(
        released::BallisticInputs {
            bc_value: 1.05,
            bc_type: released::DragModel::G1,
            bullet_mass: 0.0486,
            bullet_diameter: 0.0127,
            muzzle_velocity: 823.0,
            target_distance: 2000.0,
            temperature: COLD_ATMOSPHERE_K - 273.15,
            pressure: COLD_PRESSURE_PA / 100.0,
            humidity: 0.0,
            ..released::BallisticInputs::default()
        },
        2000.0,
        1000.0,
        released::WindConditions::default(),
        released::AtmosphericConditions {
            temperature: COLD_ATMOSPHERE_K - 273.15,
            pressure: COLD_PRESSURE_PA / 100.0,
            humidity: 0.0,
            altitude: 0.0,
        },
        released::StationAtmospherePresence {
            temperature_supplied: true,
            pressure_supplied: true,
        },
    )
    .expect_err("0.49.0 found no zero")
    .to_string();
    assert_quotes("CHANGELOG.md", &[format!("report \"{blamed}\"")]);
}

// --- the CLI -------------------------------------------------------------------------------------

#[cfg(feature = "cli")]
mod cli {
    use super::ALOFT;
    use std::process::{Command, Output};
    use std::sync::atomic::{AtomicU32, Ordering};

    const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

    /// Runs the binary with `HOME` at a fresh private directory, so no saved profile of the
    /// machine's user can reach the run.
    fn run(args: &[&str]) -> Output {
        static N: AtomicU32 = AtomicU32::new(0);
        let home = std::env::temp_dir().join(format!(
            "bx-air-temperature-{}-{}",
            std::process::id(),
            N.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&home).unwrap();
        let output = Command::new(BIN)
            .env("HOME", &home)
            .args(args)
            .output()
            .expect("spawn ballistics");
        let _ = std::fs::remove_dir_all(&home);
        output
    }

    fn trajectory(units: &str, air: &[&str]) -> Output {
        let mut args = vec![
            "trajectory",
            "--units",
            units,
            "-v",
            "823",
            "-b",
            "0.475",
            "--max-range",
            "300",
        ];
        args.extend_from_slice(match units {
            "metric" => &["-m", "10.9", "-d", "7.82"],
            _ => &["-m", "168", "-d", "0.308"],
        });
        args.extend_from_slice(air);
        run(&args)
    }

    fn stderr(output: &Output) -> String {
        String::from_utf8_lossy(&output.stderr).into_owned()
    }

    /// `--temperature` takes every reading the solve flies: 250 °C, 482 °F, -150 °C and -238 °F
    /// used to be refused as "out of range"; a --location cell or --zero-temperature with the
    /// same reading flew.
    #[test]
    fn the_temperature_flag_takes_every_reading_the_solve_flies() {
        for (units, temperature) in [
            ("metric", "250"),
            ("metric", "-150"),
            ("imperial", "482"),
            ("imperial", "-238"),
        ] {
            let flag = format!("--temperature={temperature}");
            let output = trajectory(units, &[&flag, "--humidity", "0"]);
            assert_eq!(
                output.status.code(),
                Some(0),
                "--temperature {temperature} ({units}) flies:\n{}",
                stderr(&output)
            );
        }
    }

    /// At or below absolute zero is the one temperature the flag refuses, in its own unit.
    #[test]
    fn the_temperature_flag_refuses_only_air_at_or_below_absolute_zero() {
        for (units, temperature, cause) in [
            (
                "imperial",
                "-500",
                "Error: --temperature -500 F is at or below absolute zero (-459.67 F)",
            ),
            (
                "metric",
                "-273.15",
                "Error: --temperature -273.15 C is at or below absolute zero (-273.15 C)",
            ),
            (
                "metric",
                "NaN",
                "Error: --temperature NaN C is not a finite number",
            ),
        ] {
            let flag = format!("--temperature={temperature}");
            let output = trajectory(units, &[&flag]);
            assert_eq!(output.status.code(), Some(1), "{temperature} {units}");
            let stderr = stderr(&output);
            assert!(
                stderr.lines().next().unwrap_or_default().starts_with(cause),
                "{temperature} {units}: expected `{cause}`, got:\n{stderr}"
            );
        }
    }

    /// The examples CLI_USAGE.md gives for the one temperature rule are the binary's own output.
    #[test]
    fn the_cli_usage_examples_are_engine_output() {
        let usage = super::normalized(&super::read("CLI_USAGE.md"));
        let quoted = |line: &str| {
            assert!(
                usage.contains(&super::normalized(line)),
                "CLI_USAGE.md does not quote `{line}`"
            )
        };
        let base = [
            "trajectory",
            "-v",
            "2700",
            "-b",
            "0.475",
            "-m",
            "168",
            "-d",
            "0.308",
            "--max-range",
            "300",
        ];
        let hot = run(&[&base[..], &["--temperature", "482", "--humidity", "0"]].concat());
        assert_eq!(hot.status.code(), Some(0));
        let stdout = String::from_utf8_lossy(&hot.stdout).into_owned();
        for line in stdout.lines().filter(|line| line.contains("Impact")) {
            quoted(line);
        }
        let cold = run(&[&base[..], &["--temperature=-500"]].concat());
        quoted(stderr(&cold).lines().next().unwrap_or_default());
        let lob = |angle: &str, output: &[&str]| {
            run(&[
                &[
                    "trajectory",
                    "-u",
                    "metric",
                    "-v",
                    "823",
                    "-b",
                    "1.05",
                    "-m",
                    "48.6",
                    "-d",
                    "12.7",
                    "--drag-model",
                    "g1",
                    "--max-range",
                    "20000",
                    "--angle",
                    angle,
                    "--temperature=-263.15",
                    "--pressure",
                    "250",
                    "--humidity",
                    "0",
                ][..],
                output,
            ]
            .concat())
        };
        let refused = lob("40", &[]);
        let error = stderr(&refused);
        quoted(
            error
                .lines()
                .find(|line| line.starts_with("Error: "))
                .unwrap_or_default(),
        );
        let flat = lob("1", &["-o", "json"]);
        let summary: serde_json::Value = serde_json::from_slice(&flat.stdout).expect("json");
        quoted(&format!(
            "(`--angle 1`: {:.0} m, apex {:.1} m)",
            summary["max_range"].as_f64().unwrap(),
            summary["max_height"].as_f64().unwrap()
        ));
        // The same lob under --units imperial, as CLI_USAGE.md abridges it.
        let imperial = run(&[
            "trajectory",
            "-v",
            "2700",
            "-b",
            "1.05",
            "-m",
            "750",
            "-d",
            "0.5",
            "--drag-model",
            "g1",
            "--max-range",
            "21872",
            "--angle",
            "40",
            "--temperature=-441.67",
            "--pressure",
            "8",
            "--humidity",
            "0",
        ]);
        let error = stderr(&imperial);
        let error = error
            .lines()
            .find_map(|line| line.strip_prefix("Error: "))
            .unwrap_or_default();
        let (reading, rest) = error
            .split_once(": the air aloft")
            .expect("an aloft refusal");
        let height = rest
            .split(" reaches 0 K ")
            .nth(1)
            .and_then(|rest| rest.split(", which").next())
            .unwrap_or_default();
        quoted(&format!("`{reading}: ... {height} ...`"));
        // Saturated air at 104 °F 65617 ft up, with no --pressure: blamed on the altitude, every
        // altitude in feet.
        let high = run(&[
            "trajectory",
            "-v",
            "2700",
            "-b",
            "0.475",
            "-m",
            "168",
            "-d",
            "0.308",
            "--max-range",
            "300",
            "--altitude",
            "65617",
            "--temperature",
            "104",
            "--humidity",
            "100",
        ]);
        let error = stderr(&high);
        let station = error
            .split("the station altitude of ")
            .nth(1)
            .and_then(|rest| rest.split(" is a pressure altitude of").next())
            .unwrap_or_default();
        let bounds = error
            .split("outside the troposphere ")
            .nth(1)
            .and_then(|rest| rest.split(':').next())
            .unwrap_or_default();
        assert_eq!(
            (station, bounds),
            ("65617 ft", "(-16404 ft to 36089 ft)"),
            "{error}"
        );
        quoted(&format!(
            "the station altitude of {station} is a pressure altitude of ... ft, outside the \
             troposphere {bounds} ..."
        ));
    }

    /// A lob from a 10 K station under 250 hPa reaches air below absolute zero: refused naming
    /// the temperature by its flag, in the run's units, not as a non-finite integration state.
    /// It used to read `temperature_c = -263.15` and give the height in metres whatever `--units`
    /// said.
    #[test]
    fn a_cli_lob_into_air_below_absolute_zero_is_refused_naming_the_temperature() {
        for (units, readings, refusal) in [
            (
                "metric",
                [
                    "-v",
                    "823",
                    "-m",
                    "48.6",
                    "-d",
                    "12.7",
                    "--max-range",
                    "20000",
                    "--temperature=-263.15",
                    "--pressure",
                    "250",
                ],
                "Error: --temperature = -263.15 °C: the air aloft would fall below absolute zero: \
                 from this temperature the standard lapse reaches 0 K at 1539 m above mean sea \
                 level, which this flight reaches (check the temperature and its unit)",
            ),
            (
                "imperial",
                [
                    "-v",
                    "2700",
                    "-m",
                    "750",
                    "-d",
                    "0.5",
                    "--max-range",
                    "21872",
                    "--temperature=-441.67",
                    "--pressure",
                    "8",
                ],
                "Error: --temperature = -441.67 °F: the air aloft would fall below absolute zero: \
                 from this temperature the standard lapse reaches 0 K at 5049 ft above mean sea \
                 level, which this flight reaches (check the temperature and its unit)",
            ),
        ] {
            let mut args = vec!["trajectory", "--units", units];
            args.extend_from_slice(&readings);
            args.extend_from_slice(&[
                "-b",
                "1.05",
                "--drag-model",
                "g1",
                "--angle",
                "40",
                "--humidity",
                "0",
            ]);
            let output = run(&args);
            assert_eq!(output.status.code(), Some(1), "{units}");
            // A warning that the pressure is implausible at sea level comes first; then the
            // refusal.
            let stderr = stderr(&output);
            let error = stderr
                .lines()
                .find(|line| line.starts_with("Error: "))
                .unwrap_or_default();
            assert_eq!(error, refusal, "{units}:\n{stderr}");
            assert!(error.contains(ALOFT), "{stderr}");
            assert!(!stderr.contains("non-finite"), "{stderr}");
        }
    }

    /// The CLI flies, or refuses, each named atmosphere exactly when `solve` does: `--temperature`
    /// has no range of its own.
    #[test]
    fn the_cli_flies_each_named_atmosphere_exactly_when_solve_does() {
        let mut flown = 0;
        for (temperature_k, pressure_pa, relative_humidity) in super::named_atmospheres() {
            let temperature = format!("--temperature={}", temperature_k - 273.15);
            let pressure = (pressure_pa / 100.0).to_string();
            let humidity = (relative_humidity * 100.0).to_string();
            let output = trajectory(
                "metric",
                &[
                    &temperature,
                    "--pressure",
                    &pressure,
                    "--humidity",
                    &humidity,
                ],
            );
            let solved = super::solve(&super::solve_request(
                serde_json::json!({"temperature_k": temperature_k, "pressure_pa": pressure_pa,
                                   "relative_humidity": relative_humidity}),
                None,
                "rk45",
            ));
            assert_eq!(
                output.status.success(),
                solved.is_ok(),
                "{temperature_k} K, {pressure_pa} Pa, RH {relative_humidity}: solve {solved:?}, \
                 CLI:\n{}",
                stderr(&output)
            );
            flown += usize::from(solved.is_ok());
        }
        assert_eq!(flown, super::NAMED_ATMOSPHERES_THAT_FLY);
    }
}
