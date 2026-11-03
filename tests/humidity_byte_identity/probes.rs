// The probes of `tests/humidity_byte_identity.rs`, written once and compiled twice: `include!`d
// into its module `current` with `engine` = this crate, and into `released_0_49_0` with `engine` =
// the released 0.49.0 crate (the `ballistics_engine_049` dev-dependency). The few calls whose
// signatures 0.50.0 changed (the kernels' atmosphere, the inputs' humidity field) go through each
// module's `shim`, so both engines are asked exactly the same question. Not a test target of its
// own: cargo builds only `tests/*.rs` and `tests/*/main.rs`.

use engine::atmosphere::{
    calculate_air_density_cimp, calculate_atmosphere, density_matched_altitude_m,
    faa_rule_density_altitude_ft, get_local_atmosphere_humid, moist_speed_of_sound,
    nws_pressure_altitude_ft, resolve_atmosphere_for_density_altitude, AtmoSock,
    StationAtmospherePresence,
};
use engine::fast_trajectory::{fast_integrate, fast_integrate_with_segments, FastSolution};
use engine::monte_carlo::solve_trajectory_for_monte_carlo_with_station_presence;
use engine::spin_drift::{effective_sg_from_inputs, litz_drift_meters};
use engine::stability::compute_stability_coefficient;
use engine::trajectory_integration::try_integrate_trajectory;
use engine::wind::{WindSegment, WindSock};
use engine::{
    calculate_zero_angle_with_station_presence, run_monte_carlo_with_station_presence,
    run_monte_carlo_with_wind_and_direction_std_dev_seeded, AtmosphericConditions, BallisticInputs,
    DragModel, MonteCarloParams, MonteCarloResults, TrajectoryResult, TrajectorySolver,
    WindConditions, DEFAULT_HIT_RADIUS_M,
};
use nalgebra::Vector3;
use serde_json::{json, Value};

/// `true`: `TrajectorySolver::new_with_station_presence(.., BOTH)`, the presence-aware
/// constructor; `false`: the same with `NEITHER`, which is exactly what the presence-blind
/// `TrajectorySolver::new` calls.
fn presence(station: Station) -> StationAtmospherePresence {
    if station.presence_aware {
        StationAtmospherePresence::BOTH
    } else {
        StationAtmospherePresence::NEITHER
    }
}

fn load(method: Method, station: Station, shot: Shot) -> BallisticInputs {
    let (use_rk4, use_adaptive_rk45) = method.flags();
    BallisticInputs {
        bc_value: BC_G7,
        bc_type: DragModel::G7,
        bullet_mass: MASS_KG,
        muzzle_velocity: MUZZLE_VELOCITY_MPS,
        bullet_diameter: DIAMETER_M,
        bullet_length: LENGTH_M,
        twist_rate: TWIST_IN,
        muzzle_angle: match shot {
            Shot::Flat => FLAT_MUZZLE_ANGLE_RAD,
            Shot::Lob => LOB_MUZZLE_ANGLE_RAD,
        },
        sight_height: 0.05,
        altitude: station.altitude_m,
        temperature: station.temperature_c,
        pressure: station.pressure_hpa,
        use_rk4,
        use_adaptive_rk45,
        enable_trajectory_sampling: shot == Shot::Flat,
        sample_interval: 50.0,
        ..BallisticInputs::default()
    }
}

fn crosswind() -> WindConditions {
    WindConditions {
        speed: CROSSWIND_MPS,
        direction: CROSSWIND_FROM_RAD,
        vertical_speed: 0.0,
    }
}

fn atmosphere(station: Station, humidity_percent: f64) -> AtmosphericConditions {
    AtmosphericConditions {
        temperature: station.temperature_c,
        pressure: station.pressure_hpa,
        humidity: humidity_percent,
        altitude: station.altitude_m,
    }
}

fn solver(
    inputs: BallisticInputs,
    station: Station,
    humidity_percent: f64,
    zones: Zones,
    shot: Shot,
) -> TrajectorySolver {
    let mut solver = TrajectorySolver::new_with_station_presence(
        inputs,
        crosswind(),
        atmosphere(station, humidity_percent),
        presence(station),
    );
    solver.set_max_range(match shot {
        Shot::Flat => 1000.0,
        Shot::Lob => 30_000.0,
    });
    if let Some(segments) = zones.segments(station) {
        solver.set_atmo_segments(segments);
    }
    solver
}

/// Record a `TrajectoryResult`: summary fields, five point indices, and digests over every
/// point's time and position and every sample's distance, drop, drift and time.
fn record_result(pins: &mut Pins, result: &TrajectoryResult) {
    pins.f("max_range", result.max_range);
    pins.f("max_height", result.max_height);
    pins.f("time_of_flight", result.time_of_flight);
    pins.f("impact_velocity", result.impact_velocity);
    pins.f("impact_energy", result.impact_energy);
    pins.f(
        "station_speed_of_sound_mps",
        result.station_speed_of_sound_mps,
    );
    pins.opt("mach_1_2_distance_m", result.mach_1_2_distance_m);
    pins.opt("mach_1_0_distance_m", result.mach_1_0_distance_m);
    pins.opt("mach_0_9_distance_m", result.mach_0_9_distance_m);
    pins.opt("transonic_mach", result.transonic_mach);
    pins.text("termination", &format!("{:?}", result.termination));

    let n = result.points.len();
    pins.n("points.len", n as u64);
    for index in [1, n / 4, n / 2, (3 * n) / 4, n - 1] {
        let point = &result.points[index];
        pins.f(format!("points[{index}].time"), point.time);
        pins.f(format!("points[{index}].x"), point.position.x);
        pins.f(format!("points[{index}].y"), point.position.y);
        pins.f(format!("points[{index}].z"), point.position.z);
        pins.f(format!("points[{index}].speed"), point.velocity_magnitude);
        pins.opt(
            &format!("points[{index}].drag_coefficient"),
            point.drag_coefficient,
        );
    }
    // The digest covers time and position only; see the module docs on libm.
    pins.digest(
        "points (time, x, y, z)",
        result.points.iter().flat_map(|point| {
            [
                point.time,
                point.position.x,
                point.position.y,
                point.position.z,
            ]
        }),
    );
    match &result.sampled_points {
        Some(samples) => {
            pins.n("sampled_points.len", samples.len() as u64);
            pins.digest(
                "sampled_points (distance, drop, drift, time, flag count)",
                samples.iter().flat_map(|sample| {
                    [
                        sample.distance_m,
                        sample.drop_m,
                        sample.wind_drift_m,
                        sample.time_s,
                        sample.flags.len() as f64,
                    ]
                }),
            );
        }
        None => pins.n("sampled_points.len (none)", u64::MAX),
    }
}

fn ts_case(method: Method, station: Station, rh: f64, zones: Zones, shot: Shot) -> (String, Pins) {
    let name = format!(
        "{}/{}/{}/rh{rh}/{}",
        shot.name(),
        method.name(),
        station.name,
        zones.name()
    );
    let result = solver(load(method, station, shot), station, rh, zones, shot)
        .solve()
        .unwrap_or_else(|error| panic!("{name}: solve failed: {error}"));
    let mut pins = Pins::default();
    record_result(&mut pins, &result);
    (name, pins)
}

/// The `TrajectorySolver::solve` grid for one shot, split by whether any RH involved is above
/// zero.
fn ts_grid(shot: Shot, humid: bool) -> Probe {
    let mut cases = Vec::new();
    for method in Method::ALL {
        for station in STATIONS {
            let combos: &[(f64, Zones)] = match (shot, humid) {
                (Shot::Flat, false) => &[(0.0, Zones::None), (0.0, Zones::Dry)],
                (Shot::Flat, true) => &[
                    (0.0, Zones::Humid),
                    (50.0, Zones::None),
                    (50.0, Zones::Humid),
                    (100.0, Zones::None),
                    (100.0, Zones::Humid),
                ],
                (Shot::Lob, false) => &[(0.0, Zones::None)],
                (Shot::Lob, true) => &[
                    (50.0, Zones::None),
                    (100.0, Zones::None),
                    (50.0, Zones::Humid),
                ],
            };
            for &(rh, zones) in combos {
                cases.push(ts_case(method, station, rh, zones, shot));
            }
        }
    }
    cases
}

pub(super) fn probe_ts_flat_dry() -> Probe {
    ts_grid(Shot::Flat, false)
}

pub(super) fn probe_ts_flat_humid() -> Probe {
    ts_grid(Shot::Flat, true)
}

pub(super) fn probe_ts_lob_dry() -> Probe {
    ts_grid(Shot::Lob, false)
}

pub(super) fn probe_ts_lob_humid() -> Probe {
    ts_grid(Shot::Lob, true)
}

/// Spin drift, aerodynamic jump and Magnus through `TrajectorySolver` (RK4).
fn ts_effects(humid: bool) -> Probe {
    let mut cases = Vec::new();
    let rhs: &[f64] = if humid { &[50.0, 100.0] } else { &[0.0] };
    for station in STATIONS {
        for &rh in rhs {
            for (effect, configure) in [
                (
                    "spin-drift+aerodynamic-jump",
                    (|inputs: &mut BallisticInputs| {
                        inputs.use_enhanced_spin_drift = true;
                        inputs.enable_aerodynamic_jump = true;
                    }) as fn(&mut BallisticInputs),
                ),
                ("magnus", |inputs: &mut BallisticInputs| {
                    inputs.enable_magnus = true;
                }),
            ] {
                let mut inputs = load(Method::Rk4, station, Shot::Flat);
                configure(&mut inputs);
                let name = format!("{effect}/rk4/{}/rh{rh}", station.name);
                let result = solver(inputs, station, rh, Zones::None, Shot::Flat)
                    .solve()
                    .unwrap_or_else(|error| panic!("{name}: solve failed: {error}"));
                let mut pins = Pins::default();
                record_result(&mut pins, &result);
                cases.push((name, pins));
            }
        }
    }
    cases
}

pub(super) fn probe_ts_effects_dry() -> Probe {
    ts_effects(false)
}

pub(super) fn probe_ts_effects_humid() -> Probe {
    ts_effects(true)
}

/// Zero solves: every zero trial flies `calculate_acceleration`.
fn zero_grid(humid: bool) -> Probe {
    let mut cases = Vec::new();
    let rhs: &[f64] = if humid { &[50.0, 100.0] } else { &[0.0] };
    for method in Method::ALL {
        for station in STATIONS {
            for &rh in rhs {
                let name = format!("zero-300m/{}/{}/rh{rh}", method.name(), station.name);
                let angle = calculate_zero_angle_with_station_presence(
                    load(method, station, Shot::Flat),
                    300.0,
                    0.05,
                    crosswind(),
                    atmosphere(station, rh),
                    presence(station),
                )
                .unwrap_or_else(|error| panic!("{name}: zero failed: {error}"));
                let mut pins = Pins::default();
                pins.f("zero_angle_rad", angle);
                cases.push((name, pins));
            }
        }
    }
    cases
}

pub(super) fn probe_zero_dry() -> Probe {
    zero_grid(false)
}

pub(super) fn probe_zero_humid() -> Probe {
    zero_grid(true)
}

// =================================================================================================
// Monte Carlo (TrajectorySolver sampler), seeded, at the default humidity
// =================================================================================================

fn record_mc(pins: &mut Pins, results: &MonteCarloResults) {
    pins.n("ranges.len", results.ranges.len() as u64);
    pins.f("ranges[0]", results.ranges[0]);
    pins.f("impact_velocities[0]", results.impact_velocities[0]);
    pins.f("impact_positions[0].y", results.impact_positions[0].y);
    pins.f("impact_positions[0].z", results.impact_positions[0].z);
    pins.f(
        "hit_probability(DEFAULT_HIT_RADIUS_M)",
        results.hit_probability(DEFAULT_HIT_RADIUS_M),
    );
    pins.digest("ranges", results.ranges.iter().copied());
    pins.digest(
        "impact_velocities",
        results.impact_velocities.iter().copied(),
    );
    pins.digest(
        "impact_positions",
        results
            .impact_positions
            .iter()
            .flat_map(|position| [position.x, position.y, position.z]),
    );
}

/// The inputs a CLI `monte-carlo` run builds when no humidity is given: the
/// `BallisticInputs` default (50 %), never named here so the pin survives the field's rename.
fn mc_inputs(station: Station) -> BallisticInputs {
    BallisticInputs {
        muzzle_height: 1.5,
        ground_threshold: 0.0,
        ..load(Method::Rk45, station, Shot::Flat)
    }
}

fn mc_params() -> MonteCarloParams {
    MonteCarloParams {
        num_simulations: 40,
        target_distance: Some(500.0),
        base_wind_speed: CROSSWIND_MPS,
        base_wind_direction: CROSSWIND_FROM_RAD,
        ..MonteCarloParams::default()
    }
}

pub(super) fn probe_monte_carlo_default_humidity() -> Probe {
    let mut cases = Vec::new();

    // The CLI's seeded path (presence-blind).
    let results = run_monte_carlo_with_wind_and_direction_std_dev_seeded(
        mc_inputs(SEA_LEVEL),
        crosswind(),
        mc_params(),
        0.01,
        0x5EED_0050,
    )
    .expect("seeded sea-level Monte Carlo");
    let mut pins = Pins::default();
    record_mc(&mut pins, &results);
    cases.push(("seeded/icao-sea-level/default-rh".to_string(), pins));

    // The presence-aware twin (the FFI's core), at 1500 m.
    let results = run_monte_carlo_with_station_presence(
        mc_inputs(HIGH_STATION),
        crosswind(),
        mc_params(),
        0.01,
        Some(0x5EED_1500),
        StationAtmospherePresence::BOTH,
    )
    .expect("seeded 1500 m Monte Carlo");
    let mut pins = Pins::default();
    record_mc(&mut pins, &results);
    cases.push(("seeded/1500m-25C-830hPa/default-rh".to_string(), pins));

    cases
}

// =================================================================================================
// Station air, local air, density altitude, stability
// =================================================================================================

const RH_GRID: [f64; 4] = [0.0, 37.5, 50.0, 100.0];

pub(super) fn probe_station_air() -> Probe {
    let mut cases = Vec::new();
    let overrides: [(Option<f64>, Option<f64>); 8] = [
        (None, None),
        (Some(-10.0), Some(1013.25)),
        (Some(15.0), Some(1013.25)),
        (Some(35.0), Some(1013.25)),
        (Some(25.0), Some(830.0)),
        (Some(35.0), Some(850.0)),
        (Some(30.0), None),
        (None, Some(900.0)),
    ];
    for altitude_m in [0.0, 1500.0, 3000.0] {
        for (temperature_c, pressure_hpa) in overrides {
            let name =
                format!("calculate_atmosphere/{altitude_m}m/{temperature_c:?}/{pressure_hpa:?}");
            let mut pins = Pins::default();
            for rh in RH_GRID {
                let (density, speed_of_sound) =
                    calculate_atmosphere(altitude_m, temperature_c, pressure_hpa, rh);
                pins.f(format!("rh{rh}.density"), density);
                pins.f(format!("rh{rh}.speed_of_sound"), speed_of_sound);
            }
            cases.push((name, pins));
        }
    }
    // The two humid-air primitives the station values are built from.
    for (temperature_c, pressure_hpa) in [(-10.0, 1013.25), (15.0, 1013.25), (35.0, 850.0)] {
        let name = format!("primitives/{temperature_c}C/{pressure_hpa}hPa");
        let mut pins = Pins::default();
        for rh in RH_GRID {
            pins.f(
                format!("rh{rh}.calculate_air_density_cimp"),
                calculate_air_density_cimp(temperature_c, pressure_hpa, rh),
            );
            pins.f(
                format!("rh{rh}.moist_speed_of_sound"),
                moist_speed_of_sound(temperature_c + 273.15, pressure_hpa * 100.0, rh),
            );
        }
        cases.push((name, pins));
    }
    cases
}

/// `get_local_atmosphere_humid` over heights from 200 m below the station to above the
/// tropopause, for each station base. The base ratio is the station's own CIPM ratio at that
/// RH, exactly what `TrajectorySolver` feeds it.
fn local_air(humid: bool) -> Probe {
    let mut cases = Vec::new();
    let bases = [
        (0.0, 15.0, 1013.25),
        (0.0, 35.0, 1013.25),
        (1500.0, 25.0, 830.0),
    ];
    let rhs: &[f64] = if humid { &[37.5, 50.0, 100.0] } else { &[0.0] };
    for (base_alt, base_t, base_p) in bases {
        for &rh in rhs {
            let (station_density, _) =
                calculate_atmosphere(base_alt, Some(base_t), Some(base_p), rh);
            let base_ratio = station_density / 1.225;
            let name =
                format!("get_local_atmosphere_humid/{base_alt}m/{base_t}C/{base_p}hPa/rh{rh}");
            let mut pins = Pins::default();
            pins.f("base_ratio", base_ratio);
            for offset in [-200.0, 0.0, 10.0, 500.0, 3000.0, 12_000.0, 25_000.0] {
                let (density, speed_of_sound) = get_local_atmosphere_humid(
                    base_alt + offset,
                    base_alt,
                    base_t,
                    base_p,
                    base_ratio,
                    rh,
                );
                pins.f(format!("{offset:+}m.density"), density);
                pins.f(format!("{offset:+}m.speed_of_sound"), speed_of_sound);
            }
            cases.push((name, pins));
        }
    }
    cases
}

pub(super) fn probe_local_air_dry() -> Probe {
    local_air(false)
}

pub(super) fn probe_local_air_humid() -> Probe {
    local_air(true)
}

pub(super) fn probe_density_altitude() -> Probe {
    let mut cases = Vec::new();

    let mut pins = Pins::default();
    for pressure_hpa in [700.0, 830.0, 900.0, 1013.25, 1050.0] {
        pins.f(
            format!("nws_pressure_altitude_ft({pressure_hpa})"),
            nws_pressure_altitude_ft(pressure_hpa),
        );
        for temperature_f in [14.0, 59.0, 95.0] {
            pins.f(
                format!("faa_rule_density_altitude_ft({pressure_hpa}, {temperature_f})"),
                faa_rule_density_altitude_ft(pressure_hpa, temperature_f),
            );
        }
    }
    cases.push(("nws+faa_rule".to_string(), pins));

    let mut pins = Pins::default();
    for (temperature_c, pressure_hpa) in [(15.0, 1013.25), (30.0, 1013.25), (25.0, 830.0)] {
        for rh in RH_GRID {
            pins.f(
                format!("density_matched_altitude_m({temperature_c}, {pressure_hpa}, {rh})"),
                density_matched_altitude_m(temperature_c, pressure_hpa, rh),
            );
        }
    }
    cases.push(("density_matched".to_string(), pins));

    let mut pins = Pins::default();
    for density_altitude_m in [-500.0, 0.0, 1500.0, 3000.0] {
        for temperature_c in [None, Some(30.0)] {
            let (altitude_m, resolved_t, resolved_p) =
                resolve_atmosphere_for_density_altitude(density_altitude_m, temperature_c);
            let label = format!("resolve({density_altitude_m}, {temperature_c:?})");
            pins.f(format!("{label}.altitude_m"), altitude_m);
            pins.f(format!("{label}.temperature_c"), resolved_t);
            pins.f(format!("{label}.pressure_hpa"), resolved_p);
        }
    }
    cases.push(("resolve_atmosphere_for_density_altitude".to_string(), pins));

    // The service behind the `atmosphere.density_altitude` bridge command.
    for (name, atmosphere) in [
        (
            "service/icao-sea-level/rh0.5",
            json!({"altitude_m": 0.0, "temperature_k": 288.15, "pressure_pa": 101325.0,
                   "relative_humidity": 0.5}),
        ),
        (
            "service/1500m/rh0",
            json!({"altitude_m": 1500.0, "temperature_k": 298.15, "pressure_pa": 83000.0,
                   "relative_humidity": 0.0}),
        ),
        (
            "service/1500m/rh0.9",
            json!({"altitude_m": 1500.0, "temperature_k": 298.15, "pressure_pa": 83000.0,
                   "relative_humidity": 0.9}),
        ),
    ] {
        let request = engine::atmosphere_service::decode_density_altitude_request_v1(
            &json!({ "atmosphere": atmosphere }),
        )
        .expect("decode density-altitude request");
        let response =
            engine::atmosphere_service::density_altitude_v1(&request).expect("density altitude");
        let mut pins = Pins::default();
        pins.f("pressure_altitude.m", response.pressure_altitude.m);
        pins.f("pressure_altitude.ft", response.pressure_altitude.ft);
        pins.f("faa_rule.m", response.density_altitude.faa_rule.m);
        pins.f("faa_rule.ft", response.density_altitude.faa_rule.ft);
        pins.f(
            "density_matched.m",
            response.density_altitude.density_matched.m,
        );
        pins.f(
            "density_matched.ft",
            response.density_altitude.density_matched.ft,
        );
        pins.f("air_density_kg_m3", response.air_density_kg_m3);
        cases.push((name.to_string(), pins));
    }
    cases
}

/// Muzzle Miller Sg (deliberately dry: the published Miller T/P correction), Litz spin drift
/// and aerodynamic jump.
pub(super) fn probe_stability_spin_jump() -> Probe {
    let mut cases = Vec::new();
    let inputs = load(Method::Rk4, SEA_LEVEL, Shot::Flat);

    let mut pins = Pins::default();
    for (temperature_c, pressure_hpa) in [(15.0, 1013.25), (35.0, 1013.25), (25.0, 830.0)] {
        let sg = compute_stability_coefficient(&inputs, (0.0, temperature_c, pressure_hpa, 0.0));
        pins.f(
            format!("compute_stability_coefficient({temperature_c}, {pressure_hpa})"),
            sg,
        );
        let sg_eff = effective_sg_from_inputs(&inputs, temperature_c, pressure_hpa);
        pins.f(
            format!("effective_sg_from_inputs({temperature_c}, {pressure_hpa})"),
            sg_eff,
        );
        for time_s in [0.5, 1.0, 2.0] {
            pins.f(
                format!("litz_drift_meters(sg_eff, {time_s})"),
                litz_drift_meters(sg_eff, time_s, true),
            );
        }
    }
    // The estimated-length fallback inside effective_sg_from_inputs.
    let no_length = BallisticInputs {
        bullet_length: 0.0,
        ..inputs.clone()
    };
    pins.f(
        "effective_sg_from_inputs(no length, 15, 1013.25)",
        effective_sg_from_inputs(&no_length, 15.0, 1013.25),
    );
    cases.push(("miller-sg+litz".to_string(), pins));

    // fast_integrate's aerodynamic-jump launch offset, standard and direct atmosphere.
    let jump_inputs = BallisticInputs {
        enable_aerodynamic_jump: true,
        wind_speed: CROSSWIND_MPS,
        wind_angle: CROSSWIND_FROM_RAD,
        ..inputs.clone()
    };
    let mut pins = Pins::default();
    for (label, atmo) in [
        ("standard(0, 15, 1013.25, 1.0)", (0.0, 15.0, 1013.25, 1.0)),
        ("standard(1500, 25, 830, 0.79)", (1500.0, 25.0, 830.0, 0.79)),
        ("direct(1.1, 335, 0, 0)", (1.1, 335.0, 0.0, 0.0)),
    ] {
        pins.f(
            format!("aerodynamic_jump_launch_offset_rad {label}"),
            shim::aerodynamic_jump_launch_offset_rad(&jump_inputs, atmo),
        );
    }
    cases.push(("fast-path-aerodynamic-jump".to_string(), pins));

    // TrajectorySolver's aerodynamic-jump components. They depend on the station T/P only,
    // never on RH or the flight, so every RH is pinned here as IMMUTABLE.
    for station in STATIONS {
        for rh in [0.0, 50.0, 100.0] {
            let mut jump_load = load(Method::Rk45, station, Shot::Flat);
            jump_load.enable_aerodynamic_jump = true;
            let result = solver(jump_load, station, rh, Zones::None, Shot::Flat)
                .solve()
                .expect("aerodynamic-jump solve");
            let jump = result.aerodynamic_jump.expect("jump components");
            let mut pins = Pins::default();
            pins.f("vertical_jump_moa", jump.vertical_jump_moa);
            pins.f("horizontal_jump_moa", jump.horizontal_jump_moa);
            pins.f("jump_angle_rad", jump.jump_angle_rad);
            pins.f("magnus_component_moa", jump.magnus_component_moa);
            pins.f("yaw_component_moa", jump.yaw_component_moa);
            pins.f("stabilization_factor", jump.stabilization_factor);
            cases.push((
                format!("trajectory-solver-jump/{}/rh{rh}", station.name),
                pins,
            ));
        }
    }
    cases
}

// =================================================================================================
// Kernels: direct ("uniform") mode -- bits -- and standard mode at RH 0 -- relative 1e-15
// =================================================================================================

/// The kernel's legacy direct-mode tuple: `(density, speed_of_sound, 0, 0)`.
const DIRECT: (f64, f64, f64, f64) = (1.1, 335.0, 0.0, 0.0);

/// A kernel station: the legacy standard-mode tuple `(alt, T, P, ratio)` with the ratio taken
/// from `calculate_atmosphere` at RH 0, which is exactly what a dry station built from the
/// same alt/T/P yields.
fn dry_standard_tuple(station: Station) -> (f64, f64, f64, f64) {
    let (density, _) = calculate_atmosphere(
        station.altitude_m,
        Some(station.temperature_c),
        Some(station.pressure_hpa),
        0.0,
    );
    (
        station.altitude_m,
        station.temperature_c,
        station.pressure_hpa,
        density / 1.225,
    )
}

fn kernel_inputs(station: Station) -> BallisticInputs {
    BallisticInputs {
        altitude: station.altitude_m,
        temperature: station.temperature_c,
        pressure: station.pressure_hpa,
        muzzle_angle: FLAT_MUZZLE_ANGLE_RAD,
        wind_speed: CROSSWIND_MPS,
        wind_angle: CROSSWIND_FROM_RAD,
        enable_trajectory_sampling: false,
        ..load(Method::Rk45, station, Shot::Flat)
    }
}

fn launch_state() -> [f64; 6] {
    [
        0.0,
        0.0,
        0.0,
        MUZZLE_VELOCITY_MPS * FLAT_MUZZLE_ANGLE_RAD.cos(),
        MUZZLE_VELOCITY_MPS * FLAT_MUZZLE_ANGLE_RAD.sin(),
        0.0,
    ]
}

fn crosswind_segments() -> Vec<WindSegment> {
    vec![WindSegment::new(CROSSWIND_MPS * 3.6, 90.0, 2000.0)]
}

/// Indices pinned explicitly for a kernel trajectory of `n` points.
fn kernel_indices(n: usize) -> [usize; 6] {
    [1, n / 8, n / 4, n / 2, (3 * n) / 4, n - 1]
}

fn record_states(pins: &mut Pins, label: &str, times: &[f64], states: &[[f64; 6]], digest: bool) {
    let n = times.len();
    pins.n(format!("{label}.len"), n as u64);
    for index in kernel_indices(n) {
        pins.f(format!("{label}[{index}].t"), times[index]);
        for (component, value) in ["x", "y", "z", "vx", "vy", "vz"].iter().zip(states[index]) {
            pins.f(format!("{label}[{index}].{component}"), value);
        }
    }
    if digest {
        pins.digest(
            label,
            times
                .iter()
                .zip(states)
                .flat_map(|(t, state)| std::iter::once(*t).chain(state.iter().copied())),
        );
    }
}

fn record_fast_solution(pins: &mut Pins, solution: &FastSolution, digest: bool) {
    pins.n("success", u64::from(solution.success));
    let states: Vec<[f64; 6]> = (0..solution.t.len())
        .map(|i| std::array::from_fn(|component| solution.y[component][i]))
        .collect();
    record_states(pins, "fast", &solution.t, &states, digest);
    for (slot, events) in solution.t_events.iter().enumerate() {
        pins.n(format!("t_events[{slot}].len"), events.len() as u64);
        for (index, time) in events.iter().enumerate() {
            pins.f(format!("t_events[{slot}][{index}]"), *time);
        }
    }
}

/// Run every kernel family once for an atmosphere and record it.
fn kernel_cases(
    prefix: &str,
    station: Station,
    atmos_params: (f64, f64, f64, f64),
    zones: Zones,
    digest: bool,
) -> Probe {
    let mut cases = Vec::new();
    let sock = || zones.segments(station).map(AtmoSock::new);
    let inputs = kernel_inputs(station);

    // compute_derivatives, at the station and 300 m above it.
    let mut pins = Pins::default();
    for (label, position) in [
        ("muzzle", Vector3::new(0.0, 0.0, 0.0)),
        ("500m-downrange", Vector3::new(500.0, -1.0, 0.2)),
        ("300m-above", Vector3::new(600.0, 300.0, 0.0)),
    ] {
        let derivative = shim::compute_derivatives(
            position,
            Vector3::new(700.0, 2.0, 0.5),
            &inputs,
            Vector3::new(0.0, 0.0, -CROSSWIND_MPS),
            atmos_params,
            BC_G7,
            sock().as_ref(),
        );
        for (component, value) in ["vx", "vy", "vz", "ax", "ay", "az"].iter().zip(derivative) {
            pins.f(format!("{label}.{component}"), value);
        }
    }
    cases.push((format!("{prefix}/compute_derivatives"), pins));

    // try_integrate_trajectory, RK4 and RK45.
    for (method, tolerance, max_step) in [("RK4", 1e-6, 0.001), ("RK45", 1e-6, 0.01)] {
        let trajectory = try_integrate_trajectory(
            launch_state(),
            (0.0, 3.0),
            shim::trajectory_params(atmos_params, sock()),
            method,
            tolerance,
            max_step,
        )
        .expect("valid wind segments");
        let times: Vec<f64> = trajectory.iter().map(|(t, _)| *t).collect();
        let states: Vec<[f64; 6]> = trajectory
            .iter()
            .map(|(_, state)| std::array::from_fn(|component| state[component]))
            .collect();
        let mut pins = Pins::default();
        record_states(&mut pins, "trajectory", &times, &states, digest);
        cases.push((format!("{prefix}/try_integrate_trajectory/{method}"), pins));
    }

    // try_solve_trajectory_rust (the legacy binding entry; it takes no zones).
    if zones == Zones::None {
        let points =
            shim::try_solve_trajectory_rust(launch_state(), crosswind_segments(), atmos_params)
                .expect("valid wind segments");
        let times: Vec<f64> = points.iter().map(|point| point["t"]).collect();
        let states: Vec<[f64; 6]> = points
            .iter()
            .map(|point| ["x", "y", "z", "vx", "vy", "vz"].map(|component| point[component]))
            .collect();
        let mut pins = Pins::default();
        record_states(&mut pins, "points", &times, &states, digest);
        cases.push((format!("{prefix}/try_solve_trajectory_rust/RK45"), pins));
    }

    // fast_integrate (fixed step) and fast_integrate_with_segments (RK45).
    let wind_sock = WindSock::new(crosswind_segments());
    let mut pins = Pins::default();
    record_fast_solution(
        &mut pins,
        &fast_integrate(&inputs, &wind_sock, shim::fast_params(atmos_params, sock())),
        digest,
    );
    cases.push((format!("{prefix}/fast_integrate"), pins));

    let mut pins = Pins::default();
    record_fast_solution(
        &mut pins,
        &fast_integrate_with_segments(
            &inputs,
            crosswind_segments(),
            shim::fast_params(atmos_params, sock()),
        ),
        digest,
    );
    cases.push((format!("{prefix}/fast_integrate_with_segments"), pins));

    cases
}

/// Direct mode: no station, no humidity, so every family is IMMUTABLE bit for bit.
pub(super) fn probe_kernel_direct() -> Probe {
    kernel_cases("direct(1.1,335)", SEA_LEVEL, DIRECT, Zones::None, true)
}

/// Standard mode at RH 0 (station and zones dry), plus the lean Monte Carlo at RH 0.
pub(super) fn probe_kernel_standard_rh0() -> Probe {
    let mut cases = Vec::new();
    for station in STATIONS {
        for zones in [Zones::None, Zones::Dry] {
            cases.extend(kernel_cases(
                &format!("standard/{}/rh0/{}", station.name, zones.name()),
                station,
                dry_standard_tuple(station),
                zones,
                false,
            ));
        }

        // The lean Monte Carlo kernel (fast_integrate behind one BallisticInputs) at RH 0.
        let inputs = shim::dry(BallisticInputs {
            target_distance: 800.0,
            ..kernel_inputs(station)
        });
        let output =
            solve_trajectory_for_monte_carlo_with_station_presence(&inputs, presence(station))
                .expect("lean Monte Carlo solve");
        let mut pins = Pins::default();
        pins.f("drop", output.drop);
        pins.f("wind_drift", output.wind_drift);
        pins.f("time", output.time);
        pins.f("velocity", output.velocity);
        pins.f("energy", output.energy);
        pins.f("mach", output.mach);
        pins.f("spin_drift", output.spin_drift);
        pins.f("distance", output.distance);
        cases.push((
            format!(
                "standard/{}/rh0/solve_trajectory_for_monte_carlo",
                station.name
            ),
            pins,
        ));
    }
    cases
}

// =================================================================================================
// Front ends: solve_json -> solve_v1, the bridge `solve`, the FFI single solve, card_service
// =================================================================================================

fn solve_request(
    method: &str,
    altitude_m: f64,
    temperature_k: f64,
    pressure_pa: f64,
    rh: f64,
) -> Value {
    json!({
        "schema_version": 1,
        "projectile": {
            "mass_kg": MASS_KG,
            "diameter_m": DIAMETER_M,
            "length_m": LENGTH_M,
            "drag_model": "G7",
            "ballistic_coefficient": BC_G7
        },
        "rifle": {
            "muzzle_velocity_mps": MUZZLE_VELOCITY_MPS,
            "sight_height_m": 0.05,
            "muzzle_height_m": 1.2,
            "twist_rate_m_per_turn": 0.254,
            "twist_direction": "right"
        },
        "shot": {
            "max_range_m": 600.0,
            "zero_distance_m": 100.0,
            "shooting_angle_rad": 0.0,
            "cant_angle_rad": 0.0,
            "ground_threshold_m": -100.0
        },
        "atmosphere": {
            "altitude_m": altitude_m,
            "temperature_k": temperature_k,
            "pressure_pa": pressure_pa,
            "relative_humidity": rh
        },
        "wind": {
            "speed_mps": CROSSWIND_MPS,
            "direction_from_rad": CROSSWIND_FROM_RAD
        },
        "solver": { "method": method, "time_step_s": 0.001 },
        "effects": {
            "magnus": false,
            "coriolis": false,
            "enhanced_spin_drift": true,
            "aerodynamic_jump": true
        },
        "sampling": { "interval_m": 50.0 }
    })
}

/// The solve-json requests pinned on every JSON front end.
fn solve_requests(humid: bool) -> Vec<(&'static str, Value)> {
    if humid {
        vec![
            (
                "rk4/icao-sea-level/rh0.5",
                solve_request("rk4", 0.0, 288.15, 101_325.0, 0.5),
            ),
            (
                "rk45/1500m-25C-830hPa/rh0.9",
                solve_request("rk45", 1500.0, 298.15, 83_000.0, 0.9),
            ),
        ]
    } else {
        vec![
            (
                "euler/icao-sea-level/rh0",
                solve_request("euler", 0.0, 288.15, 101_325.0, 0.0),
            ),
            (
                "rk45/1500m-25C-830hPa/rh0",
                solve_request("rk45", 1500.0, 298.15, 83_000.0, 0.0),
            ),
        ]
    }
}

fn number(value: &Value, path: &str) -> f64 {
    value
        .pointer(path)
        .and_then(Value::as_f64)
        .unwrap_or_else(|| panic!("missing number at {path}"))
}

/// Pin a solve-json v1 success document (as JSON): the summary and every sample.
fn record_solve_json(pins: &mut Pins, success: &Value) {
    for field in [
        "actual_range_m",
        "maximum_height_m",
        "time_of_flight_s",
        "terminal_speed_mps",
        "terminal_energy_j",
        "stability_factor",
        "spin_drift_m",
        "aerodynamic_jump_moa",
    ] {
        pins.f(
            format!("summary.{field}"),
            number(success, &format!("/summary/{field}")),
        );
    }
    pins.f(
        "resolved_request.shot.muzzle_angle_rad",
        number(success, "/resolved_request/shot/muzzle_angle_rad"),
    );
    let samples = success["samples"].as_array().expect("samples array");
    pins.n("samples.len", samples.len() as u64);
    for index in [1, samples.len() / 2, samples.len() - 1] {
        for field in [
            "distance_m",
            "time_s",
            "speed_mps",
            "drop_m",
            "windage_m",
            "mach",
        ] {
            pins.f(
                format!("samples[{index}].{field}"),
                number(&samples[index], &format!("/{field}")),
            );
        }
    }
    pins.digest(
        "samples (distance, time, drop, windage)",
        samples.iter().flat_map(|sample| {
            ["distance_m", "time_s", "drop_m", "windage_m"]
                .map(|field| number(sample, &format!("/{field}")))
        }),
    );
}

fn front_end_solve_json(humid: bool) -> Probe {
    let mut cases = Vec::new();
    for (name, request) in solve_requests(humid) {
        let decoded = engine::decode_solve_request_v1(&request.to_string())
            .expect("valid solve-json request");
        let success = engine::solve_v1(decoded).expect("solve_v1 succeeds");
        let success = serde_json::to_value(&success).expect("serialize success");
        let mut pins = Pins::default();
        record_solve_json(&mut pins, &success);
        cases.push((format!("solve_v1/{name}"), pins));
    }
    cases
}

#[cfg(feature = "bridge")]
fn front_end_bridge(humid: bool) -> Probe {
    let mut cases = Vec::new();
    for (name, request) in solve_requests(humid) {
        let envelope = json!({"api_version": 1, "command": "solve", "request": request});
        let raw = engine::bridge::bridge_call(&envelope.to_string());
        let out: Value = serde_json::from_str(&raw).expect("bridge output parses");
        assert_eq!(out["ok"], true, "bridge solve failed: {out}");
        let mut pins = Pins::default();
        record_solve_json(&mut pins, &out["result"]);
        cases.push((format!("bridge-solve/{name}"), pins));
    }
    cases
}

fn front_end_card(humid: bool) -> Probe {
    let mut cases = Vec::new();
    let rhs: &[f64] = if humid { &[50.0, 100.0] } else { &[0.0] };
    for (station, altitude_m, temperature_f, pressure_inhg) in [
        ("icao-sea-level", 0.0, 59.0, 29.92),
        ("1500m-77F-24.51inHg", 1500.0, 77.0, 24.51),
    ] {
        for &rh in rhs {
            let request: engine::card_service::CardRequestV1 = serde_json::from_value(json!({
                "units": "imperial",
                "muzzle_velocity": 2700.0,
                "ballistic_coefficient": BC_G7,
                "mass": 175.0,
                "diameter": 0.308,
                "drag_model": "g7",
                "sight_height": 1.5,
                "zero_distance": 100.0,
                "temperature": temperature_f,
                "pressure": pressure_inhg,
                "humidity": rh,
                "altitude": altitude_m,
                "wind_speed": 10.0,
                "wind_direction_deg": 90.0,
                "start": 100.0,
                "end": 800.0,
                "step": 100.0,
                "adjustment_unit": "mil"
            }))
            .expect("card request");
            let response = engine::card_service::come_ups_v1(&request).expect("come-ups card");
            let response = serde_json::to_value(&response).expect("serialize card");
            let rows = response["rows"].as_array().expect("rows");
            let mut pins = Pins::default();
            pins.n("rows.len", rows.len() as u64);
            let fields = [
                "range",
                "drop_linear",
                "drop_adj",
                "come_up",
                "wind_linear",
                "wind_adj",
                "velocity",
                "energy",
                "time",
            ];
            for (index, row) in rows.iter().enumerate() {
                for field in fields {
                    pins.opt(&format!("rows[{index}].{field}"), row[field].as_f64());
                }
            }
            cases.push((format!("card.come_ups/{station}/rh{rh}"), pins));
        }
    }
    cases
}

pub(super) fn probe_front_ends_dry() -> Probe {
    let mut cases = front_end_solve_json(false);
    cases.extend(front_end_card(false));
    cases
}

pub(super) fn probe_front_ends_humid() -> Probe {
    let mut cases = front_end_solve_json(true);
    cases.extend(front_end_card(true));
    cases
}

#[cfg(feature = "bridge")]
pub(super) fn probe_bridge_solve_dry() -> Probe {
    front_end_bridge(false)
}

#[cfg(feature = "bridge")]
pub(super) fn probe_bridge_solve_humid() -> Probe {
    front_end_bridge(true)
}
