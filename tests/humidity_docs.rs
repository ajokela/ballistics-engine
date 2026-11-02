//! The public docs say what the humidity code does (humidity-consistency train, step S7, plan
//! C6: H19 stale comments, H20 docs that contradict the code).
//!
//! Every number the docs give for the humidity change is recomputed here by running the engine,
//! so a doc can neither quote an estimate nor drift from the code:
//!
//! - **Sound speed.** The kernels' station sound speed before 0.50.0 (dry, 0 %) and now (the
//!   station's own humidity), at 15 °C / 50 %, 30 °C / 50 % and 35 °C / 90 %, and the change.
//! - **Trajectories.** `try_integrate_trajectory` RK4 at 1 ms: a 2 s flat shot and a 30-degree,
//!   30 s lob, before and after. "Before" is 0.49.0's kernel air rebuilt on this code,
//!   `KernelAtmosphere::from_legacy_params` with the humid density ratio and
//!   `RelativeHumidity::DRY`; that reproduction was checked against the 0.49.0 build itself bit
//!   for bit, and [`legacy_kernel_air_reproduces_0_49_0`] keeps it honest.
//! - **The aloft vapour model (S8).** The air at +1 km and +3 km, and three `TrajectorySolver`
//!   shots (flat to 1000 m, a 30-degree lob, and both from a 1500 m station), before and after.
//!   "Before" is 0.49.0, which the train before S8 matched bit for bit there; its values were
//!   flown by the 0.49.0 and 4281ccf builds themselves (the S8 measurement probe), and
//!   [`the_air_before_s8_is_the_old_aloft_model`] rebuilds the air from the old formulas.
//! - **The Miller exception.** The factor a humid density would put on the muzzle Miller Sg,
//!   which 0.50.0 deliberately does not apply. A pin shows the solve's Sg does not move with
//!   humidity.
//! - **Constants.** The engine's own dry sound speed at 15 °C, which `SPEED_OF_SOUND_MPS`'s doc
//!   states.
//!
//! It also checks the CHANGELOG's 0.50.0 section carries the breaking API list, the sensor band
//! from its constant and the downstream follow-ups, and that the stale wording H19/H20 named is
//! gone.

use std::fs;
use std::path::Path;

use ballistics_engine::atmosphere::{calculate_air_density_cimp, moist_speed_of_sound};
use ballistics_engine::trajectory_integration::{try_integrate_trajectory, TrajectoryParams};
use ballistics_engine::{
    air_at, calculate_zero_angle_with_station_presence, decode_solve_request_v1, solve_v1,
    AtmosphericConditions, BallisticInputs, DragModel, KernelAtmosphere, LocalAir,
    RelativeHumidity, StationAir, StationAtmospherePresence, TrajectorySolver, WindConditions,
    HUMIDITY_SENSOR_TOLERANCE_PERCENT,
};

fn read(relative: &str) -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join(relative);
    fs::read_to_string(&path).unwrap_or_else(|error| panic!("read {}: {error}", path.display()))
}

/// The CHANGELOG's `## [Unreleased]` section, up to the next release heading.
fn unreleased() -> String {
    let changelog = read("CHANGELOG.md");
    let start = changelog
        .find("\n## [Unreleased]")
        .expect("CHANGELOG.md has no `## [Unreleased]` section for 0.50.0");
    let body = &changelog[start + 1..];
    let end = body[3..].find("\n## [").map_or(body.len(), |at| at + 3);
    body[..end].to_string()
}

fn assert_contains(haystack: &str, needle: &str, what: &str) {
    assert!(haystack.contains(needle), "{what}: missing `{needle}`");
}

// =================================================================================================
// The CHANGELOG's 0.50.0 section
// =================================================================================================

#[test]
fn the_changelog_has_an_unreleased_section_above_the_last_release() {
    let changelog = read("CHANGELOG.md");
    let unreleased = changelog
        .find("\n## [Unreleased]")
        .expect("no `## [Unreleased]` section");
    let last_release = changelog.find("\n## [0.49.0]").expect("no 0.49.0 section");
    assert!(
        unreleased < last_release,
        "[Unreleased] must sit above 0.49.0"
    );
}

#[test]
fn the_unreleased_section_lists_every_breaking_rust_api_change() {
    let section = unreleased();
    for name in [
        // TrajectoryParams.atmos_params / FastIntegrationParams.atmo_params -> atmosphere.
        "atmos_params",
        "atmo_params",
        "KernelAtmosphere",
        "TrajectoryParams",
        "FastIntegrationParams",
        // The signatures that take the kernel atmosphere.
        "compute_derivatives",
        "try_solve_trajectory_rust",
        "aerodynamic_jump_launch_offset_rad",
        // The checked kernel's error type.
        "try_integrate_trajectory",
        "KernelInputError",
        // The field rename and the removed method.
        "humidity_percent",
        "humidity_percent()",
        // Deprecated, with its sound-speed shift.
        "get_local_atmosphere",
        "4.23e-7",
        // Height measured from the kernel atmosphere's station, no longer from the inputs.
        "inputs.altitude",
        // The new crate-root names, and the glob-import ambiguity.
        "RelativeHumidity",
        "StationAir",
        "LocalAir",
        "air_at",
        "AtmosphereError",
        "DEFAULT_HUMIDITY_PERCENT",
        "InputAxis",
    ] {
        assert_contains(&section, name, "CHANGELOG [Unreleased] breaking list");
    }
}

#[test]
fn the_unreleased_section_states_the_sensor_band_from_its_constant() {
    let section = unreleased();
    let band = HUMIDITY_SENSOR_TOLERANCE_PERCENT;
    assert_contains(
        &section,
        "HUMIDITY_SENSOR_TOLERANCE_PERCENT",
        "the band constant",
    );
    for edge in [100.0 + band, -band, 1.0 + band / 100.0, -band / 100.0] {
        assert_contains(&section, &format!("{edge}"), "a sensor-band edge");
    }
}

#[test]
fn the_unreleased_section_lists_the_downstream_follow_ups() {
    let section = unreleased();
    for repo in [
        "rockets-engine",
        "ballistics-engine-py",
        "ballistics-engine-rb",
        "Apps and services that embed the engine",
        "handbook",
    ] {
        assert_contains(&section, repo, "a downstream follow-up");
    }
}

// =================================================================================================
// Numbers the CHANGELOG quotes, from the engine
// =================================================================================================

const SEA_LEVEL_HPA: f64 = 1013.25;
/// The humid conditions the CHANGELOG quotes, at sea level: (temperature C, RH percent).
const QUOTED: [(f64, f64); 3] = [(15.0, 50.0), (30.0, 50.0), (35.0, 90.0)];

fn station(temperature_c: f64, humidity_percent: f64) -> StationAir {
    StationAir::new(
        0.0,
        temperature_c,
        SEA_LEVEL_HPA,
        RelativeHumidity::try_percent(humidity_percent).unwrap(),
    )
    .unwrap()
}

#[test]
fn the_unreleased_sound_speed_shifts_are_engine_output() {
    let section = unreleased();
    for (temperature_c, humidity_percent) in QUOTED {
        // 0.49.0's kernels flew the station's sound speed at 0 %; 0.50.0's fly its own humidity.
        let dry = station(temperature_c, 0.0).speed_of_sound_mps();
        let moist = station(temperature_c, humidity_percent).speed_of_sound_mps();
        let what = format!("sound speed at {temperature_c} C / {humidity_percent} %");
        assert_contains(&section, &format!("{dry:.3}"), &what);
        assert_contains(&section, &format!("{moist:.3}"), &what);
        assert_contains(
            &section,
            &format!("{:+.3}%", (moist / dry - 1.0) * 100.0),
            &what,
        );
    }
}

const MASS_KG: f64 = 0.011_34; // .308 175 gr
const DIAMETER_M: f64 = 0.007_82;
const LENGTH_M: f64 = 0.031_5; // 1.24 in
const MUZZLE_VELOCITY_MPS: f64 = 792.48; // 2600 fps
const BC_G7: f64 = 0.243;

fn params(atmosphere: KernelAtmosphere) -> TrajectoryParams {
    TrajectoryParams {
        mass_kg: MASS_KG,
        bc: BC_G7,
        drag_model: DragModel::G7,
        wind_segments: vec![],
        atmosphere,
        omega_vector: None,
        enable_spin_drift: false,
        enable_magnus: false,
        enable_coriolis: false,
        target_distance_m: 1.0e6,
        enable_wind_shear: false,
        wind_shear_model: "none".into(),
        shooter_altitude_m: 0.0,
        is_twist_right: true,
        shooting_angle: 0.0,
        bullet_diameter: DIAMETER_M,
        bullet_length: LENGTH_M,
        twist_rate: 10.0,
        custom_drag_table: None,
        cd_scale: 1.0,
        custom_drag_reference: false,
        bc_segments: None,
        use_bc_segments: false,
        bc_segments_data: None,
        ground_threshold: -1.0e5,
        atmo_sock: None,
    }
}

/// 0.49.0's kernel air at a humid station: the humid density ratio with a 0 % sound speed.
fn legacy_kernel_air(temperature_c: f64, humidity_percent: f64) -> KernelAtmosphere {
    let ratio = station(temperature_c, humidity_percent).density_ratio();
    KernelAtmosphere::from_legacy_params(
        (0.0, temperature_c, SEA_LEVEL_HPA, ratio),
        RelativeHumidity::DRY,
    )
    .unwrap()
}

/// `(downrange m, height m, speed m/s)` at the end of an RK4 1 ms flight.
fn fly(atmosphere: KernelAtmosphere, angle_rad: f64, seconds: f64) -> (f64, f64, f64) {
    let launch = [
        0.0,
        0.0,
        0.0,
        MUZZLE_VELOCITY_MPS * angle_rad.cos(),
        MUZZLE_VELOCITY_MPS * angle_rad.sin(),
        0.0,
    ];
    let steps = try_integrate_trajectory(
        launch,
        (0.0, seconds),
        params(atmosphere),
        "RK4",
        1e-6,
        0.001,
    )
    .expect("valid kernel input");
    let (time, state) = steps.last().unwrap();
    assert_eq!(*time, seconds);
    let speed = (state[3] * state[3] + state[4] * state[4] + state[5] * state[5]).sqrt();
    (state[0], state[1], speed)
}

const FLAT: (f64, f64) = (0.002, 2.0);
const LOB: (f64, f64) = (std::f64::consts::FRAC_PI_6, 30.0);

/// The rebuilt "before" is 0.49.0's kernel: the 0.49.0 build itself (main 261c1eb) flew these
/// at 15 C / 50 % (the S7 measurement probe). 1e-12 relative leaves room for another libm.
#[test]
fn legacy_kernel_air_reproduces_0_49_0() {
    let close = |got: f64, want: f64| (got - want).abs() <= 1e-12 * want.abs();
    let (x, y, v) = fly(legacy_kernel_air(15.0, 50.0), FLAT.0, FLAT.1);
    assert!(close(x, 983.029_417_355_035_3), "flat x {x:?}");
    assert!(close(y, -12.733_321_151_623_656), "flat y {y:?}");
    assert!(close(v, 318.355_350_529_326_13), "flat speed {v:?}");
    let (x, y, v) = fly(legacy_kernel_air(15.0, 50.0), LOB.0, LOB.1);
    assert!(close(x, 4_666.747_674_621_04), "lob x {x:?}");
    assert!(close(y, -143.326_524_870_322_97), "lob y {y:?}");
    assert!(close(v, 140.858_365_275_925_73), "lob speed {v:?}");
    // At 0 % the two airs are the same air, so the flights are the same flight.
    assert_eq!(
        fly(legacy_kernel_air(15.0, 0.0), FLAT.0, FLAT.1),
        fly(
            KernelAtmosphere::Station(station(15.0, 0.0)),
            FLAT.0,
            FLAT.1
        )
    );
}

#[test]
fn the_unreleased_trajectory_shifts_are_engine_output() {
    let section = unreleased();
    for (temperature_c, humidity_percent) in QUOTED {
        let before = legacy_kernel_air(temperature_c, humidity_percent);
        let after = KernelAtmosphere::Station(station(temperature_c, humidity_percent));
        let what = format!("trajectory shift at {temperature_c} C / {humidity_percent} %");

        let (x0, y0, v0) = fly(before, FLAT.0, FLAT.1);
        let (x1, y1, v1) = fly(after, FLAT.0, FLAT.1);
        assert_contains(&section, &format!("{:+.3} m", x1 - x0), &what);
        assert_contains(&section, &format!("{:+.1} mm", (y1 - y0) * 1000.0), &what);
        assert_contains(&section, &format!("{:+.3} m/s", v1 - v0), &what);

        let (x0, y0, _) = fly(before, LOB.0, LOB.1);
        let (x1, y1, _) = fly(after, LOB.0, LOB.1);
        assert_contains(&section, &format!("{:+.3} m", x1 - x0), &what);
        assert_contains(&section, &format!("{:+.3} m", y1 - y0), &what);
    }
}

// =================================================================================================
// S8: one vapour model aloft, on every surface
// =================================================================================================

const SEA_LEVEL_STATION_HPA: f64 = 1013.25;
/// The conditions the aloft table quotes: (temperature C, RH percent).
const ALOFT_QUOTED: [(f64, f64); 2] = [(15.0, 50.0), (35.0, 90.0)];

/// The air over a sea-level station at 1013.25 hPa before S8, as the 0.49.0 build (main 261c1eb)
/// and the 4281ccf build flew it, bit for bit alike: `(T C, RH %, height m, density kg/m3,
/// sound speed m/s)`.
const AIR_BEFORE_S8: [(f64, f64, f64, f64, f64); 4] = [
    (15.0, 50.0, 1000.0, 1.1086026577295789, 336.76382227330026),
    (15.0, 50.0, 3000.0, 0.9067541155034896, 328.7463422984116),
    (35.0, 90.0, 1000.0, 1.0267223810704649, 350.29942794649406),
    (35.0, 90.0, 3000.0, 0.8513011892339569, 341.7818362119831),
];

fn air_over_sea_level(temperature_c: f64, humidity_percent: f64, height_m: f64) -> LocalAir {
    let station = StationAir::new(
        0.0,
        temperature_c,
        SEA_LEVEL_STATION_HPA,
        RelativeHumidity::try_percent(humidity_percent).unwrap(),
    )
    .unwrap();
    air_at(&KernelAtmosphere::Station(station), None, 0.0, height_m)
}

/// The constants above are the old model: the dry pressure lapse, the density scaled by
/// `(T_0 P) / (P_0 T)`, and the sound speed re-derived from the held RELATIVE humidity at the
/// local temperature and pressure. Dry air's lapse did not change, so it supplies the old
/// temperature and pressure.
#[test]
fn the_air_before_s8_is_the_old_aloft_model() {
    for (temperature_c, humidity_percent, height_m, density, sound_speed) in AIR_BEFORE_S8 {
        let dry = air_over_sea_level(temperature_c, 0.0, height_m);
        let station = StationAir::new(
            0.0,
            temperature_c,
            SEA_LEVEL_STATION_HPA,
            RelativeHumidity::try_percent(humidity_percent).unwrap(),
        )
        .unwrap();
        let old_sound_speed =
            moist_speed_of_sound(dry.temperature_k, dry.pressure_pa, humidity_percent);
        let old_density = station.density_ratio()
            * ((temperature_c + 273.15) * (dry.pressure_pa / 100.0))
            / (SEA_LEVEL_STATION_HPA * dry.temperature_k)
            * 1.225;
        let close = |got: f64, want: f64| (got - want).abs() <= 1e-12 * want.abs();
        let what = format!("{temperature_c} C / {humidity_percent} % at +{height_m} m");
        assert!(
            close(old_sound_speed, sound_speed),
            "c {what}: {old_sound_speed:?}"
        );
        assert!(
            close(old_density, density),
            "density {what}: {old_density:?}"
        );
    }
}

/// `(station, altitude m, pressure hPa, T C, RH %, flat drop at 1000 m m, flat speed at
/// 1000 m m/s, lob range m, lob apex m, lob impact speed m/s)` before S8, flown by the 0.49.0
/// and 4281ccf builds, bit for bit alike.
#[allow(clippy::type_complexity)]
const SHOTS_BEFORE_S8: [(&str, f64, f64, f64, f64, f64, f64, f64, f64, f64); 4] = [
    (
        "sea level",
        0.0,
        1013.25,
        15.0,
        50.0,
        14.283070868744938,
        315.91839047884184,
        4643.827151319146,
        1107.2037843284834,
        140.22246913957096,
    ),
    (
        "sea level",
        0.0,
        1013.25,
        35.0,
        90.0,
        13.456344837416587,
        333.5755196962042,
        4967.3410410665165,
        1174.9812842154524,
        145.52515063650773,
    ),
    (
        "1500 m",
        1500.0,
        845.6,
        15.0,
        50.0,
        12.234306304109765,
        362.120280940684,
        5241.572104789031,
        1226.5146272252885,
        151.53885237740963,
    ),
    (
        "1500 m",
        1500.0,
        845.6,
        35.0,
        90.0,
        11.666142478948375,
        384.60831144456057,
        5613.809630005237,
        1302.0284360797184,
        157.45728165755804,
    ),
];

const SIGHT_HEIGHT_M: f64 = 0.0381; // 1.5 in

fn shot_inputs(
    altitude_m: f64,
    temperature_c: f64,
    pressure_hpa: f64,
    angle_rad: f64,
) -> BallisticInputs {
    BallisticInputs {
        bc_value: BC_G7,
        bc_type: DragModel::G7,
        bullet_mass: MASS_KG,
        muzzle_velocity: MUZZLE_VELOCITY_MPS,
        bullet_diameter: DIAMETER_M,
        bullet_length: LENGTH_M,
        twist_rate: 10.0,
        muzzle_angle: angle_rad,
        sight_height: SIGHT_HEIGHT_M,
        altitude: altitude_m,
        temperature: temperature_c,
        pressure: pressure_hpa,
        use_rk4: true,
        ..BallisticInputs::default()
    }
}

/// `(drop at 1000 m, speed at 1000 m, lob range, lob apex, lob impact speed)` from
/// `TrajectorySolver` (RK4): the flat shot zeroed at 100 m, the lob at 30 degrees until it falls
/// 100 m below the muzzle.
fn solver_shots(
    altitude_m: f64,
    pressure_hpa: f64,
    temperature_c: f64,
    humidity_percent: f64,
) -> (f64, f64, f64, f64, f64) {
    let air = AtmosphericConditions {
        temperature: temperature_c,
        pressure: pressure_hpa,
        humidity: humidity_percent,
        altitude: altitude_m,
    };
    let presence = StationAtmospherePresence::BOTH;
    let zero = calculate_zero_angle_with_station_presence(
        shot_inputs(altitude_m, temperature_c, pressure_hpa, 0.0),
        100.0,
        SIGHT_HEIGHT_M,
        WindConditions::default(),
        air.clone(),
        presence,
    )
    .unwrap();
    let mut flat = shot_inputs(altitude_m, temperature_c, pressure_hpa, zero);
    flat.enable_trajectory_sampling = true;
    flat.sample_interval = 100.0;
    let mut solver = TrajectorySolver::new_with_station_presence(
        flat,
        WindConditions::default(),
        air.clone(),
        presence,
    );
    solver.set_max_range(1000.0);
    let flat = solver.solve().unwrap();
    let last = flat
        .sampled_points
        .as_ref()
        .unwrap()
        .last()
        .unwrap()
        .clone();
    assert_eq!(last.distance_m, 1000.0);
    let mut solver = TrajectorySolver::new_with_station_presence(
        shot_inputs(
            altitude_m,
            temperature_c,
            pressure_hpa,
            std::f64::consts::FRAC_PI_6,
        ),
        WindConditions::default(),
        air,
        presence,
    );
    solver.set_max_range(30_000.0);
    let lob = solver.solve().unwrap();
    (
        last.drop_m,
        last.velocity_mps,
        lob.max_range,
        lob.max_height,
        lob.impact_velocity,
    )
}

#[test]
fn the_aloft_vapour_model_shifts_are_engine_output() {
    let section = unreleased();
    for (temperature_c, humidity_percent, height_m, density, sound_speed) in AIR_BEFORE_S8 {
        assert!(ALOFT_QUOTED.contains(&(temperature_c, humidity_percent)));
        let air = air_over_sea_level(temperature_c, humidity_percent, height_m);
        let what = format!("air at {temperature_c} C / {humidity_percent} % +{height_m} m");
        assert_contains(
            &section,
            &format!(
                "{sound_speed:.3} → {:.3} m/s ({:+.3}%)",
                air.speed_of_sound_mps,
                (air.speed_of_sound_mps / sound_speed - 1.0) * 100.0
            ),
            &what,
        );
        assert_contains(
            &section,
            &format!(
                "{density:.4} → {:.4} kg/m³ ({:+.3}%)",
                air.density_kg_m3,
                (air.density_kg_m3 / density - 1.0) * 100.0
            ),
            &what,
        );
    }
    for (
        station,
        altitude_m,
        pressure_hpa,
        temperature_c,
        humidity_percent,
        drop,
        speed,
        range,
        apex,
        impact,
    ) in SHOTS_BEFORE_S8
    {
        let (drop1, speed1, range1, apex1, impact1) =
            solver_shots(altitude_m, pressure_hpa, temperature_c, humidity_percent);
        let what = format!("{station} {temperature_c} C / {humidity_percent} %");
        assert_contains(
            &section,
            &format!("{:+.3} mm", (drop1 - drop) * 1000.0),
            &what,
        );
        assert_contains(&section, &format!("{:+.4} m/s", speed1 - speed), &what);
        assert_contains(&section, &format!("{:+.3} m (", range1 - range), &what);
        assert_contains(&section, &format!("{range:.2} → {range1:.2} m"), &what);
        if station == "sea level" {
            assert_contains(&section, &format!("{:+.3} m |", apex1 - apex), &what);
            assert_contains(&section, &format!("{:+.3} m/s |", impact1 - impact), &what);
        }
    }
}

// =================================================================================================
// The muzzle Miller Sg: dry by design
// =================================================================================================

/// The factor a humid density correction would put on the muzzle Miller Sg at `(T, RH)`: the
/// dry CIPM-2007 density over the moist one at the same temperature and pressure.
fn humid_miller_factor(temperature_c: f64, humidity_percent: f64) -> f64 {
    calculate_air_density_cimp(temperature_c, SEA_LEVEL_HPA, 0.0)
        / calculate_air_density_cimp(temperature_c, SEA_LEVEL_HPA, humidity_percent)
}

#[test]
fn the_miller_exception_is_documented_with_engine_numbers() {
    let stability = read("src/stability.rs");
    let section = unreleased();
    let cli_usage = read("CLI_USAGE.md");
    for (temperature_c, humidity_percent) in QUOTED {
        let factor = format!(
            "{:.4}",
            humid_miller_factor(temperature_c, humidity_percent)
        );
        let what = format!("the Miller factor at {temperature_c} C / {humidity_percent} %");
        assert_contains(&stability, &factor, &format!("src/stability.rs: {what}"));
        assert_contains(&section, &factor, &format!("CHANGELOG: {what}"));
        assert_contains(&cli_usage, &factor, &format!("CLI_USAGE.md: {what}"));
    }
}

/// Pin (passes on the code as it stands): the solve's muzzle Sg is the published dry Miller
/// correction, so it is the same at every humidity, in the band included.
#[test]
fn the_solves_muzzle_sg_does_not_move_with_humidity() {
    let request = read("tests/fixtures/solve_json_v1/representative.request.json");
    let sg_at = |humidity: f64| {
        let mut json: serde_json::Value = serde_json::from_str(&request).unwrap();
        json["atmosphere"]["relative_humidity"] = humidity.into();
        let decoded = decode_solve_request_v1(&json.to_string()).expect("valid request");
        let success = serde_json::to_value(solve_v1(decoded).expect("solves")).unwrap();
        success["summary"]["stability_factor"]
            .as_f64()
            .expect("the solve reports its stability factor")
    };
    let dry = sg_at(0.0);
    for humidity in [0.5, 1.0, 1.03] {
        assert_eq!(sg_at(humidity).to_bits(), dry.to_bits(), "RH {humidity}");
    }
}

// =================================================================================================
// Model descriptions and stale wording (H19, H20)
// =================================================================================================

#[test]
fn cli_usage_describes_the_humid_air_model_the_engine_uses() {
    let cli_usage = read("CLI_USAGE.md");
    assert!(
        !cli_usage.contains("Arden Buck"),
        "the engine's vapour pressure is the Wagner-Pruss saturation equation, not Arden Buck"
    );
    for needle in [
        "Wagner",
        "CIPM-2007",
        // The first-order moist-air coefficients moist_speed_of_sound uses.
        "0.062",
        "0.378",
        "#### Humidity: units, defaults and validation",
    ] {
        assert_contains(&cli_usage, needle, "CLI_USAGE.md humid-air model");
    }
}

#[test]
fn the_saturation_pressure_is_labelled_by_its_source() {
    // The coefficients are the Wagner-Pruss saturation-pressure equation (the IAPWS
    // supplementary release on saturation properties); IAPWS-IF97's region-4 equation is a
    // different formula.
    let atmosphere = read("src/atmosphere.rs");
    assert!(!atmosphere.contains("IAPWS-IF97"), "p_sv is not IAPWS-IF97");
    assert_contains(&atmosphere, "Wagner", "src/atmosphere.rs p_sv label");
}

#[test]
fn the_sound_speed_constant_doc_gives_the_models_dry_value() {
    let constants = read("src/constants.rs");
    assert!(
        !constants.contains("331.3"),
        "the engine has no 331.3 sqrt(T/273.15) model"
    );
    let dry_15c = moist_speed_of_sound(288.15, 101_325.0, 0.0);
    assert_contains(&constants, &format!("{dry_15c:.3}"), "src/constants.rs");
}

#[test]
fn the_migration_notes_and_readme_cover_the_kernel_atmosphere() {
    let migration = read("MCCOY_MIGRATION.md");
    for needle in ["KernelAtmosphere", "humidity_percent", "from_legacy_params"] {
        assert_contains(&migration, needle, "MCCOY_MIGRATION.md");
    }
    let readme = read("README.md");
    for needle in ["RelativeHumidity", "KernelAtmosphere", "relative_humidity"] {
        assert_contains(&readme, needle, "README.md");
    }
}

#[test]
fn each_migration_section_keeps_its_own_bullets() {
    // The 0.50.0 humidity notes were added after the SI Unit Unification notes. The SI
    // section's bug-fix list must stay whole, and none of its bullets may read as a 0.50.0
    // humidity change.
    let migration = read("MCCOY_MIGRATION.md");
    let si_start = migration
        .find("\n# SI Unit Unification")
        .expect("MCCOY_MIGRATION.md has the SI Unit Unification notes");
    let humidity_start = migration
        .find("\n# Humidity Consistency (0.50.0)")
        .expect("MCCOY_MIGRATION.md has the 0.50.0 humidity notes");
    assert!(
        si_start < humidity_start,
        "the 0.50.0 humidity notes follow the SI notes"
    );
    let si = &migration[si_start..humidity_start];
    let humidity_body = &migration[humidity_start + 1..];
    let humidity = &humidity_body[..humidity_body.find("\n# ").unwrap_or(humidity_body.len())];
    for bullet in [
        "- **Enhanced spin drift**",
        "- **cluster_bc**",
        "- **Monte Carlo trajectory**",
        "- **Reynolds correction**",
    ] {
        assert_contains(
            si,
            bullet,
            "MCCOY_MIGRATION.md SI Unit Unification bug fixes",
        );
        assert!(
            !humidity.contains(bullet),
            "MCCOY_MIGRATION.md: the SI section's `{bullet}` bullet sits in the 0.50.0 humidity notes"
        );
    }
}

/// `ballistics <path> --help`, which must succeed.
#[cfg(feature = "cli")]
fn help_of(path: &[String]) -> String {
    let output = std::process::Command::new(env!("CARGO_BIN_EXE_ballistics"))
        .args(path)
        .arg("--help")
        .output()
        .expect("spawn ballistics");
    assert!(
        output.status.success(),
        "`ballistics {} --help` failed: {}",
        path.join(" "),
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).expect("help is UTF-8")
}

/// The subcommand names a clap help lists under `Commands:`, without `help`.
#[cfg(feature = "cli")]
fn subcommands_in(help: &str) -> Vec<String> {
    help.lines()
        .skip_while(|line| line.trim_end() != "Commands:")
        .skip(1)
        .take_while(|line| !line.trim().is_empty())
        .filter(|line| line.starts_with("  ") && !line.starts_with("   "))
        .filter_map(|line| line.split_whitespace().next())
        .filter(|name| *name != "help")
        .map(str::to_string)
        .collect()
}

/// Whether a clap help declares the long option `--<flag> <VALUE>` (not a prose mention).
#[cfg(feature = "cli")]
fn declares_option(help: &str, flag: &str) -> bool {
    let spec = format!("--{flag} <");
    help.lines().map(str::trim_start).any(|line| {
        line.starts_with(&spec) || (line.starts_with('-') && line.contains(&format!(", {spec}")))
    })
}

/// The backticked names in the first cell of the humidity table's row for `field`.
#[cfg(feature = "cli")]
fn humidity_table_commands(field: &str) -> std::collections::BTreeSet<String> {
    let cli_usage = read("CLI_USAGE.md");
    let heading = "#### Humidity: units, defaults and validation";
    let start = cli_usage
        .find(heading)
        .expect("CLI_USAGE.md has the humidity section");
    let section = &cli_usage[start + heading.len()..];
    let section = &section[..section.find("\n#").unwrap_or(section.len())];
    let row = section
        .lines()
        .filter(|line| line.starts_with('|'))
        .find(|line| line.split('|').nth(2).map(str::trim) == Some(field))
        .unwrap_or_else(|| panic!("CLI_USAGE.md humidity table has no {field} row"));
    let where_cell = row.split('|').nth(1).expect("the row has a first cell");
    where_cell
        .split('`')
        .skip(1)
        .step_by(2)
        // `trajectory --auto-zero` names the command `trajectory`.
        .map(|name| name.split(" --").next().unwrap_or(name).to_string())
        .collect()
}

#[cfg(feature = "cli")]
#[test]
fn the_cli_usage_humidity_table_names_every_command_that_takes_a_humidity_flag() {
    // Walk the binary's whole command tree, so a command that takes --humidity cannot be left
    // out of the table (nor can the table name one that does not take it).
    let mut humidity = std::collections::BTreeSet::new();
    let mut zero_humidity = std::collections::BTreeSet::new();
    let mut pending: Vec<Vec<String>> = vec![Vec::new()];
    while let Some(path) = pending.pop() {
        let help = help_of(&path);
        if !path.is_empty() {
            if declares_option(&help, "humidity") {
                humidity.insert(path.join(" "));
            }
            if declares_option(&help, "zero-humidity") {
                zero_humidity.insert(path.join(" "));
            }
        }
        for name in subcommands_in(&help) {
            let mut child = path.clone();
            child.push(name);
            pending.push(child);
        }
    }
    // The walk itself must find the commands known to take the flag.
    for known in ["trajectory", "monte-carlo", "reticle hold", "profile save"] {
        assert!(
            humidity.contains(known),
            "the help walk missed `{known}`: {humidity:?}"
        );
    }
    assert_eq!(
        humidity_table_commands("`--humidity`"),
        humidity,
        "CLI_USAGE.md's humidity table must name exactly the commands that take --humidity"
    );
    assert_eq!(
        humidity_table_commands("`--zero-humidity`"),
        zero_humidity,
        "CLI_USAGE.md's humidity table must name exactly the commands that take --zero-humidity"
    );
}

#[test]
fn stale_humidity_comments_are_gone() {
    for (path, stale, why) in [
        (
            "src/main.rs",
            "and `fast_integrate` overwrites the inputs'",
            "run_trajectory never calls fast_integrate",
        ),
        (
            "src/cli_api.rs",
            "reports a per-point Cd derived from these same two values",
            "the per-point Cd is evaluated at the station sound speed, not this step's",
        ),
        (
            "src/pdf_dope_card.rs",
            "the solver's air density comes from",
            "every solver's air comes from air_at, not only one solver's",
        ),
        (
            "src/atmosphere.rs",
            "Every solver (cli_api / monte_carlo / ffi / fast_trajectory) reaches this one formula",
            "zones reach the CIPM density directly, stations through StationAir::new",
        ),
    ] {
        assert!(!read(path).contains(stale), "{path}: stale comment ({why})");
    }
}
