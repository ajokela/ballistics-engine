//! Byte-identity pins for the humidity-consistency train (0.50.0), step S1 (plan step C0a).
//!
//! # Why this file exists
//!
//! The humidity train routes every solver family through one air lookup and gives the kernel
//! API a type that cannot describe an atmosphere without its relative humidity. Its Phase 1
//! promise is that every `TrajectorySolver`-backed surface stays BYTE-IDENTICAL at every
//! valid RH, and so do direct ("uniform") mode on every family, station density and sound
//! speed, the muzzle Miller Sg, Litz spin drift, aerodynamic jump and density altitude. These
//! pins hold the train to that promise with `f64::to_bits` equality. A refactor that moves
//! one of them has changed physics, whatever its description says.
//!
//! # Rules for later steps
//!
//! - A pin is never edited. It may only be MIGRATED: when an API it calls is renamed or
//!   retyped (for example the kernel atmosphere tuple becoming `KernelAtmosphere`), the call
//!   changes and the pinned constant does not.
//! - Each group's header says which later step may change it. Two words are used:
//!   - `IMMUTABLE`: must pass unchanged through every step, S8 included.
//!   - `MAY CHANGE ONLY IN S8`: the unified aloft vapour model (plan C7) may re-pin it, and
//!     no other step may. Every such group involves a relative humidity above zero, at the
//!     station or in a zone; at RH 0 the plan claims S8 is bit-identical, so the RH-0 cases
//!     live in separate `IMMUTABLE` groups.
//! - Kernel outputs at RH > 0 (`compute_derivatives`, `try_integrate_trajectory`,
//!   `try_solve_trajectory_rust`, `fast_integrate*`, the lean Monte Carlo) are deliberately NOT
//!   pinned: today they fly the dry speed of sound at any RH, and S3 fixes that.
//! - Kernel STANDARD mode at RH 0 is pinned to a relative 1e-15, not bits: the plan allows
//!   last-ulp movement there (fast_integrate loses a `ratio * 1.225 / 1.225` round trip). For
//!   scale, a one-ulp change in the base density ratio was measured to move these kernels by
//!   at most 5.4e-16 relative, so the tolerance absorbs about two such ulps and no more.
//!
//! In-crate pins for the `pub(crate)` gravity-aware entry points live in
//! `src/humidity_byte_identity_tests.rs` under the same rules.
//!
//! # How the pinned values were generated
//!
//! By running this file's own probes against the engine at 261c1eb (ballistics-engine
//! 0.49.0, the commit this train branched from) and pasting the printed constants verbatim:
//!
//! ```text
//! cargo test --test humidity_byte_identity -- --ignored --nocapture print_current_pin_values
//! ```
//!
//! They were captured on macOS (aarch64), and the identical file then passed in Linux glibc
//! containers on aarch64 and x86_64 (the CI platform).
//!
//! Each value is the IEEE-754 bit pattern of an `f64` (`f64::to_bits`), or an exact count.
//! A `*.is_some` entry is 1 or 0. A `text` entry is the FNV-1a hash of a string such as a
//! termination reason. A `digest` entry is a 64-bit FNV-1a hash over the bit patterns of a
//! whole trajectory, so a pin covers every point and not just the sampled ones.
//!
//! # Why the digests cover time and position only
//!
//! The flight loops call `powf`/`exp` (the hydrostatic lapse, the Litz drift `t^1.83`), and
//! macOS's libm and glibc round a small fraction of those calls to different last bits. Such a
//! one-ulp difference in an acceleration moves the velocity's last bit for a few steps and is
//! then rounded away, so time and position stay bit-identical across the two libms while the
//! per-point speed, energy and Cd do not, at some points of a long flight. The first
//! version of this file digested speed too and failed on Linux for exactly that reason. Speed,
//! energy, Mach and Cd are therefore pinned at explicit indices, every one of which agrees on
//! all three platforms, and the FFI case leaves Litz drift off (a per-point `powf` added to
//! the lateral position); `TS_EFFECTS_*` and `STABILITY_SPIN_JUMP` pin the drift instead.

use ballistics_engine::atmosphere::{
    calculate_air_density_cimp, calculate_atmosphere, density_matched_altitude_m,
    faa_rule_density_altitude_ft, get_local_atmosphere_humid, moist_speed_of_sound,
    nws_pressure_altitude_ft, resolve_atmosphere_for_density_altitude, AtmoSock,
    StationAtmospherePresence,
};
use ballistics_engine::derivatives::compute_derivatives;
use ballistics_engine::fast_trajectory::{
    aerodynamic_jump_launch_offset_rad, fast_integrate, fast_integrate_with_segments,
    FastIntegrationParams, FastSolution,
};
use ballistics_engine::monte_carlo::solve_trajectory_for_monte_carlo_with_station_presence;
use ballistics_engine::spin_drift::{effective_sg_from_inputs, litz_drift_meters};
use ballistics_engine::stability::compute_stability_coefficient;
use ballistics_engine::trajectory_integration::{
    try_integrate_trajectory, try_solve_trajectory_rust, TrajectoryParams,
};
use ballistics_engine::wind::{WindSegment, WindSock};
use ballistics_engine::{
    calculate_zero_angle_with_station_presence, run_monte_carlo_with_station_presence,
    run_monte_carlo_with_wind_and_direction_std_dev_seeded, AtmosphericConditions, BallisticInputs,
    DragModel, MonteCarloParams, MonteCarloResults, TrajectoryResult, TrajectorySolver,
    WindConditions, DEFAULT_HIT_RADIUS_M,
};
use nalgebra::Vector3;
use serde_json::{json, Value};

// =================================================================================================
// Pin machinery
// =================================================================================================

/// One pinned quantity: the bits of an `f64`, or an exact integer (a count, a flag, a digest).
#[derive(Clone, Copy)]
enum Pin {
    F(f64),
    N(u64),
}

impl Pin {
    fn bits(self) -> u64 {
        match self {
            Pin::F(value) => value.to_bits(),
            Pin::N(value) => value,
        }
    }
}

const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

fn fnv1a(bytes: impl IntoIterator<Item = u8>) -> u64 {
    bytes.into_iter().fold(FNV_OFFSET, |hash, byte| {
        (hash ^ u64::from(byte)).wrapping_mul(FNV_PRIME)
    })
}

/// FNV-1a over the little-endian bit patterns of a sequence of `f64`s.
fn digest_f64s(values: impl IntoIterator<Item = f64>) -> u64 {
    fnv1a(
        values
            .into_iter()
            .flat_map(|value| value.to_bits().to_le_bytes()),
    )
}

/// The pinned quantities of one case, each with a label used in failure messages and in the
/// generator's output.
#[derive(Default)]
struct Pins {
    entries: Vec<(String, Pin)>,
}

impl Pins {
    fn f(&mut self, label: impl Into<String>, value: f64) {
        self.entries.push((label.into(), Pin::F(value)));
    }

    fn n(&mut self, label: impl Into<String>, value: u64) {
        self.entries.push((label.into(), Pin::N(value)));
    }

    /// An optional value as two entries: presence (1/0), then the value (0.0 when absent).
    fn opt(&mut self, label: &str, value: Option<f64>) {
        self.n(format!("{label}.is_some"), u64::from(value.is_some()));
        self.f(label, value.unwrap_or(0.0));
    }

    fn text(&mut self, label: &str, text: &str) {
        self.n(format!("{label} text {text:?}"), fnv1a(text.bytes()));
    }

    fn digest(&mut self, label: &str, values: impl IntoIterator<Item = f64>) {
        self.n(format!("{label} digest"), digest_f64s(values));
    }

    fn bits(&self) -> Vec<u64> {
        self.entries.iter().map(|(_, pin)| pin.bits()).collect()
    }
}

/// A group of named cases, in a fixed order.
type Probe = Vec<(String, Pins)>;
/// The pinned form of a [`Probe`].
type Pinned = &'static [(&'static str, &'static [u64])];

fn describe(pin: Pin) -> String {
    match pin {
        Pin::F(value) => format!("{value:?}"),
        Pin::N(value) => format!("{value}"),
    }
}

fn case_shape_or_panic(group: &str, actual: &Probe, expected: Pinned) {
    let actual_names: Vec<&str> = actual.iter().map(|(name, _)| name.as_str()).collect();
    let expected_names: Vec<&str> = expected.iter().map(|(name, _)| *name).collect();
    assert_eq!(
        actual_names, expected_names,
        "{group}: the probe's cases no longer match the pinned cases"
    );
    for ((name, pins), (_, bits)) in actual.iter().zip(expected) {
        assert_eq!(
            pins.entries.len(),
            bits.len(),
            "{group}/{name}: the probe records {} quantities, the pin holds {}",
            pins.entries.len(),
            bits.len()
        );
    }
}

/// Every quantity of every case must equal its pin bit for bit.
fn assert_bits(group: &str, actual: &Probe, expected: Pinned) {
    case_shape_or_panic(group, actual, expected);
    let mut moved = Vec::new();
    let mut total = 0;
    for ((name, pins), (_, bits)) in actual.iter().zip(expected) {
        for ((label, pin), &want) in pins.entries.iter().zip(bits.iter()) {
            total += 1;
            if pin.bits() != want {
                let pinned = match pin {
                    Pin::F(_) => format!("{:?}", f64::from_bits(want)),
                    Pin::N(_) => format!("{want}"),
                };
                moved.push(format!(
                    "  {name} :: {label}: pinned {want:#018x} ({pinned}), now {:#018x} ({})",
                    pin.bits(),
                    describe(*pin)
                ));
            }
        }
    }
    assert!(
        moved.is_empty(),
        "{group}: {} of {total} byte-identity pins moved:\n{}",
        moved.len(),
        moved.join("\n")
    );
}

/// Counts must match exactly; each float must match within `rel_tol` of its pin's magnitude
/// (a pinned 0.0 must stay exactly 0.0).
fn assert_relative(group: &str, actual: &Probe, expected: Pinned, rel_tol: f64) {
    case_shape_or_panic(group, actual, expected);
    let mut moved = Vec::new();
    let mut total = 0;
    for ((name, pins), (_, bits)) in actual.iter().zip(expected) {
        for ((label, pin), &want) in pins.entries.iter().zip(bits.iter()) {
            total += 1;
            match *pin {
                Pin::N(value) => {
                    if value != want {
                        moved.push(format!("  {name} :: {label}: pinned {want}, now {value}"));
                    }
                }
                Pin::F(value) => {
                    let pinned = f64::from_bits(want);
                    let ok =
                        value.to_bits() == want || (value - pinned).abs() <= rel_tol * pinned.abs();
                    if !ok {
                        moved.push(format!(
                            "  {name} :: {label}: pinned {pinned:?}, now {value:?} \
                             (relative {:e}, allowed {rel_tol:e})",
                            ((value - pinned) / pinned).abs()
                        ));
                    }
                }
            }
        }
    }
    assert!(
        moved.is_empty(),
        "{group}: {} of {total} pins moved beyond {rel_tol:e} relative:\n{}",
        moved.len(),
        moved.join("\n")
    );
}

/// Print a probe as the Rust constant this file pins it with.
fn print_probe(const_name: &str, probe: &Probe) {
    println!("const {const_name}: Pinned = &[");
    for (name, pins) in probe {
        println!("    (");
        println!("        {name:?},");
        println!("        &[");
        for chunk in pins.bits().chunks(4) {
            let line: Vec<String> = chunk.iter().map(|bits| format!("{bits:#018x}")).collect();
            println!("            {},", line.join(", "));
        }
        println!("        ],");
        println!("    ),");
    }
    println!("];");
    println!();
}

/// How a group's floats are compared to their pins.
#[derive(Clone, Copy)]
enum Compare {
    Bits,
    /// Relative tolerance for kernel standard mode at RH 0.
    Relative,
}

const KERNEL_STANDARD_RH0_REL_TOL: f64 = 1e-15;

// =================================================================================================
// Shared fixture: a .308 175 gr-class load (G7 0.243, 823 m/s, 1:10")
// =================================================================================================

const MASS_KG: f64 = 0.01134;
const DIAMETER_M: f64 = 0.00782;
const LENGTH_M: f64 = 0.031;
const MUZZLE_VELOCITY_MPS: f64 = 823.0;
const BC_G7: f64 = 0.243;
const TWIST_IN: f64 = 10.0;
const FLAT_MUZZLE_ANGLE_RAD: f64 = 0.002;
/// 30 degrees: a high-angle lob that climbs more than a kilometre above the station, so the
/// aloft (lapse) half of the air model carries real weight.
const LOB_MUZZLE_ANGLE_RAD: f64 = std::f64::consts::FRAC_PI_6;
/// 4 m/s from the right, so lateral position is never trivially zero.
const CROSSWIND_MPS: f64 = 4.0;
const CROSSWIND_FROM_RAD: f64 = std::f64::consts::FRAC_PI_2;

#[derive(Clone, Copy)]
enum Method {
    Euler,
    Rk4,
    Rk45,
}

impl Method {
    const ALL: [Method; 3] = [Method::Euler, Method::Rk4, Method::Rk45];

    fn name(self) -> &'static str {
        match self {
            Method::Euler => "euler",
            Method::Rk4 => "rk4",
            Method::Rk45 => "rk45",
        }
    }

    fn flags(self) -> (bool, bool) {
        match self {
            Method::Euler => (false, false),
            Method::Rk4 => (true, false),
            Method::Rk45 => (true, true),
        }
    }
}

/// A station the solver is constructed for.
#[derive(Clone, Copy)]
struct Station {
    name: &'static str,
    altitude_m: f64,
    temperature_c: f64,
    pressure_hpa: f64,
    /// `true`: `TrajectorySolver::new_with_station_presence(.., BOTH)`, the presence-aware
    /// constructor; `false`: the same with `NEITHER`, which is exactly what the presence-blind
    /// `TrajectorySolver::new` calls.
    presence_aware: bool,
}

impl Station {
    fn presence(self) -> StationAtmospherePresence {
        if self.presence_aware {
            StationAtmospherePresence::BOTH
        } else {
            StationAtmospherePresence::NEITHER
        }
    }
}

/// ICAO standard sea level through the presence-blind constructor.
const SEA_LEVEL: Station = Station {
    name: "icao-sea-level",
    altitude_m: 0.0,
    temperature_c: 15.0,
    pressure_hpa: 1013.25,
    presence_aware: false,
};

/// 1500 m with a temperature and pressure far from ICAO's 5.25 C / 845.6 hPa there, through
/// the presence-aware constructor (both supplied, so both fly as given).
const HIGH_STATION: Station = Station {
    name: "1500m-25C-830hPa",
    altitude_m: 1500.0,
    temperature_c: 25.0,
    pressure_hpa: 830.0,
    presence_aware: true,
};

const STATIONS: [Station; 2] = [SEA_LEVEL, HIGH_STATION];

/// Downrange atmosphere zones for a station.
#[derive(Clone, Copy, PartialEq)]
enum Zones {
    None,
    /// Two zones whose RH (20 %, 80 %) differs from each other and from every station RH.
    Humid,
    /// The same two zones, both at RH 0.
    Dry,
}

impl Zones {
    fn name(self) -> &'static str {
        match self {
            Zones::None => "no-zones",
            Zones::Humid => "zones-rh20-rh80",
            Zones::Dry => "zones-rh0-rh0",
        }
    }

    fn segments(self, station: Station) -> Option<Vec<(f64, f64, f64, f64)>> {
        let (near_rh, far_rh) = match self {
            Zones::None => return None,
            Zones::Humid => (20.0, 80.0),
            Zones::Dry => (0.0, 0.0),
        };
        Some(vec![
            (
                station.temperature_c + 2.0,
                station.pressure_hpa - 3.0,
                near_rh,
                400.0,
            ),
            (
                station.temperature_c - 3.0,
                station.pressure_hpa + 2.0,
                far_rh,
                100_000.0,
            ),
        ])
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Shot {
    /// Near-flat fire to 1000 m with sampling on.
    Flat,
    /// The 30-degree lob, until it falls 100 m below the muzzle (the default ground threshold).
    Lob,
}

impl Shot {
    fn name(self) -> &'static str {
        match self {
            Shot::Flat => "flat",
            Shot::Lob => "lob30",
        }
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
        station.presence(),
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

fn probe_ts_flat_dry() -> Probe {
    ts_grid(Shot::Flat, false)
}

fn probe_ts_flat_humid() -> Probe {
    ts_grid(Shot::Flat, true)
}

fn probe_ts_lob_dry() -> Probe {
    ts_grid(Shot::Lob, false)
}

fn probe_ts_lob_humid() -> Probe {
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

fn probe_ts_effects_dry() -> Probe {
    ts_effects(false)
}

fn probe_ts_effects_humid() -> Probe {
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
                    station.presence(),
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

fn probe_zero_dry() -> Probe {
    zero_grid(false)
}

fn probe_zero_humid() -> Probe {
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

fn probe_monte_carlo_default_humidity() -> Probe {
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

fn probe_station_air() -> Probe {
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

fn probe_local_air_dry() -> Probe {
    local_air(false)
}

fn probe_local_air_humid() -> Probe {
    local_air(true)
}

fn probe_density_altitude() -> Probe {
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
        let request = ballistics_engine::atmosphere_service::decode_density_altitude_request_v1(
            &json!({ "atmosphere": atmosphere }),
        )
        .expect("decode density-altitude request");
        let response = ballistics_engine::atmosphere_service::density_altitude_v1(&request)
            .expect("density altitude");
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
fn probe_stability_spin_jump() -> Probe {
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
            aerodynamic_jump_launch_offset_rad(&jump_inputs, atmo),
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

fn trajectory_params(
    atmos_params: (f64, f64, f64, f64),
    atmo_sock: Option<AtmoSock>,
) -> TrajectoryParams {
    TrajectoryParams {
        mass_kg: MASS_KG,
        bc: BC_G7,
        drag_model: DragModel::G7,
        wind_segments: crosswind_segments(),
        atmos_params,
        omega_vector: None,
        enable_spin_drift: false,
        enable_magnus: false,
        enable_coriolis: false,
        target_distance_m: 800.0,
        enable_wind_shear: false,
        wind_shear_model: "none".to_string(),
        shooter_altitude_m: 0.0,
        is_twist_right: true,
        shooting_angle: 0.0,
        bullet_diameter: DIAMETER_M,
        bullet_length: LENGTH_M,
        twist_rate: TWIST_IN,
        custom_drag_table: None,
        cd_scale: 1.0,
        custom_drag_reference: false,
        bc_segments: None,
        use_bc_segments: false,
        bc_segments_data: None,
        ground_threshold: -1000.0,
        atmo_sock,
    }
}

fn fast_params(
    atmo_params: (f64, f64, f64, f64),
    atmo_sock: Option<AtmoSock>,
) -> FastIntegrationParams {
    FastIntegrationParams {
        horiz: 800.0,
        vert: 0.0,
        initial_state: launch_state(),
        t_span: (0.0, 30.0),
        atmo_params,
        atmo_sock,
    }
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
        let derivative = compute_derivatives(
            position,
            Vector3::new(700.0, 2.0, 0.5),
            &inputs,
            Vector3::new(0.0, 0.0, -CROSSWIND_MPS),
            atmos_params,
            BC_G7,
            None,
            0.0,
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
            trajectory_params(atmos_params, sock()),
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
        let points = try_solve_trajectory_rust(
            launch_state(),
            (0.0, 3.0),
            MASS_KG,
            BC_G7,
            DragModel::G7,
            crosswind_segments(),
            atmos_params,
            None,
            false,
            false,
            false,
            "RK45".to_string(),
            1e-6,
            0.01,
            800.0,
        )
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
        &fast_integrate(&inputs, &wind_sock, fast_params(atmos_params, sock())),
        digest,
    );
    cases.push((format!("{prefix}/fast_integrate"), pins));

    let mut pins = Pins::default();
    record_fast_solution(
        &mut pins,
        &fast_integrate_with_segments(
            &inputs,
            crosswind_segments(),
            fast_params(atmos_params, sock()),
        ),
        digest,
    );
    cases.push((format!("{prefix}/fast_integrate_with_segments"), pins));

    cases
}

/// Direct mode: no station, no humidity, so every family is IMMUTABLE bit for bit.
fn probe_kernel_direct() -> Probe {
    kernel_cases("direct(1.1,335)", SEA_LEVEL, DIRECT, Zones::None, true)
}

/// Standard mode at RH 0 (station and zones dry), plus the lean Monte Carlo at RH 0.
fn probe_kernel_standard_rh0() -> Probe {
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
        let inputs = BallisticInputs {
            humidity: 0.0,
            target_distance: 800.0,
            ..kernel_inputs(station)
        };
        let output =
            solve_trajectory_for_monte_carlo_with_station_presence(&inputs, station.presence())
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
        let decoded = ballistics_engine::decode_solve_request_v1(&request.to_string())
            .expect("valid solve-json request");
        let success = ballistics_engine::solve_v1(decoded).expect("solve_v1 succeeds");
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
        let raw = ballistics_engine::bridge::bridge_call(&envelope.to_string());
        let out: Value = serde_json::from_str(&raw).expect("bridge output parses");
        assert_eq!(out["ok"], true, "bridge solve failed: {out}");
        let mut pins = Pins::default();
        record_solve_json(&mut pins, &out["result"]);
        cases.push((format!("bridge-solve/{name}"), pins));
    }
    cases
}

#[cfg(feature = "ffi")]
fn front_end_ffi(humid: bool) -> Probe {
    use ballistics_engine::ffi::{
        ballistics_calculate_trajectory, ballistics_free_trajectory_result,
        FFIAtmosphericConditions, FFIBallisticInputs, FFIWindConditions,
    };

    let mut cases = Vec::new();
    let rhs: &[f64] = if humid { &[50.0, 100.0] } else { &[0.0] };
    for station in STATIONS {
        for &rh in rhs {
            let inputs = FFIBallisticInputs {
                muzzle_velocity: MUZZLE_VELOCITY_MPS,
                muzzle_angle: FLAT_MUZZLE_ANGLE_RAD,
                bc_value: BC_G7,
                bullet_mass: MASS_KG,
                bullet_diameter: DIAMETER_M,
                bc_type: 1, // G7
                sight_height: 0.05,
                target_distance: 600.0,
                temperature: station.temperature_c,
                twist_rate: TWIST_IN,
                is_twist_right: 1,
                shooting_angle: 0.0,
                altitude: station.altitude_m,
                latitude: f64::NAN,
                azimuth_angle: 0.0,
                use_rk4: 1,
                use_adaptive_rk45: 0,
                enable_wind_shear: 0,
                enable_trajectory_sampling: 1,
                sample_interval: 50.0,
                enable_pitch_damping: 0,
                enable_precession_nutation: 0,
                // Off: Litz drift is a per-point powf, whose last bit differs between libms
                // (see the module docs). TS_EFFECTS_* and STABILITY_SPIN_JUMP pin it.
                enable_spin_drift: 0,
                enable_magnus: 0,
                enable_coriolis: 0,
                shot_azimuth: 0.0,
                cant_angle: 0.0,
                zero_poi_vertical: 0.0,
                zero_poi_horizontal: 0.0,
                sight_offset_lateral: 0.0,
            };
            let wind = FFIWindConditions {
                speed: CROSSWIND_MPS,
                direction: CROSSWIND_FROM_RAD,
                vertical_speed: 0.0,
            };
            let atmosphere = FFIAtmosphericConditions {
                temperature: station.temperature_c,
                pressure: station.pressure_hpa,
                humidity: rh,
                altitude: station.altitude_m,
            };
            let name = format!("ffi-trajectory/rk4/{}/rh{rh}", station.name);
            let mut pins = Pins::default();
            // SAFETY: every pointer is to a live local, and the result is freed exactly once.
            unsafe {
                let result =
                    ballistics_calculate_trajectory(&inputs, &wind, &atmosphere, 600.0, 1.0);
                assert!(!result.is_null(), "{name}: FFI solve returned null");
                let r = &*result;
                pins.f("max_range", r.max_range);
                pins.f("max_height", r.max_height);
                pins.f("time_of_flight", r.time_of_flight);
                pins.f("impact_velocity", r.impact_velocity);
                pins.f("impact_energy", r.impact_energy);
                pins.f("transonic_mach (NaN when absent)", r.transonic_mach);
                let points = std::slice::from_raw_parts(r.points, r.point_count as usize);
                let samples =
                    std::slice::from_raw_parts(r.sampled_points, r.sampled_point_count as usize);
                pins.n("point_count", points.len() as u64);
                let last = points.last().expect("points");
                pins.f("points[last].time", last.time);
                pins.f("points[last].position_x", last.position_x);
                pins.f("points[last].position_y", last.position_y);
                pins.f("points[last].position_z", last.position_z);
                pins.f("points[last].velocity_magnitude", last.velocity_magnitude);
                pins.digest(
                    "points (time, x, y, z)",
                    points.iter().flat_map(|point| {
                        [
                            point.time,
                            point.position_x,
                            point.position_y,
                            point.position_z,
                        ]
                    }),
                );
                pins.n("sampled_point_count", samples.len() as u64);
                pins.digest(
                    "sampled_points (distance, time, drop, windage)",
                    samples.iter().flat_map(|sample| {
                        [
                            sample.distance,
                            sample.time,
                            sample.drop_meters,
                            sample.windage_meters,
                        ]
                    }),
                );
                let middle = &samples[samples.len() / 2];
                pins.f("sampled_points[mid].velocity_mps", middle.velocity_mps);
                pins.f("sampled_points[mid].mach", middle.mach);
                ballistics_free_trajectory_result(result);
            }
            cases.push((name, pins));
        }
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
            let request: ballistics_engine::card_service::CardRequestV1 =
                serde_json::from_value(json!({
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
            let response =
                ballistics_engine::card_service::come_ups_v1(&request).expect("come-ups card");
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

fn probe_front_ends_dry() -> Probe {
    let mut cases = front_end_solve_json(false);
    cases.extend(front_end_card(false));
    cases
}

fn probe_front_ends_humid() -> Probe {
    let mut cases = front_end_solve_json(true);
    cases.extend(front_end_card(true));
    cases
}

#[cfg(feature = "bridge")]
fn probe_bridge_solve_dry() -> Probe {
    front_end_bridge(false)
}

#[cfg(feature = "bridge")]
fn probe_bridge_solve_humid() -> Probe {
    front_end_bridge(true)
}

#[cfg(feature = "ffi")]
fn probe_ffi_trajectory_dry() -> Probe {
    front_end_ffi(false)
}

#[cfg(feature = "ffi")]
fn probe_ffi_trajectory_humid() -> Probe {
    front_end_ffi(true)
}

// =================================================================================================
// The groups, their pins, and the tests
// =================================================================================================

/// One pinned group, as the generator prints it.
struct Group {
    /// The constant the group is pinned in.
    name: &'static str,
    /// The doc comment printed above the constant: what is pinned, and which later step may
    /// change it.
    header: &'static str,
    probe: fn() -> Probe,
    compare: Compare,
    /// The cargo feature the group needs, if any.
    feature: Option<&'static str>,
}

const IMMUTABLE: &str = "IMMUTABLE: no later step may change it, S8 included.";
const S8_ONLY: &str = "MAY CHANGE ONLY IN S8 (the unified aloft vapour model, plan C7): \
                       byte-identical through every other step.";

/// Every pinned group, in file order.
fn groups() -> Vec<Group> {
    let group = |name, header, probe, compare| Group {
        name,
        header,
        probe,
        compare,
        feature: None,
    };
    #[allow(unused_mut)]
    let mut groups = vec![
        group(
            "TS_FLAT_DRY",
            "`TrajectorySolver::solve`, flat fire to 1000 m with sampling: Euler, RK4 and RK45 \
             at ICAO sea level (presence-blind) and at 1500 m, 25 C, 830 hPa (presence-aware), \
             without zones and with two dry zones. RH 0 throughout.",
            probe_ts_flat_dry,
            Compare::Bits,
        ),
        group(
            "TS_FLAT_HUMID",
            "`TrajectorySolver::solve`, flat fire, the same methods and stations: station RH 0 \
             under humid zones (RH 20 then 80), and station RH 50 and 100 with and without \
             those zones.",
            probe_ts_flat_humid,
            Compare::Bits,
        ),
        group(
            "TS_LOB_DRY",
            "`TrajectorySolver::solve`, the 30-degree lob to ground impact, Euler, RK4 and RK45 \
             at both stations, RH 0, no zones.",
            probe_ts_lob_dry,
            Compare::Bits,
        ),
        group(
            "TS_LOB_HUMID",
            "`TrajectorySolver::solve`, the 30-degree lob: RH 50 and 100 without zones, and RH 50 \
             under humid zones.",
            probe_ts_lob_humid,
            Compare::Bits,
        ),
        group(
            "TS_EFFECTS_DRY",
            "`TrajectorySolver::solve` (RK4, flat) with Litz spin drift plus aerodynamic jump, \
             and with Magnus, at both stations, RH 0.",
            probe_ts_effects_dry,
            Compare::Bits,
        ),
        group(
            "TS_EFFECTS_HUMID",
            "As `TS_EFFECTS_DRY`, at RH 50 and 100.",
            probe_ts_effects_humid,
            Compare::Bits,
        ),
        group(
            "ZERO_DRY",
            "`calculate_zero_angle_with_station_presence` to 300 m (every zero trial flies \
             `calculate_acceleration`), Euler, RK4 and RK45 at both stations, RH 0. The \
             bisection stops at 1e-7 rad, so a zero pin moves only when a trial height changes \
             sign; the flight groups are the sensitive ones.",
            probe_zero_dry,
            Compare::Bits,
        ),
        group(
            "ZERO_HUMID",
            "As `ZERO_DRY`, at RH 50 and 100.",
            probe_zero_humid,
            Compare::Bits,
        ),
        group(
            "MONTE_CARLO_DEFAULT_HUMIDITY",
            "Seeded `TrajectorySolver` Monte Carlo at the `BallisticInputs` default humidity \
             (50 %): the CLI's seeded entry at ICAO sea level, and the presence-aware entry the \
             FFI Monte Carlo runs on, at 1500 m.",
            probe_monte_carlo_default_humidity,
            Compare::Bits,
        ),
        group(
            "STATION_AIR",
            "Station density and sound speed: `calculate_atmosphere` over altitude, T/P override \
             and RH, plus `calculate_air_density_cimp` and `moist_speed_of_sound`. Station air \
             is not part of the aloft model, so this holds through S8 at every RH.",
            probe_station_air,
            Compare::Bits,
        ),
        group(
            "LOCAL_AIR_DRY",
            "`get_local_atmosphere_humid` from 200 m below each station base to 25 km above it, \
             RH 0.",
            probe_local_air_dry,
            Compare::Bits,
        ),
        group(
            "LOCAL_AIR_HUMID",
            "As `LOCAL_AIR_DRY`, at RH 37.5, 50 and 100. The legacy wrapper shares `air_at`'s \
             composition, so the aloft model reaches it.",
            probe_local_air_humid,
            Compare::Bits,
        ),
        group(
            "DENSITY_ALTITUDE",
            "NWS pressure altitude, the FAA-rule and density-matched density altitudes, the \
             density-altitude entry inverse, and the `atmosphere.density_altitude` service.",
            probe_density_altitude,
            Compare::Bits,
        ),
        group(
            "STABILITY_SPIN_JUMP",
            "Muzzle Miller Sg (the published dry T/P correction, kept deliberately: plan C8 was \
             rejected), Litz spin drift, the fast path's aerodynamic-jump offset, and \
             `TrajectorySolver`'s aerodynamic-jump components, which depend on station T/P only \
             and so are pinned at every RH.",
            probe_stability_spin_jump,
            Compare::Bits,
        ),
        group(
            "KERNEL_DIRECT",
            "Kernel direct (uniform) mode, density 1.1 kg/m3 and sound speed 335 m/s: \
             `compute_derivatives`, `try_integrate_trajectory` RK4 and RK45, \
             `try_solve_trajectory_rust`, `fast_integrate` and `fast_integrate_with_segments`. \
             S3 migrates the tuple to the uniform `KernelAtmosphere`; the values do not move.",
            probe_kernel_direct,
            Compare::Bits,
        ),
        group(
            "KERNEL_STANDARD_RH0",
            "Kernel standard mode at RH 0, station and zones dry, the same families plus the lean \
             Monte Carlo, at both stations. Compared to a relative 1e-15, not bits: the plan \
             allows last-ulp movement here when S3 moves the kernels onto `air_at`.",
            probe_kernel_standard_rh0,
            Compare::Relative,
        ),
        group(
            "FRONT_ENDS_DRY",
            "One in-process call through each JSON front end, RH 0: solve-json decoded and run \
             by `solve_v1` (Euler at sea level, RK45 at 1500 m), and `card_service::come_ups_v1` \
             at both stations.",
            probe_front_ends_dry,
            Compare::Bits,
        ),
        group(
            "FRONT_ENDS_HUMID",
            "As `FRONT_ENDS_DRY` with humidity: solve-json at RH 0.5 (RK4, sea level) and 0.9 \
             (RK45, 1500 m), and the come-ups card at 50 and 100 %.",
            probe_front_ends_humid,
            Compare::Bits,
        ),
    ];
    #[cfg(feature = "bridge")]
    groups.extend([
        Group {
            feature: Some("bridge"),
            ..group(
                "BRIDGE_SOLVE_DRY",
                "The bridge `solve` command (`bridge_call`) on the `FRONT_ENDS_DRY` requests.",
                probe_bridge_solve_dry,
                Compare::Bits,
            )
        },
        Group {
            feature: Some("bridge"),
            ..group(
                "BRIDGE_SOLVE_HUMID",
                "The bridge `solve` command on the `FRONT_ENDS_HUMID` requests.",
                probe_bridge_solve_humid,
                Compare::Bits,
            )
        },
    ]);
    #[cfg(feature = "ffi")]
    groups.extend([
        Group {
            feature: Some("ffi"),
            ..group(
                "FFI_TRAJECTORY_DRY",
                "The FFI single solve `ballistics_calculate_trajectory` (RK4, sampling and spin \
                 drift on) at both stations, RH 0.",
                probe_ffi_trajectory_dry,
                Compare::Bits,
            )
        },
        Group {
            feature: Some("ffi"),
            ..group(
                "FFI_TRAJECTORY_HUMID",
                "As `FFI_TRAJECTORY_DRY`, at RH 50 and 100.",
                probe_ffi_trajectory_humid,
                Compare::Bits,
            )
        },
    ]);
    groups
}

/// Whether a group's constant may change in S8: every group whose cases involve an RH above
/// zero, at the station or in a zone, except the station-air and stability groups.
fn s8_may_change(name: &str) -> bool {
    name.ends_with("_HUMID") || name == "MONTE_CARLO_DEFAULT_HUMIDITY"
}

/// Wrap a header at 96 columns as `///` lines.
fn print_doc(text: &str) {
    let mut line = String::new();
    for word in text.split_whitespace() {
        if !line.is_empty() && line.len() + 1 + word.len() > 96 {
            println!("/// {line}");
            line.clear();
        }
        if !line.is_empty() {
            line.push(' ');
        }
        line.push_str(word);
    }
    if !line.is_empty() {
        println!("/// {line}");
    }
}

/// Regenerate the constants below from the current engine. Never run to "fix" a failing pin:
/// a moved pin is a physics change to explain, not a constant to refresh.
#[test]
#[ignore = "generator: prints the pinned constants from the current engine"]
fn print_current_pin_values() {
    for group in groups() {
        let probe = (group.probe)();
        print_doc(group.header);
        println!("///");
        print_doc(if s8_may_change(group.name) {
            S8_ONLY
        } else {
            IMMUTABLE
        });
        if let Compare::Relative = group.compare {
            println!("///");
            print_doc(&format!(
                "Floats compared to a relative {KERNEL_STANDARD_RH0_REL_TOL:e}; counts exactly."
            ));
        }
        if let Some(feature) = group.feature {
            println!("#[cfg(feature = {feature:?})]");
        }
        // rustfmt would put every 18-character literal on its own line.
        println!("#[rustfmt::skip]");
        print_probe(group.name, &probe);
    }
}

macro_rules! pin_test {
    ($test:ident, $probe:ident, $pinned:ident) => {
        #[test]
        fn $test() {
            assert_bits(stringify!($pinned), &$probe(), $pinned);
        }
    };
    ($test:ident, $probe:ident, $pinned:ident, relative) => {
        #[test]
        fn $test() {
            assert_relative(
                stringify!($pinned),
                &$probe(),
                $pinned,
                KERNEL_STANDARD_RH0_REL_TOL,
            );
        }
    };
}

pin_test!(
    trajectory_solver_flat_dry_is_byte_identical,
    probe_ts_flat_dry,
    TS_FLAT_DRY
);
pin_test!(
    trajectory_solver_flat_humid_is_byte_identical,
    probe_ts_flat_humid,
    TS_FLAT_HUMID
);
pin_test!(
    trajectory_solver_lob_dry_is_byte_identical,
    probe_ts_lob_dry,
    TS_LOB_DRY
);
pin_test!(
    trajectory_solver_lob_humid_is_byte_identical,
    probe_ts_lob_humid,
    TS_LOB_HUMID
);
pin_test!(
    trajectory_solver_effects_dry_are_byte_identical,
    probe_ts_effects_dry,
    TS_EFFECTS_DRY
);
pin_test!(
    trajectory_solver_effects_humid_are_byte_identical,
    probe_ts_effects_humid,
    TS_EFFECTS_HUMID
);
pin_test!(zero_angles_dry_are_byte_identical, probe_zero_dry, ZERO_DRY);
pin_test!(
    zero_angles_humid_are_byte_identical,
    probe_zero_humid,
    ZERO_HUMID
);
pin_test!(
    seeded_monte_carlo_at_default_humidity_is_byte_identical,
    probe_monte_carlo_default_humidity,
    MONTE_CARLO_DEFAULT_HUMIDITY
);
pin_test!(
    station_air_is_byte_identical,
    probe_station_air,
    STATION_AIR
);
pin_test!(
    local_air_dry_is_byte_identical,
    probe_local_air_dry,
    LOCAL_AIR_DRY
);
pin_test!(
    local_air_humid_is_byte_identical,
    probe_local_air_humid,
    LOCAL_AIR_HUMID
);
pin_test!(
    density_altitude_is_byte_identical,
    probe_density_altitude,
    DENSITY_ALTITUDE
);
pin_test!(
    stability_spin_drift_and_jump_are_byte_identical,
    probe_stability_spin_jump,
    STABILITY_SPIN_JUMP
);
pin_test!(
    kernel_direct_mode_is_byte_identical,
    probe_kernel_direct,
    KERNEL_DIRECT
);
pin_test!(
    kernel_standard_mode_at_rh0_holds_to_1e15,
    probe_kernel_standard_rh0,
    KERNEL_STANDARD_RH0,
    relative
);
pin_test!(
    front_ends_dry_are_byte_identical,
    probe_front_ends_dry,
    FRONT_ENDS_DRY
);
pin_test!(
    front_ends_humid_are_byte_identical,
    probe_front_ends_humid,
    FRONT_ENDS_HUMID
);
#[cfg(feature = "bridge")]
pin_test!(
    bridge_solve_dry_is_byte_identical,
    probe_bridge_solve_dry,
    BRIDGE_SOLVE_DRY
);
#[cfg(feature = "bridge")]
pin_test!(
    bridge_solve_humid_is_byte_identical,
    probe_bridge_solve_humid,
    BRIDGE_SOLVE_HUMID
);
#[cfg(feature = "ffi")]
pin_test!(
    ffi_trajectory_dry_is_byte_identical,
    probe_ffi_trajectory_dry,
    FFI_TRAJECTORY_DRY
);
#[cfg(feature = "ffi")]
pin_test!(
    ffi_trajectory_humid_is_byte_identical,
    probe_ffi_trajectory_humid,
    FFI_TRAJECTORY_HUMID
);

// ----- PINNED CONSTANTS (generated; see the module docs) -----------------------------------------

/// `TrajectorySolver::solve`, flat fire to 1000 m with sampling: Euler, RK4 and RK45 at ICAO sea
/// level (presence-blind) and at 1500 m, 25 C, 830 hPa (presence-aware), without zones and with two
/// dry zones. RH 0 throughout.
///
/// IMMUTABLE: no later step may change it, S8 included.
#[rustfmt::skip]
const TS_FLAT_DRY: Pinned = &[
    (
        "flat/euler/icao-sea-level/rh0/no-zones",
        &[
            0x408f400000000000, 0x3fc06b98d7967e21, 0x3fff76297d72230e, 0x407430ebb4f7475c,
            0x40827e0cc14be09c, 0x407544b4bbbd89be, 0x0000000000000001, 0x4088244f84fd18b2,
            0x0000000000000001, 0x408d392d9d8a5a35, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007b0,
            0x3f50624dd2f1a9fc, 0x3fea51b773315f67, 0x3f5aca4c8688729e, 0xbec554a6c13e5dde,
            0x4089b3d07aac0133, 0x0000000000000001, 0x3fd30dc66fa4c7de, 0x3fdf7ced916872b7,
            0x4075e449a7d29d20, 0xbfd87bfac0055c97, 0xbfd0ff9987d11b44, 0x408346d61179c286,
            0x0000000000000001, 0x3fd594d366d1d6d4, 0x3fef7ced916872b7, 0x40835059d9a2d975,
            0xc0062a3a3019a489, 0xbfedd40df4fd7b36, 0x407e03c5796d4013, 0x0000000000000001,
            0x3fd8c1aadddbadbd, 0x3ff79db22d0e551b, 0x4089e8c5d0703d80, 0xc01b192ac462d84d,
            0xbffdfd61d7d836f2, 0x4077f80f564e8e3a, 0x0000000000000001, 0x3fdbae108d8ab609,
            0x3fff76297d72230e, 0x408f400000000000, 0xc028557d9f64ff43, 0xc0080a7a0dc4fdff,
            0x407430ebb4f7475c, 0x0000000000000001, 0x3fcc4caabfbfb6d5, 0x1aac9a4b637ca177,
            0x0000000000000015, 0x83eba179114eb2cd,
        ],
    ),
    (
        "flat/euler/icao-sea-level/rh0/zones-rh0-rh0",
        &[
            0x408f400000000000, 0x3fc06dd61c346a0f, 0x3fff73969366a763, 0x40741ae5c9094414,
            0x408255cb91e733fa, 0x407544b4bbbd89be, 0x0000000000000001, 0x408823ffd5b28fb8,
            0x0000000000000001, 0x408d2709e1c48da7, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007af,
            0x3f50624dd2f1a9fc, 0x3fea51c07cb26ccb, 0x3f5aca55c7907a59, 0xbec527ac9f17c6af,
            0x4089b3d94df52e2e, 0x0000000000000001, 0x3fd30dc3579a8801, 0x3fdf6c8b4395810d,
            0x4075e0ad9cd1ff82, 0xbfd84addac105e01, 0xbfd0d00734b70298, 0x4083547e30e13c8c,
            0x0000000000000001, 0x3fd58dc1cbaaa6ad, 0x3fef74bc6a7ef9e2, 0x408353287d2aa1ba,
            0xc0061ad22c7d11aa, 0xbfeda550bd1347ce, 0x407e0ab2e61e9a1f, 0x0000000000000001,
            0x3fd8be6c150a519b, 0x3ff79999999998b1, 0x4089eb46b6d35ed9, 0xc01b06f2271e7b02,
            0xbffde6c4af1f6020, 0x4077f01691067a49, 0x0000000000000001, 0x3fdbb12a21940320,
            0x3fff73969366a763, 0x408f400000000000, 0xc0284534ba13f53e, 0xc008055439ae0691,
            0x40741ae5c9094414, 0x0000000000000001, 0x3fcafff5671916ca, 0x537717c2d126dcba,
            0x0000000000000015, 0xb9ba68e34787a7ec,
        ],
    ),
    (
        "flat/euler/1500m-25C-830hPa/rh0/no-zones",
        &[
            0x408f400000000000, 0x3fc0a350080d9c9b, 0x3ffbfe14213b4ab4, 0x4078d06254767b67,
            0x408bee024596814d, 0x4075a2608de6e31f, 0x0000000000000001, 0x408da53b79493a75,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d7,
            0x3f50624dd2f1a9fc, 0x3fea529531f092df, 0x3f5acb2f97c26211, 0xbec1050afce375ad,
            0x4089b4a9070eea46, 0x0000000000000001, 0x3fd33490ab7111c8, 0x3fdbf7ced9168731,
            0x40743ca53883280c, 0xbfcd4f311cf42519, 0xbfc64f1c48a9acd9, 0x4084dde4daf0eb92,
            0x0000000000000001, 0x3fd50009c67099ca, 0x3fec000000000006, 0x4082702152a2b5f5,
            0xc0010787607d22fa, 0xbfe43c45baa111dd, 0x40814c898f29a28e, 0x0000000000000001,
            0x3fd708d018685549, 0x3ff5020c49ba5d9d, 0x40896314ba2052c5, 0xc015c44bd075a5a7,
            0xbff4db7a1d54b4d8, 0x407d1bc10b6404f2, 0x0000000000000001, 0x3fd97392635a8384,
            0x3ffbfe14213b4ab4, 0x408f400000000000, 0xc023ee2d01fa8ff2, 0xc0011a4f55574c4d,
            0x4078d06254767b67, 0x0000000000000001, 0x3fdb80aa677c4282, 0x5b883586de32a6c4,
            0x0000000000000015, 0x254c68a6f261b0fb,
        ],
    ),
    (
        "flat/euler/1500m-25C-830hPa/rh0/zones-rh0-rh0",
        &[
            0x408f400000000000, 0x3fc0a546616af160, 0x3ffbfa7435b1aa25, 0x4078c2d9e9c2bef5,
            0x408bcf93bd674694, 0x4075a2608de6e31f, 0x0000000000000001, 0x408d950d343438f4,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d6,
            0x3f50624dd2f1a9fc, 0x3fea529cd8eca017, 0x3f5acb376dc31ee2, 0xbec0def545676cae,
            0x4089b4b08022086d, 0x0000000000000001, 0x3fd3348e1768371d, 0x3fdbf7ced9168731,
            0x4074412d92324962, 0xbfcd4cb4f12d991e, 0xbfc621fe55724bc5, 0x4084e7376f352436,
            0x0000000000000001, 0x3fd4fc45d9536b7b, 0x3fec000000000006, 0x408275bac50decb3,
            0xc00106d6a1dc6845, 0xbfe4206792ad0542, 0x40814f9b259217b4, 0x0000000000000001,
            0x3fd7068675e86fdb, 0x3ff4fdf3b645a133, 0x408964b709cf4e82, 0xc015b636c01b443a,
            0xbff4c706dbe4863b, 0x407d176431086026, 0x0000000000000001, 0x3fd975af10f6534b,
            0x3ffbfa7435b1aa25, 0x408f400000000000, 0xc023e1b1d6927107, 0xc001130f7e440b56,
            0x4078c2d9e9c2bef5, 0x0000000000000001, 0x3fdb8680b2b003a7, 0x4c22fcafa07e5293,
            0x0000000000000015, 0xcd265b6e48e03ce5,
        ],
    ),
    (
        "flat/rk4/icao-sea-level/rh0/no-zones",
        &[
            0x408f400000000000, 0x3fc087dde518e09a, 0x3fff710706a4eb9d, 0x40743327745a4e9b,
            0x4082822445abdee0, 0x407544b4bbbd89be, 0x0000000000000001, 0x40882a05d1a01309,
            0x0000000000000001, 0x408d3ef253a771b8, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007af,
            0x3f50624dd2f1a9fc, 0x3fea53dc3977ed14, 0x3f5ae1107de50453, 0xbeb552e290a8a792,
            0x4089b3d0ffc19145, 0x0000000000000001, 0x3fd30dc640fcfb5d, 0x3fdf6c8b4395810d,
            0x4075dc6a8374e748, 0xbfd820173f7facbf, 0xbfd0e53d85d33a99, 0x40834a1a82ae0e9f,
            0x0000000000000001, 0x3fd59320bd9fde15, 0x3fef74bc6a7ef9e2, 0x40834e660fd52e27,
            0xc006122af243ee6f, 0xbfedbd00c1d76347, 0x407e08fcc8c184fc, 0x0000000000000001,
            0x3fd8bf3940c7b445, 0x3ff79999999998b1, 0x4089e85c08b8f26d, 0xc01b079c83b0c3d8,
            0xbffdee06c4dc4b28, 0x4077fc3971e5169c, 0x0000000000000001, 0x3fdbac6eb73ad5bd,
            0x3fff710706a4eb9d, 0x408f400000000000, 0xc02848423d96c88a, 0xc0080035202a8f1b,
            0x40743327745a4e9b, 0x0000000000000001, 0x3fcc72a456a11b12, 0x86c65b32859a15e5,
            0x0000000000000015, 0xcf18b8e5b4398abf,
        ],
    ),
    (
        "flat/rk4/icao-sea-level/rh0/zones-rh0-rh0",
        &[
            0x408f400000000000, 0x3fc08a1ac20a15f3, 0x3fff6e7a92a4822b, 0x40741d23a7de38ce,
            0x408259e27f9963e8, 0x407544b4bbbd89be, 0x0000000000000001, 0x4088266b09abf43f,
            0x0000000000000001, 0x408d2cc4759184c5, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ae,
            0x3f50624dd2f1a9fc, 0x3fea53e0bd8822ac, 0x3f5ae1151b5a9df2, 0xbeb525ef4989d81b,
            0x4089b3d9d105bf74, 0x0000000000000001, 0x3fd30dc329a7f956, 0x3fdf6c8b4395810d,
            0x4075e2a90d8abbb7, 0xbfd81fc278b2016d, 0xbfd0c629b6d44428, 0x4083551641289941,
            0x0000000000000001, 0x3fd58d738199795a, 0x3fef74bc6a7ef9e2, 0x4083550761865394,
            0xc0061044cf72f282, 0xbfed9c01545b5010, 0x407e0c26b8010763, 0x0000000000000001,
            0x3fd8bdbdfe4b94eb, 0x3ff795810624dc47, 0x4089ead6ccd014c3, 0xc01af56d0ceab66d,
            0xbffdd778dfc326a4, 0x4077f43d0446e15c, 0x0000000000000001, 0x3fdbaf8e0396f785,
            0x3fff6e7a92a4822b, 0x408f400000000000, 0xc028380f8ad7fce3, 0xc007fb1c3829bc0f,
            0x40741d23a7de38ce, 0x0000000000000001, 0x3fcb1e98c9b6385a, 0xc7b52f1d14f49005,
            0x0000000000000015, 0xc5c7977ce3cd4073,
        ],
    ),
    (
        "flat/rk4/1500m-25C-830hPa/rh0/no-zones",
        &[
            0x408f400000000000, 0x3fc0bf5b66092374, 0x3ffbfacd69a594c0, 0x4078d362f3ddd347,
            0x408bf4c4f2821158, 0x4075a2608de6e31f, 0x0000000000000001, 0x408da79f31e0457c,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d6,
            0x3f50624dd2f1a9fc, 0x3fea544b088d78cb, 0x3f5ae181bc23bf88, 0xbeb103e8f0646846,
            0x4089b4a95c6c7265, 0x0000000000000001, 0x3fd334908dfd512e, 0x3fdbf7ced9168731,
            0x40743e151d24a139, 0xbfcd037fbece5a1f, 0xbfc640ce550a0a41, 0x4084de4a01e3d32a,
            0x0000000000000001, 0x3fd4ffe0d228497d, 0x3fec000000000006, 0x40827184e7b19440,
            0xc000fe2fe6e4fdd0, 0xbfe4355c05f43ee7, 0x40814d1763cf17f6, 0x0000000000000001,
            0x3fd7086637e24f4c, 0x3ff4fdf3b645a133, 0x4089615a15eb1f48, 0xc015b42b79832c11,
            0xbff4cf654e90ed8c, 0x407d1fcdd9715983, 0x0000000000000001, 0x3fd9719be44b0fbe,
            0x3ffbfacd69a594c0, 0x408f400000000000, 0xc023e4e3a6711b64, 0xc00113c1e62be087,
            0x4078d362f3ddd347, 0x0000000000000001, 0x3fdb7f5c3c8ad9d2, 0xc11789495d22f71a,
            0x0000000000000015, 0x719d8bddea2f2b86,
        ],
    ),
    (
        "flat/rk4/1500m-25C-830hPa/rh0/zones-rh0-rh0",
        &[
            0x408f400000000000, 0x3fc0c14cf31b67ed, 0x3ffbf735944bebcb, 0x4078c5d1498c3947,
            0x408bd63df1baef00, 0x4075a2608de6e31f, 0x0000000000000001, 0x408d9aba5c7927c7,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d5,
            0x3f50624dd2f1a9fc, 0x3fea544edb8ca76c, 0x3f5ae185a4a48e5e, 0xbeb0ddd8277d7b22,
            0x4089b4b0d40bc377, 0x0000000000000001, 0x3fd3348dfa74c36c, 0x3fdbf7ced9168731,
            0x4074429a54bc1af0, 0xbfcd010ecf60374d, 0xbfc613cf91ba6cec, 0x4084e79b2b8390ab,
            0x0000000000000001, 0x3fd4fc1daab9587d, 0x3febf7ced9168731, 0x408272abb5f2b464,
            0xc000f12b220a9f10, 0xbfe40edc309c0a2e, 0x408151f6558f6a8c, 0x0000000000000001,
            0x3fd704c60ca4ea4a, 0x3ff4f9db22d0e4c9, 0x408962f6c4f722f9, 0xc015a61d8083d810,
            0xbff4bb000e18dc0c, 0x407d1b6c6e70604d, 0x0000000000000001, 0x3fd973bb5f3c75dd,
            0x3ffbf735944bebcb, 0x408f400000000000, 0xc023d878ac3fed62, 0xc0010c923b788e74,
            0x4078c5d1498c3947, 0x0000000000000001, 0x3fdb853ad6de0980, 0x5a06dd3153a14f8d,
            0x0000000000000015, 0xfda2bac7b04dbcc0,
        ],
    ),
    (
        "flat/rk45/icao-sea-level/rh0/no-zones",
        &[
            0x408f400000000000, 0x3fc087134a1b51e7, 0x3fff71092697ef90, 0x4074332de0869631,
            0x408282300aabd72b, 0x407544b4bbbd89be, 0x0000000000000001, 0x408843ae3cfe07e2,
            0x0000000000000001, 0x408d415170d0e628, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000c9,
            0x3f50624dd2f1a9fc, 0x3fea53dc3977ec61, 0x3f5ae1107de503f6, 0xbeb552e290af8a29,
            0x4089b3d0ffc19145, 0x0000000000000001, 0x3fd30dc640fcfb5d, 0x3fde33269df97ab0,
            0x40751e868614f03e, 0xbfd49717a1edf6d4, 0xbfcf617a32478f3e, 0x40837d83f7020916,
            0x0000000000000001, 0x3fd578ecf9ed8976, 0x3fef19934efcbd5c, 0x408323819e2e2383,
            0xc0057ce1ee97bcfe, 0xbfed25d528b10142, 0x407e32655691a28f, 0x0000000000000001,
            0x3fd8abe1af8d4c0a, 0x3ff78cc9a77e5eb0, 0x4089dec0bbed76db, 0xc01ae8ed3c8031c5,
            0xbffdd2aedc383320, 0x4078047f76e39c2e, 0x0000000000000001, 0x3fdba92a038ff66d,
            0x3fff71092697ef90, 0x408f400000000000, 0xc028485039f65c70, 0xc008003960109324,
            0x4074332de0869631, 0x0000000000000001, 0x3fcc7312c34c6c75, 0xb5cabf62ef0684f0,
            0x0000000000000015, 0x286921b765f8c763,
        ],
    ),
    (
        "flat/rk45/icao-sea-level/rh0/zones-rh0-rh0",
        &[
            0x408f400000000000, 0x3fc08959c92df0bc, 0x3fff6e826ec0fedd, 0x40741d2697e1af0a,
            0x408259e7dbd1cb8a, 0x407544b4bbbd89be, 0x0000000000000001, 0x408829415f99b47e,
            0x0000000000000001, 0x408d44d8494f84cc, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000c9,
            0x3f50624dd2f1a9fc, 0x3fea53e0bd88220a, 0x3f5ae1151b5a9db5, 0xbeb525ef49902ed6,
            0x4089b3d9d105bf74, 0x0000000000000001, 0x3fd30dc329a7f956, 0x3fde33269df97ab0,
            0x4075245abad6ac85, 0xbfd49655bff187ce, 0xbfcf2775013c935b, 0x4083883eba961642,
            0x0000000000000001, 0x3fd573905e6efa45, 0x3fef19934efcbd5c, 0x40832a1c9ee7c7ea,
            0xc0057b2680a57da2, 0xbfed04f52aba28ab, 0x407e35f492287609, 0x0000000000000001,
            0x3fd8aa39451065ca, 0x3ff78cc9a77e5eb0, 0x4089e44d5b473c29, 0xc01ae098bd29894d,
            0xbffdc4dfb2085c02, 0x4077f9e51c96a81f, 0x0000000000000001, 0x3fdbad58b5dd64b6,
            0x3fff6e826ec0fedd, 0x408f400000000000, 0xc02838241f407985, 0xc007fb2bf062b1ba,
            0x40741d2697e1af0a, 0x0000000000000001, 0x3fcb1ec162fe8b42, 0x02f13fa5645c49d8,
            0x0000000000000015, 0xa997a95b3938ef17,
        ],
    ),
    (
        "flat/rk45/1500m-25C-830hPa/rh0/no-zones",
        &[
            0x408f400000000000, 0x3fc0bf489fb95561, 0x3ffbfad1572b59b3, 0x4078d3644023e104,
            0x408bf4c7deda0f3f, 0x4075a2608de6e31f, 0x0000000000000001, 0x408dc50b020c75d1,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000b3,
            0x3f50624dd2f1a9fc, 0x3fea544b088d788b, 0x3f5ae181bc23bf7c, 0xbeb103e8f066e978,
            0x4089b4a95c6c7265, 0x0000000000000001, 0x3fd334908dfd512e, 0x3fda5c1c6088d6d8,
            0x407330132e8b787a, 0xbfc5322c4fd750de, 0xbfc3e2d7b9dec68b, 0x40851b4221ab1314,
            0x0000000000000001, 0x3fd4e779e923a841, 0x3feb947496aad1d6, 0x408237389439abad,
            0xc0005d6aa27e8588, 0xbfe3a9539bf56dab, 0x408164fa0880cf98, 0x0000000000000001,
            0x3fd6f6c8250c3e72, 0x3ff4fd6d7e889c20, 0x408960dff34fb17c, 0xc015b2fd0d1c77c1,
            0xbff4ce7c5d5f1815, 0x407d202ac73642fd, 0x0000000000000001, 0x3fd9716ed352d6b8,
            0x3ffbfad1572b59b3, 0x408f400000000000, 0xc023e4f15dedd7ee, 0xc00113c9c1376775,
            0x4078d3644023e104, 0x0000000000000001, 0x3fdb7f5babe07137, 0xbf5e375138f4c81a,
            0x0000000000000015, 0xf391851593e3055a,
        ],
    ),
    (
        "flat/rk45/1500m-25C-830hPa/rh0/zones-rh0-rh0",
        &[
            0x408f400000000000, 0x3fc0c13e164f5331, 0x3ffbf73e984bb93c, 0x4078c5cc17b65839,
            0x408bd63245165630, 0x4075a2608de6e31f, 0x0000000000000001, 0x408da790664cf3a8,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000b3,
            0x3f50624dd2f1a9fc, 0x3fea544edb8ca732, 0x3f5ae185a4a48e60, 0xbeb0ddd8277fc85a,
            0x4089b4b0d40bc376, 0x0000000000000001, 0x3fd3348dfa74c36c, 0x3fda5c1c6088d6d8,
            0x40733422f67e43b2, 0xbfc52f26bc95bdea, 0xbfc3ba69e779d110, 0x40852431c3f892ea,
            0x0000000000000001, 0x3fd4e3efe3b2b40d, 0x3feb947496aad1d6, 0x40823cc25db2ab7e,
            0xc0005cd82f914fa4, 0xbfe38dc364397435, 0x40816842c80b839f, 0x0000000000000001,
            0x3fd6f463cb0dd270, 0x3ff4fd6d7e889c20, 0x40896633a6b0de3d, 0xc015ae23c8a3852b,
            0xbff4c13ad8fa83f0, 0x407d18e8960df8ce, 0x0000000000000001, 0x3fd974f31c997c8e,
            0x3ffbf73e984bb93c, 0x408f400000000000, 0xc023d88d903fb8d0, 0xc0010ca443782689,
            0x4078c5cc17b65839, 0x0000000000000001, 0x3fdb853d12580baa, 0x878500d03a31a800,
            0x0000000000000015, 0x4d2911556051cc32,
        ],
    ),
];

/// `TrajectorySolver::solve`, flat fire, the same methods and stations: station RH 0 under humid
/// zones (RH 20 then 80), and station RH 50 and 100 with and without those zones.
///
/// MAY CHANGE ONLY IN S8 (the unified aloft vapour model, plan C7): byte-identical through every
/// other step.
#[rustfmt::skip]
const TS_FLAT_HUMID: Pinned = &[
    (
        "flat/euler/icao-sea-level/rh0/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc06e26691f431f, 0x3fff6998eb8528b6, 0x407427080b40a714,
            0x40826bf3da3247c9, 0x407544b4bbbd89be, 0x0000000000000001, 0x408830fa6eed113a,
            0x0000000000000001, 0x408d39ebe3af1746, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ad,
            0x3f50624dd2f1a9fc, 0x3fea51c1c06bf31a, 0x3f5aca57130f1490, 0xbec5216177e29fc7,
            0x4089b3da8a188477, 0x0000000000000001, 0x3fd30dc2e8c659d5, 0x3fdf6c8b4395810d,
            0x4075e18ed3196eba, 0xbfd84ad189a7329c, 0xbfd0cba657c49df1, 0x4083560a1f0fd0d9,
            0x0000000000000001, 0x3fd58cf605fbb7a7, 0x3fef6c8b4395810d, 0x40835108b65b8313,
            0xc0060ec45678c4a4, 0xbfed8f1e7717a25b, 0x407e15d3f4596c9f, 0x0000000000000001,
            0x3fd8b936f89b3a21, 0x3ff7916872b01fdd, 0x4089e8fe441c9c83, 0xc01af68c0b1fbb81,
            0xbffdcbae73ef59f9, 0x4077feff38505d64, 0x0000000000000001, 0x3fdbab574cd477d7,
            0x3fff6998eb8528b6, 0x408f400000000000, 0xc028394b9db1464a, 0xc007f158e9eb093d,
            0x407427080b40a714, 0x0000000000000001, 0x3fcbae3e44c6c694, 0xd552ce8d8dec436d,
            0x0000000000000015, 0x8f3d3f4c29539c3e,
        ],
    ),
    (
        "flat/euler/icao-sea-level/rh50/no-zones",
        &[
            0x408f400000000000, 0x3fc06c4a2bc26ca5, 0x3fff69ffc4732261, 0x40743c236b2fb2c3,
            0x4082929eae800c6e, 0x40754bf73b88949e, 0x0000000000000001, 0x40882c3e8af6a63f,
            0x0000000000000001, 0x408d42b559a85302, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ad,
            0x3f50624dd2f1a9fc, 0x3fea51ba3e4e3ed6, 0x3f5aca4f62cf3933, 0xbec546bfb16f49c3,
            0x4089b3d335068e3e, 0x0000000000000001, 0x3fd310d7d819611a, 0x3fdf6c8b4395810d,
            0x4075dc5bc676cc3d, 0xbfd84b1ea44d6d69, 0xbfd0e586dfb8a9f9, 0x40834ce90c4fe198,
            0x0000000000000001, 0x3fd59516b21824b3, 0x3fef6c8b4395810d, 0x40834ba7d6a834c6,
            0xc00610a3c6783bf8, 0xbfeda9e30fa9566d, 0x407e14632b3ccf12, 0x0000000000000001,
            0x3fd8beafd5671de5, 0x3ff7916872b01fdd, 0x4089e4eb50d6ce66, 0xc01afe49b58e9d34,
            0xbffdd5d1cb7ba7ba, 0x407809d418e06082, 0x0000000000000001, 0x3fdbaa4c9d2ae0bc,
            0x3fff69ffc4732261, 0x408f400000000000, 0xc028450d5071e1b5, 0xc007f2269bc6fc86,
            0x40743c236b2fb2c3, 0x0000000000000001, 0x3fcc96e62a38b8a5, 0x90e01a862cd9485a,
            0x0000000000000015, 0x721aee84cf0ca201,
        ],
    ),
    (
        "flat/euler/icao-sea-level/rh50/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc06e26691f431f, 0x3fff6998eb8528b6, 0x407427080b40a714,
            0x40826bf3da3247c9, 0x40754bf73b88949e, 0x0000000000000001, 0x4088272da017dd7a,
            0x0000000000000001, 0x408d2f08da7732a7, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ad,
            0x3f50624dd2f1a9fc, 0x3fea51c1c06bf31a, 0x3f5aca57130f1490, 0xbec5216177e29fc7,
            0x4089b3da8a188477, 0x0000000000000001, 0x3fd310d546ce290d, 0x3fdf6c8b4395810d,
            0x4075e18ed3196eba, 0xbfd84ad189a7329c, 0xbfd0cba657c49df1, 0x4083560a1f0fd0d9,
            0x0000000000000001, 0x3fd5905c17234798, 0x3fef6c8b4395810d, 0x40835108b65b8313,
            0xc0060ec45678c4a4, 0xbfed8f1e7717a25b, 0x407e15d3f4596c9f, 0x0000000000000001,
            0x3fd8be0361eba2cf, 0x3ff7916872b01fdd, 0x4089e8fe441c9c83, 0xc01af68c0b1fbb81,
            0xbffdcbae73ef59f9, 0x4077feff38505d64, 0x0000000000000001, 0x3fdbae8cd2aac9ed,
            0x3fff6998eb8528b6, 0x408f400000000000, 0xc028394b9db1464a, 0xc007f158e9eb093d,
            0x407427080b40a714, 0x0000000000000001, 0x3fcb48ffb8563578, 0xd552ce8d8dec436d,
            0x0000000000000015, 0x8f3d3f4c29539c3e,
        ],
    ),
    (
        "flat/euler/icao-sea-level/rh100/no-zones",
        &[
            0x408f400000000000, 0x3fc06cfb6cc104de, 0x3fff5dde6339e657, 0x40744765c9655fdd,
            0x4082a74fa0cdebaf, 0x40755334ef1fc968, 0x0000000000000001, 0x408830eabda929a8,
            0x0000000000000001, 0x408d4c41173c3f66, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007aa,
            0x3f50624dd2f1a9fc, 0x3fea51bd09223bce, 0x3f5aca523ecb5d7b, 0xbec538da0c5e0036,
            0x4089b3d5ef19ee17, 0x0000000000000001, 0x3fd313e577ef400d, 0x3fdf5c28f5c28f63,
            0x4075d468fe23d0ac, 0xbfd81a610c05664a, 0xbfd0cb8c9d7e4ce4, 0x408352faab74ecda,
            0x0000000000000001, 0x3fd5955848e81f16, 0x3fef645a1cac0838, 0x40834acab6aa5c6e,
            0xc00604905c621627, 0xbfed8d6af9aae05d, 0x407e214de9c255b0, 0x0000000000000001,
            0x3fd8bd6e4ce8119a, 0x3ff789374bc6a709, 0x4089e41cb3b3d2b7, 0xc01aed3da0409837,
            0xbffdb70f566325a4, 0x407818fab4af23d2, 0x0000000000000001, 0x3fdba788a23bc34f,
            0x3fff5dde6339e657, 0x408f400000000000, 0xc02834b05280ec67, 0xc007d9e3d9548471,
            0x40744765c9655fdd, 0x0000000000000001, 0x3fcce5eac1d00cfc, 0xc3cbfeda5c53cae7,
            0x0000000000000015, 0x14a63a57601e3d7e,
        ],
    ),
    (
        "flat/euler/icao-sea-level/rh100/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc06e26691f431f, 0x3fff6998eb8528b6, 0x407427080b40a714,
            0x40826bf3da3247c9, 0x40755334ef1fc968, 0x0000000000000001, 0x40881d5d6cd3d545,
            0x0000000000000001, 0x408d26dbd32e5677, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ad,
            0x3f50624dd2f1a9fc, 0x3fea51c1c06bf31a, 0x3f5aca57130f1490, 0xbec5216177e29fc7,
            0x4089b3da8a188477, 0x0000000000000001, 0x3fd313e3db5dfc03, 0x3fdf6c8b4395810d,
            0x4075e18ed3196eba, 0xbfd84ad189a7329c, 0xbfd0cba657c49df1, 0x4083560a1f0fd0d9,
            0x0000000000000001, 0x3fd593c1e4a0095a, 0x3fef6c8b4395810d, 0x40835108b65b8313,
            0xc0060ec45678c4a4, 0xbfed8f1e7717a25b, 0x407e15d3f4596c9f, 0x0000000000000001,
            0x3fd8c2caff1777b8, 0x3ff7916872b01fdd, 0x4089e8fe441c9c83, 0xc01af68c0b1fbb81,
            0xbffdcbae73ef59f9, 0x4077feff38505d64, 0x0000000000000001, 0x3fdbb1b57685b2cd,
            0x3fff6998eb8528b6, 0x408f400000000000, 0xc028394b9db1464a, 0xc007f158e9eb093d,
            0x407427080b40a714, 0x0000000000000001, 0x3fcaeae80ee4cc53, 0xd552ce8d8dec436d,
            0x0000000000000015, 0x8f3d3f4c29539c3e,
        ],
    ),
    (
        "flat/euler/1500m-25C-830hPa/rh0/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc0a5dc4b34a33b, 0x3ffbebc74a7c295f, 0x4078e241424554fa,
            0x408c164b80128175, 0x4075a2608de6e31f, 0x0000000000000001, 0x408dc0315f9b2271,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d3,
            0x3f50624dd2f1a9fc, 0x3fea529f2159466f, 0x3f5acb39c436ada5, 0xbec0d398a2a9eea0,
            0x4089b4b2badc6d3b, 0x0000000000000001, 0x3fd3348d527fcc08, 0x3fdbe76c8b439587,
            0x407437d3ddce2dfa, 0xbfccf8af763078a8, 0xbfc5fbfe78f7c96a, 0x4084ec66fdef8fe1,
            0x0000000000000001, 0x3fd4fa2f982ae1c0, 0x3febef9db22d0e5c, 0x40826f9b18cd97fd,
            0xc000efee7e1db1b7, 0xbfe3fd588b8d1685, 0x40815a5ef5a3e066, 0x0000000000000001,
            0x3fd6fe8e226210d9, 0x3ff4f5c28f5c285f, 0x408963d6f7c01710, 0xc015a853282f0aa3,
            0xbff4a86fd78a6fdf, 0x407d3158348cb33e, 0x0000000000000001, 0x3fd96914bd1c14fc,
            0x3ffbebc74a7c295f, 0x408f400000000000, 0xc023d21aa11a5830, 0xc000f5b5a7d909c0,
            0x4078e241424554fa, 0x0000000000000001, 0x3fdb78d78bf36554, 0x55bee184c924210c,
            0x0000000000000015, 0x6cf896ea3b0e70be,
        ],
    ),
    (
        "flat/euler/1500m-25C-830hPa/rh50/no-zones",
        &[
            0x408f400000000000, 0x3fc0a49e7b99b739, 0x3ffbebdd724e55b5, 0x4078f09bdaea54ad,
            0x408c36bc18c31d35, 0x4075b318f42cf8e0, 0x0000000000000001, 0x408db55490505c15,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d3,
            0x3f50624dd2f1a9fc, 0x3fea529a4a656351, 0x3f5acb34cf863b98, 0xbec0ebaecef0e832,
            0x4089b4ae00f1a3cd, 0x0000000000000001, 0x3fd33b66fe995b33, 0x3fdbe76c8b439587,
            0x407434f9706af68d, 0xbfccfa490219ecc3, 0xbfc6186508ca8c43, 0x4084e685f74ec2d8,
            0x0000000000000001, 0x3fd5031266555170, 0x3febef9db22d0e5c, 0x40826c3e4e59f52c,
            0xc000f0acd8303890, 0xbfe40e14a2a7f550, 0x40815968c38576b6, 0x0000000000000001,
            0x3fd7092efffbb713, 0x3ff4f5c28f5c285f, 0x40896142409eba87, 0xc015ac58ab71a454,
            0xbff4aedc01c4334f, 0x407d398b0ff8f824, 0x0000000000000001, 0x3fd9700bec7e0728,
            0x3ffbebdd724e55b5, 0x408f400000000000, 0xc023d80b47bb233e, 0xc000f5e1f77d6239,
            0x4078f09bdaea54ad, 0x0000000000000001, 0x3fdb7afbd7a5053f, 0x9c6275f1f8a828be,
            0x0000000000000015, 0xef9f6d9a10279fd7,
        ],
    ),
    (
        "flat/euler/1500m-25C-830hPa/rh50/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc0a5dc4b34a33b, 0x3ffbebc74a7c295f, 0x4078e241424554fa,
            0x408c164b80128175, 0x4075b318f42cf8e0, 0x0000000000000001, 0x408da2437d92d752,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d3,
            0x3f50624dd2f1a9fc, 0x3fea529f2159466f, 0x3f5acb39c436ada5, 0xbec0d398a2a9eea0,
            0x4089b4b2badc6d3b, 0x0000000000000001, 0x3fd33b655d058486, 0x3fdbe76c8b439587,
            0x407437d3ddce2dfa, 0xbfccf8af763078a8, 0xbfc5fbfe78f7c96a, 0x4084ec66fdef8fe1,
            0x0000000000000001, 0x3fd500b126fc82ad, 0x3febef9db22d0e5c, 0x40826f9b18cd97fd,
            0xc000efee7e1db1b7, 0xbfe3fd588b8d1685, 0x40815a5ef5a3e066, 0x0000000000000001,
            0x3fd70877b5358f22, 0x3ff4f5c28f5c285f, 0x408963d6f7c01710, 0xc015a853282f0aa3,
            0xbff4a86fd78a6fdf, 0x407d3158348cb33e, 0x0000000000000001, 0x3fd974026bf3db51,
            0x3ffbebc74a7c295f, 0x408f400000000000, 0xc023d21aa11a5830, 0xc000f5b5a7d909c0,
            0x4078e241424554fa, 0x0000000000000001, 0x3fdb813b22872bcf, 0x55bee184c924210c,
            0x0000000000000015, 0x6cf896ea3b0e70be,
        ],
    ),
    (
        "flat/euler/1500m-25C-830hPa/rh100/no-zones",
        &[
            0x408f400000000000, 0x3fc0a5ed3ba10f19, 0x3ffbd9bbb3be977b, 0x407910f0865a688e,
            0x408c8011dab0f1c2, 0x4075c3b860e9c586, 0x0000000000000001, 0x408dc57fcd6f93e1,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006ce,
            0x3f50624dd2f1a9fc, 0x3fea529f6372da24, 0x3f5acb3a07e66570, 0xbec0d24fa9443d64,
            0x4089b4b2fb696fe3, 0x0000000000000001, 0x3fd34230642e05f3, 0x3fdbd70a3d70a3dd,
            0x40742d47626fce8e, 0xbfcca59faeb97cbd, 0xbfc5e1ec39ab61f8, 0x4084ef25d2cbc62c,
            0x0000000000000001, 0x3fd50610d79c7c9b, 0x3febdf3b645a1cb2, 0x408268531f8ab9dc,
            0xc000d9da38cf00b8, 0xbfe3e00c2f1162e1, 0x4081664cb853022a, 0x0000000000000001,
            0x3fd7097b053a731e, 0x3ff4e560418936b7, 0x40895ba2dfc51c98, 0xc0158b32d9cb1ae6,
            0xbff47b50a1015d19, 0x407d5a43a92009dd, 0x0000000000000001, 0x3fd96b0da12b8868,
            0x3ffbd9bbb3be977b, 0x408f400000000000, 0xc023c20bee77aeb5, 0xc000d19e7a5de5e0,
            0x407910f0865a688e, 0x0000000000000001, 0x3fdb752ed089b847, 0xc9ae979a164a2091,
            0x0000000000000015, 0xe83ab5908a642098,
        ],
    ),
    (
        "flat/euler/1500m-25C-830hPa/rh100/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc0a5dc4b34a33b, 0x3ffbebc74a7c295f, 0x4078e241424554fa,
            0x408c164b80128175, 0x4075c3b860e9c586, 0x0000000000000001, 0x408d8793f7ab8b60,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d3,
            0x3f50624dd2f1a9fc, 0x3fea529f2159466f, 0x3f5acb39c436ada5, 0xbec0d398a2a9eea0,
            0x4089b4b2badc6d3b, 0x0000000000000001, 0x3fd342307a83c9c7, 0x3fdbe76c8b439587,
            0x407437d3ddce2dfa, 0xbfccf8af763078a8, 0xbfc5fbfe78f7c96a, 0x4084ec66fdef8fe1,
            0x0000000000000001, 0x3fd5072ebcc7a6ed, 0x3febef9db22d0e5c, 0x40826f9b18cd97fd,
            0xc000efee7e1db1b7, 0xbfe3fd588b8d1685, 0x40815a5ef5a3e066, 0x0000000000000001,
            0x3fd7126632c1beb3, 0x3ff4f5c28f5c285f, 0x408963d6f7c01710, 0xc015a853282f0aa3,
            0xbff4a86fd78a6fdf, 0x407d3158348cb33e, 0x0000000000000001, 0x3fd97ebe9d617cc4,
            0x3ffbebc74a7c295f, 0x408f400000000000, 0xc023d21aa11a5830, 0xc000f5b5a7d909c0,
            0x4078e241424554fa, 0x0000000000000001, 0x3fdb895f2d6ec93c, 0x55bee184c924210c,
            0x0000000000000015, 0x6cf896ea3b0e70be,
        ],
    ),
    (
        "flat/rk4/icao-sea-level/rh0/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc08a6b00b0f847, 0x3fff64812240da88, 0x40742948bcb96758,
            0x4082701269538b02, 0x407544b4bbbd89be, 0x0000000000000001, 0x4088336671679d60,
            0x0000000000000001, 0x408d3fa7b25ff175, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ac,
            0x3f50624dd2f1a9fc, 0x3fea53e15f4bb302, 0x3f5ae115c0abe601, 0xbeb51fa51d272d29,
            0x4089b3db0cdf2f77, 0x0000000000000001, 0x3fd30dc2baedb34a, 0x3fdf6c8b4395810d,
            0x4075e389b39d46bf, 0xbfd81fb74b9f5616, 0xbfd0c1cba7966cbe, 0x408356a1eb3e60ee,
            0x0000000000000001, 0x3fd58ca7f58a1e17, 0x3fef6c8b4395810d, 0x408352e6d8f5f827,
            0xc0060439df6053b2, 0xbfed85d2d2ad3d7d, 0x407e17488423ef91, 0x0000000000000001,
            0x3fd8b888c964985f, 0x3ff7916872b01fdd, 0x4089eb9dcc63ac24, 0xc01aeed49501302a,
            0xbffdc5275ecac1fb, 0x4078007c3a13da39, 0x0000000000000001, 0x3fdbaac0e7821e92,
            0x3fff64812240da88, 0x408f400000000000, 0xc0282c308f3b2b28, 0xc007e72957626ce6,
            0x40742948bcb96758, 0x0000000000000001, 0x3fcbd0fa091851a0, 0x4cd2dc22107de48c,
            0x0000000000000015, 0x073d879945d3b354,
        ],
    ),
    (
        "flat/rk4/icao-sea-level/rh50/no-zones",
        &[
            0x408f400000000000, 0x3fc0888f199305b8, 0x3fff64e33a7f7a13, 0x40743e62a1af51c0,
            0x408296bed5e82875, 0x40754bf73b88949e, 0x0000000000000001, 0x40882eafe01a25b6,
            0x0000000000000001, 0x408d4879e413f238, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ac,
            0x3f50624dd2f1a9fc, 0x3fea53dd9ed231b2, 0x3f5ae111eb18f016, 0xbeb544fd88203852,
            0x4089b3d3b98319b8, 0x0000000000000001, 0x3fd310d7a9b4cff1, 0x3fdf6c8b4395810d,
            0x4075de59fe463a7d, 0xbfd81ffeb5bf0732, 0xbfd0db9b8ecd39c7, 0x40834d826cbd06c8,
            0x0000000000000001, 0x3fd594c6fc22d726, 0x3fef6c8b4395810d, 0x40834d89f92716e0,
            0xc00606197d43c95a, 0xbfeda08383697eb4, 0x407e15e4981cd1d0, 0x0000000000000001,
            0x3fd8bdfb9a2b745d, 0x3ff7916872b01fdd, 0x4089e7909a72f915, 0xc01af691b0b6b64c,
            0xbffdcf3c64022e3b, 0x40780b59555895c1, 0x0000000000000001, 0x3fdba9b2a935b110,
            0x3fff64e33a7f7a13, 0x408f400000000000, 0xc02837e14888a07f, 0xc007e7ed87dfac02,
            0x40743e62a1af51c0, 0x0000000000000001, 0x3fccbef559d3181e, 0x7eb3d9d83b87d0c9,
            0x0000000000000015, 0x5ad8cbc328fd5fc7,
        ],
    ),
    (
        "flat/rk4/icao-sea-level/rh50/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc08a6b00b0f847, 0x3fff64812240da88, 0x40742948bcb96758,
            0x4082701269538b02, 0x40754bf73b88949e, 0x0000000000000001, 0x408829987e9e0ec8,
            0x0000000000000001, 0x408d34c47598eed3, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ac,
            0x3f50624dd2f1a9fc, 0x3fea53e15f4bb302, 0x3f5ae115c0abe601, 0xbeb51fa51d272d29,
            0x4089b3db0cdf2f77, 0x0000000000000001, 0x3fd310d51902f0dc, 0x3fdf6c8b4395810d,
            0x4075e389b39d46bf, 0xbfd81fb74b9f5616, 0xbfd0c1cba7966cbe, 0x408356a1eb3e60ee,
            0x0000000000000001, 0x3fd5900dbfdd90b1, 0x3fef6c8b4395810d, 0x408352e6d8f5f827,
            0xc0060439df6053b2, 0xbfed85d2d2ad3d7d, 0x407e17488423ef91, 0x0000000000000001,
            0x3fd8bd5532cb817a, 0x3ff7916872b01fdd, 0x4089eb9dcc63ac24, 0xc01aeed49501302a,
            0xbffdc5275ecac1fb, 0x4078007c3a13da39, 0x0000000000000001, 0x3fdbadf83078d135,
            0x3fff64812240da88, 0x408f400000000000, 0xc0282c308f3b2b28, 0xc007e72957626ce6,
            0x40742948bcb96758, 0x0000000000000001, 0x3fcb696587d2e251, 0x4cd2dc22107de48c,
            0x0000000000000015, 0x073d879945d3b354,
        ],
    ),
    (
        "flat/rk4/icao-sea-level/rh100/no-zones",
        &[
            0x408f400000000000, 0x3fc089403ae45153, 0x3fff58c7c423eb94, 0x407449a8b9d107c7,
            0x4082ab78ef7327d9, 0x40755334ef1fc968, 0x0000000000000001, 0x408836a0f791e35b,
            0x0000000000000001, 0x408d520575f1105a, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007a9,
            0x3f50624dd2f1a9fc, 0x3fea53df04067223, 0x3f5ae11358260256, 0xbeb53719fa0055cd,
            0x4089b3d672f8da47, 0x0000000000000001, 0x3fd313e549c9ed28, 0x3fdf5c28f5c28f63,
            0x4075d664fb2dcab2, 0xbfd7ef59454044bd, 0xbfd0c1ac6542479b, 0x408353934db01afb,
            0x0000000000000001, 0x3fd5950908e1af3c, 0x3fef5c28f5c28f63, 0x408348cfb6a8dec4,
            0xc005ec91648a9cda, 0xbfed7681aaed6281, 0x407e26835acef243, 0x0000000000000001,
            0x3fd8bafff4469ba1, 0x3ff7851eb851ea9f, 0x4089e3aa296d2883, 0xc01adbbabdf3d3c6,
            0xbffda7ca0feda737, 0x40781d2594d57df3, 0x0000000000000001, 0x3fdba5dedd6b4514,
            0x3fff58c7c423eb94, 0x408f400000000000, 0xc02827939058c461, 0xc007cfb69b288f19,
            0x407449a8b9d107c7, 0x0000000000000001, 0x3fcd1060378a604b, 0xe4ce7350fb10a8ab,
            0x0000000000000015, 0xfc69b6ec549ff149,
        ],
    ),
    (
        "flat/rk4/icao-sea-level/rh100/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc08a6b00b0f847, 0x3fff64812240da88, 0x40742948bcb96758,
            0x4082701269538b02, 0x40755334ef1fc968, 0x0000000000000001, 0x4088230d4f922c84,
            0x0000000000000001, 0x408d29dd077457b5, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ac,
            0x3f50624dd2f1a9fc, 0x3fea53e15f4bb302, 0x3f5ae115c0abe601, 0xbeb51fa51d272d29,
            0x4089b3db0cdf2f77, 0x0000000000000001, 0x3fd313e3ad9ac072, 0x3fdf6c8b4395810d,
            0x4075e389b39d46bf, 0xbfd81fb74b9f5616, 0xbfd0c1cba7966cbe, 0x408356a1eb3e60ee,
            0x0000000000000001, 0x3fd59373439e3474, 0x3fef6c8b4395810d, 0x408352e6d8f5f827,
            0xc0060439df6053b2, 0xbfed85d2d2ad3d7d, 0x407e17488423ef91, 0x0000000000000001,
            0x3fd8c21ccfc5ddfc, 0x3ff7916872b01fdd, 0x4089eb9dcc63ac24, 0xc01aeed49501302a,
            0xbffdc5275ecac1fb, 0x4078007c3a13da39, 0x0000000000000001, 0x3fdbb1229897387d,
            0x3fff64812240da88, 0x408f400000000000, 0xc0282c308f3b2b28, 0xc007e72957626ce6,
            0x40742948bcb96758, 0x0000000000000001, 0x3fcb0924251d7777, 0x4cd2dc22107de48c,
            0x0000000000000015, 0x073d879945d3b354,
        ],
    ),
    (
        "flat/rk4/1500m-25C-830hPa/rh0/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc0c1e16e61ad48, 0x3ffbe88eeca12bfa, 0x4078e535fc62e906,
            0x408c1cf82d4940f4, 0x4075a2608de6e31f, 0x0000000000000001, 0x408dc28ebe7353e8,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d2,
            0x3f50624dd2f1a9fc, 0x3fea544fff9bee4f, 0x3f5ae186cf1e4d86, 0xbeb0d27d096d8882,
            0x4089b4b30e53aec3, 0x0000000000000001, 0x3fd3348d35b3d8ec, 0x3fdbe76c8b439587,
            0x4074393ee18d2105, 0xbfccad387e292bba, 0xbfc5ede1148db8e0, 0x4084ecca29a8d823,
            0x0000000000000001, 0x3fd4fa07bc1bff43, 0x3febef9db22d0e5c, 0x408270f98262bb9e,
            0xc000e69ccd056c86, 0xbfe3f68891b66d88, 0x40815ae611a949cb, 0x0000000000000001,
            0x3fd6fe2aa4d4ba7c, 0x3ff4f1a9fbe76bf5, 0x408962121ceca0c8, 0xc015983cef021717,
            0xbff49c747327554b, 0x407d355fdfb536be, 0x0000000000000001, 0x3fd9671d7433945a,
            0x3ffbe88eeca12bfa, 0x408f400000000000, 0xc023c8eceae3239c, 0xc000ef44ec230eea,
            0x4078e535fc62e906, 0x0000000000000001, 0x3fdb778990cc7f57, 0x24fbd41142cba0c4,
            0x0000000000000015, 0x898ef9b4cdb30a0e,
        ],
    ),
    (
        "flat/rk4/1500m-25C-830hPa/rh50/no-zones",
        &[
            0x408f400000000000, 0x3fc0c0a6a7b622a1, 0x3ffbe8a0960a356d, 0x4078f3967dc77655,
            0x408c3d79ffa0a91c, 0x4075b318f42cf8e0, 0x0000000000000001, 0x408db7b5afc5e6c6,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d2,
            0x3f50624dd2f1a9fc, 0x3fea544d94744a92, 0x3f5ae184565cb185, 0xbeb0ea900276bcff,
            0x4089b4ae555a1abb, 0x0000000000000001, 0x3fd33b66e177a249, 0x3fdbe76c8b439587,
            0x407436666e4f39dd, 0xbfccaecaefbf6cff, 0xbfc60a33f6477fcf, 0x4084e6ea08d94ff4,
            0x0000000000000001, 0x3fd502e9d251f699, 0x3febef9db22d0e5c, 0x40826d9f3e09f0bb,
            0xc000e75b4ea009ee, 0xbfe407381933f895, 0x408159f55d3cf3bf, 0x0000000000000001,
            0x3fd708c64cb1e201, 0x3ff4f1a9fbe76bf5, 0x40895f803c623e55, 0xc0159c3f2a7ffa62,
            0xbff4a2d98d49d55e, 0x407d3d95e45e2405, 0x0000000000000001, 0x3fd96e16dc653d47,
            0x3ffbe8a0960a356d, 0x408f400000000000, 0xc023ced4665df114, 0xc000ef683ef521d9,
            0x4078f3967dc77655, 0x0000000000000001, 0x3fdb79ad3f2a4e2e, 0xe67c1cb08a40aba4,
            0x0000000000000015, 0x58594b1c2ef2578a,
        ],
    ),
    (
        "flat/rk4/1500m-25C-830hPa/rh50/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc0c1e16e61ad48, 0x3ffbe88eeca12bfa, 0x4078e535fc62e906,
            0x408c1cf82d4940f4, 0x4075b318f42cf8e0, 0x0000000000000001, 0x408da7f2f494f19b,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d2,
            0x3f50624dd2f1a9fc, 0x3fea544fff9bee4f, 0x3f5ae186cf1e4d86, 0xbeb0d27d096d8882,
            0x4089b4b30e53aec3, 0x0000000000000001, 0x3fd33b6540371308, 0x3fdbe76c8b439587,
            0x4074393ee18d2105, 0xbfccad387e292bba, 0xbfc5ede1148db8e0, 0x4084ecca29a8d823,
            0x0000000000000001, 0x3fd50089155a017f, 0x3febef9db22d0e5c, 0x408270f98262bb9e,
            0xc000e69ccd056c86, 0xbfe3f68891b66d88, 0x40815ae611a949cb, 0x0000000000000001,
            0x3fd70813339cc23a, 0x3ff4f1a9fbe76bf5, 0x408962121ceca0c8, 0xc015983cef021717,
            0xbff49c747327554b, 0x407d355fdfb536be, 0x0000000000000001, 0x3fd972100b9bd822,
            0x3ffbe88eeca12bfa, 0x408f400000000000, 0xc023c8eceae3239c, 0xc000ef44ec230eea,
            0x4078e535fc62e906, 0x0000000000000001, 0x3fdb7ff38b1cecb4, 0x24fbd41142cba0c4,
            0x0000000000000015, 0x898ef9b4cdb30a0e,
        ],
    ),
    (
        "flat/rk4/1500m-25C-830hPa/rh100/no-zones",
        &[
            0x408f400000000000, 0x3fc0c1f23553721f, 0x3ffbd6888827ed32, 0x407913e5211b3b0b,
            0x408c86cac5778ed4, 0x4075c3b860e9c586, 0x0000000000000001, 0x408dc7de50dae6b5,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006cd,
            0x3f50624dd2f1a9fc, 0x3fea545020a2f883, 0x3f5ae186f0df1515, 0xbeb0d13449406f25,
            0x4089b4b34ecfc51d, 0x0000000000000001, 0x3fd342304752aabf, 0x3fdbd70a3d70a3dd,
            0x40742eb17b5e1b0e, 0xbfcc5a54e5b0ece3, 0xbfc5d3d7f6917791, 0x4084ef88cfb18a98,
            0x0000000000000001, 0x3fd505e89f23a3fa, 0x3febd70a3d70a3dd, 0x4082653cd09233c0,
            0xc000c44359696014, 0xbfe3cea4e17f5634, 0x408168a6e85838e8, 0x0000000000000001,
            0x3fd707bc2cf2c353, 0x3ff4e147ae147a4d, 0x408959d8c59dd081, 0xc0157b2327340f33,
            0xbff46f624b6ff823, 0x407d5e4cd5ce97c5, 0x0000000000000001, 0x3fd969194a5e5f23,
            0x3ffbd6888827ed32, 0x408f400000000000, 0xc023b8e713f2cfcf, 0xc000cb382330913f,
            0x407913e5211b3b0b, 0x0000000000000001, 0x3fdb73e1fd500fd9, 0x6e5767ccc5f4d631,
            0x0000000000000015, 0x4099b062ca3bb3ee,
        ],
    ),
    (
        "flat/rk4/1500m-25C-830hPa/rh100/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc0c1e16e61ad48, 0x3ffbe88eeca12bfa, 0x4078e535fc62e906,
            0x408c1cf82d4940f4, 0x4075c3b860e9c586, 0x0000000000000001, 0x408d89ec1c519629,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d2,
            0x3f50624dd2f1a9fc, 0x3fea544fff9bee4f, 0x3f5ae186cf1e4d86, 0xbeb0d27d096d8882,
            0x4089b4b30e53aec3, 0x0000000000000001, 0x3fd342305da292f2, 0x3fdbe76c8b439587,
            0x4074393ee18d2105, 0xbfccad387e292bba, 0xbfc5ede1148db8e0, 0x4084ecca29a8d823,
            0x0000000000000001, 0x3fd507065ccb764f, 0x3febef9db22d0e5c, 0x408270f98262bb9e,
            0xc000e69ccd056c86, 0xbfe3f68891b66d88, 0x40815ae611a949cb, 0x0000000000000001,
            0x3fd7120098ab2f94, 0x3ff4f1a9fbe76bf5, 0x408962121ceca0c8, 0xc015983cef021717,
            0xbff49c747327554b, 0x407d355fdfb536be, 0x0000000000000001, 0x3fd97cd06c20114a,
            0x3ffbe88eeca12bfa, 0x408f400000000000, 0xc023c8eceae3239c, 0xc000ef44ec230eea,
            0x4078e535fc62e906, 0x0000000000000001, 0x3fdb881e8ae8f2d6, 0x24fbd41142cba0c4,
            0x0000000000000015, 0x898ef9b4cdb30a0e,
        ],
    ),
    (
        "flat/rk45/icao-sea-level/rh0/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc089ab6178349e, 0x3fff64868865b520, 0x40742946834630fa,
            0x4082700e57c93494, 0x407544b4bbbd89be, 0x0000000000000001, 0x40884d0e32fe4d99,
            0x0000000000000001, 0x408d4a2ab2480d21, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000c9,
            0x3f50624dd2f1a9fc, 0x3fea53e15f4bb261, 0x3f5ae115c0abe5c9, 0xbeb51fa51d2d6e35,
            0x4089b3db0cdf2f77, 0x0000000000000001, 0x3fd30dc2baedb34a, 0x3fde33269df97ab0,
            0x4075252c6a4a7b8f, 0xbfd4963b31f97b7f, 0xbfcf1f4dd4292d4e, 0x408389c145f947fa,
            0x0000000000000001, 0x3fd572cfdb164ca6, 0x3fef19934efcbd5c, 0x40832bcb79e9d0f8,
            0xc0057c838a699351, 0xbfecfc94d4fb7047, 0x407e3d4695855b60, 0x0000000000000001,
            0x3fd8a6d107df1689, 0x3ff78cc9a77e5eb0, 0x4089e8242983983e, 0xc01ae3c8256ca92b,
            0xbffdbb52056c6334, 0x40780379b92d181e, 0x0000000000000001, 0x3fdba991e264a04d,
            0x3fff64868865b520, 0x408f400000000000, 0xc0282c38ec791fbd, 0xc007e73423ac1e4c,
            0x40742946834630fa, 0x0000000000000001, 0x3fcbd0d75acdc80a, 0x8a2a187cfa360d1c,
            0x0000000000000015, 0x45f380cbc4557ef0,
        ],
    ),
    (
        "flat/rk45/icao-sea-level/rh50/no-zones",
        &[
            0x408f400000000000, 0x3fc087c77972058e, 0x3fff64e377e9db29, 0x40743e635a393a86,
            0x408296c028d050d0, 0x40754bf73b88949e, 0x0000000000000001, 0x4088486225051459,
            0x0000000000000001, 0x408d4820a195189b, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000c9,
            0x3f50624dd2f1a9fc, 0x3fea53dd9ed23105, 0x3f5ae111eb18efc6, 0xbeb544fd8826e82f,
            0x4089b3d3b98319b6, 0x0000000000000001, 0x3fd310d7a9b4cff1, 0x3fde33269df97ab0,
            0x40752055071fd938, 0xbfd496dd265f688c, 0xbfcf4f7e805f8229, 0x408380d7db1947b1,
            0x0000000000000001, 0x3fd57a977e4b5e2e, 0x3fef19934efcbd5c, 0x40832671f3cc289d,
            0xc0057e3f33426ac7, 0xbfed1734d90c6e36, 0x407e3b90f557bf5b, 0x0000000000000001,
            0x3fd8ac691e381f59, 0x3ff78cc9a77e5eb0, 0x4089e417f15be877, 0xc01aeb812f5afb85,
            0xbffdc5649d0f15f5, 0x40780e559ccef35f, 0x0000000000000001, 0x3fdba883844ef71c,
            0x3fff64e377e9db29, 0x408f400000000000, 0xc02837e2d528a29f, 0xc007e7ee02b46a5b,
            0x40743e635a393a86, 0x0000000000000001, 0x3fccbf0262fb2975, 0x875d54bb36b95529,
            0x0000000000000015, 0xd2f2d8cf184133f8,
        ],
    ),
    (
        "flat/rk45/icao-sea-level/rh50/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc089ab6178349e, 0x3fff64868865b520, 0x40742946834630fa,
            0x4082700e57c93494, 0x40754bf73b88949e, 0x0000000000000001, 0x40882c7028e77988,
            0x0000000000000001, 0x408d4a2ab2480d21, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000c9,
            0x3f50624dd2f1a9fc, 0x3fea53e15f4bb261, 0x3f5ae115c0abe5c9, 0xbeb51fa51d2d6e35,
            0x4089b3db0cdf2f77, 0x0000000000000001, 0x3fd310d51902f0dc, 0x3fde33269df97ab0,
            0x4075252c6a4a7b8f, 0xbfd4963b31f97b7f, 0xbfcf1f4dd4292d4e, 0x408389c145f947fa,
            0x0000000000000001, 0x3fd57622360145b0, 0x3fef19934efcbd5c, 0x40832bcb79e9d0f8,
            0xc0057c838a699351, 0xbfecfc94d4fb7047, 0x407e3d4695855b60, 0x0000000000000001,
            0x3fd8ab9d81ea9999, 0x3ff78cc9a77e5eb0, 0x4089e8242983983e, 0xc01ae3c8256ca92b,
            0xbffdbb52056c6334, 0x40780379b92d181e, 0x0000000000000001, 0x3fdbacccb3ca2e2a,
            0x3fff64868865b520, 0x408f400000000000, 0xc0282c38ec791fbd, 0xc007e73423ac1e4c,
            0x40742946834630fa, 0x0000000000000001, 0x3fcb69452eab28b1, 0x8a2a187cfa360d1c,
            0x0000000000000015, 0xe1722d8e4ee2826c,
        ],
    ),
    (
        "flat/rk45/icao-sea-level/rh100/no-zones",
        &[
            0x408f400000000000, 0x3fc0887b957fadef, 0x3fff58c9b8dbb0f2, 0x407449af7bf1eb0e,
            0x4082ab855fc3b2eb, 0x40755334ef1fc968, 0x0000000000000001, 0x40884d16e9a847f7,
            0x0000000000000001, 0x408d6a2ecbff02b3, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000c8,
            0x3f50624dd2f1a9fc, 0x3fea53df0406717d, 0x3f5ae11358260211, 0xbeb53719fa06d43c,
            0x4089b3d672f8da46, 0x0000000000000001, 0x3fd313e549c9ed29, 0x3fde33269df97ab0,
            0x407522237aede357, 0xbfd496a2941c7151, 0xbfcf3d83523bc1d6, 0x4083842bd8882208,
            0x0000000000000001, 0x3fd57c3f75f5eadb, 0x3fef19934efcbd5c, 0x4083296291730095,
            0xc0057f9c96440470, 0xbfed089322e548af, 0x407e44bef4dbecde, 0x0000000000000001,
            0x3fd8acebf0880b25, 0x3ff78cc9a77e5eb0, 0x4089e97045ee04ac, 0xc01aee15bc35d72d,
            0xbffdb817935cf2ed, 0x4078182ff0a25a10, 0x0000000000000001, 0x3fdba7d94eba04cd,
            0x3fff58c9b8dbb0f2, 0x408f400000000000, 0xc02827a021647d7e, 0xc007cfba849815ee,
            0x407449af7bf1eb0e, 0x0000000000000001, 0x3fcd10e0c1ad3acc, 0x5cf3ccb629f9da4f,
            0x0000000000000015, 0x91d097b439b176a2,
        ],
    ),
    (
        "flat/rk45/icao-sea-level/rh100/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc089ab6178349e, 0x3fff64868865b520, 0x40742946834630fa,
            0x4082700e57c93494, 0x40755334ef1fc968, 0x0000000000000001, 0x40882c7028e77988,
            0x0000000000000001, 0x408d2ef4df087b9c, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000c9,
            0x3f50624dd2f1a9fc, 0x3fea53e15f4bb261, 0x3f5ae115c0abe5c9, 0xbeb51fa51d2d6e35,
            0x4089b3db0cdf2f77, 0x0000000000000001, 0x3fd313e3ad9ac072, 0x3fde33269df97ab0,
            0x4075252c6a4a7b8f, 0xbfd4963b31f97b7f, 0xbfcf1f4dd4292d4e, 0x408389c145f947fa,
            0x0000000000000001, 0x3fd5797343939a6d, 0x3fef19934efcbd5c, 0x40832bcb79e9d0f8,
            0xc0057c838a699351, 0xbfecfc94d4fb7047, 0x407e3d4695855b60, 0x0000000000000001,
            0x3fd8b065280414b3, 0x3ff78cc9a77e5eb0, 0x4089e8242983983e, 0xc01ae3c8256ca92b,
            0xbffdbb52056c6334, 0x40780379b92d181e, 0x0000000000000001, 0x3fdbaffaa6a628a1,
            0x3fff64868865b520, 0x408f400000000000, 0xc0282c38ec791fbd, 0xc007e73423ac1e4c,
            0x40742946834630fa, 0x0000000000000001, 0x3fcb0905f5597609, 0x8a2a187cfa360d1c,
            0x0000000000000015, 0xe1722d8e4ee2826c,
        ],
    ),
    (
        "flat/rk45/1500m-25C-830hPa/rh0/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc0c1d3bcaa6f05, 0x3ffbe89591fa49c3, 0x4078e5323e39f549,
            0x408c1cefb9472b9b, 0x4075a2608de6e31f, 0x0000000000000001, 0x408dd2b71a12a847,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000b3,
            0x3f50624dd2f1a9fc, 0x3fea544fff9bee16, 0x3f5ae186cf1e4d8d, 0xbeb0d27d096fc3c6,
            0x4089b4b30e53aec3, 0x0000000000000001, 0x3fd3348d35b3d8ec, 0x3fda5c1c6088d6d8,
            0x4073355a4c8b6e10, 0xbfc52e415be39120, 0xbfc3ae4ee96bc6c4, 0x408526e07a6c91a0,
            0x0000000000000001, 0x3fd4e2e0f1b5b8c7, 0x3feb947496aad1d6, 0x40823f6ba7cbba53,
            0xc0005e7ee9a3f4d3, 0xbfe38084aa7a3f76, 0x40816f433c331f3b, 0x0000000000000001,
            0x3fd6ef501f91e2c9, 0x3ff4fd6d7e889c20, 0x40896ccb0d3e02c0, 0xc015b2b34989bd03,
            0xbff4b0d3d0ca7411, 0x407d2d2b0aa1bf31, 0x0000000000000001, 0x3fd96b1d9db209f8,
            0x3ffbe89591fa49c3, 0x408f400000000000, 0xc023c8fc87d114c3, 0xc000ef5236d5478e,
            0x4078e5323e39f549, 0x0000000000000001, 0x3fdb778b37e32967, 0x1a36ec8e303acd3f,
            0x0000000000000015, 0x78bf622636fd312b,
        ],
    ),
    (
        "flat/rk45/1500m-25C-830hPa/rh50/no-zones",
        &[
            0x408f400000000000, 0x3fc0c0967c480088, 0x3ffbe8a3ea11dae1, 0x4078f397963aa872,
            0x408c3d7c7a746c8b, 0x4075b318f42cf8e0, 0x0000000000000001, 0x408dd1e5b1093dde,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000b3,
            0x3f50624dd2f1a9fc, 0x3fea544d94744a55, 0x3f5ae184565cb186, 0xbeb0ea9002791418,
            0x4089b4ae555a1abb, 0x0000000000000001, 0x3fd33b66e177a249, 0x3fda5c1c6088d6d8,
            0x407332c90419287a, 0xbfc5302ce2df48f3, 0xbfc3c7dd6809a701, 0x4085213b76fde2dd,
            0x0000000000000001, 0x3fd4eb942deb09c4, 0x3feb947496aad1d6, 0x40823c1525f9df50,
            0xc0005f27bc9f728a, 0xbfe391217b8b4508, 0x40816e2a44a34f3c, 0x0000000000000001,
            0x3fd6f9e1cfa2ee08, 0x3ff4fd6d7e889c20, 0x40896a3e07d7185e, 0xc015b6c0277be547,
            0xbff4b72cd56c3961, 0x407d357213aa9236, 0x0000000000000001, 0x3fd972073f38a43f,
            0x3ffbe8a3ea11dae1, 0x408f400000000000, 0xc023cee010aad1ee, 0xc000ef6ee70469c9,
            0x4078f397963aa872, 0x0000000000000001, 0x3fdb79acc3f448ad, 0x5a3a8163b6b7a841,
            0x0000000000000015, 0x47954add471380b7,
        ],
    ),
    (
        "flat/rk45/1500m-25C-830hPa/rh50/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc0c1d3bcaa6f05, 0x3ffbe89591fa49c3, 0x4078e5323e39f549,
            0x408c1cefb9472b9b, 0x4075b318f42cf8e0, 0x0000000000000001, 0x408db1804637b3a6,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000b3,
            0x3f50624dd2f1a9fc, 0x3fea544fff9bee16, 0x3f5ae186cf1e4d8d, 0xbeb0d27d096fc3c6,
            0x4089b4b30e53aec3, 0x0000000000000001, 0x3fd33b6540371308, 0x3fda5c1c6088d6d8,
            0x4073355a4c8b6e10, 0xbfc52e415be39120, 0xbfc3ae4ee96bc6c4, 0x408526e07a6c91a0,
            0x0000000000000001, 0x3fd4e95682fe56bd, 0x3feb947496aad1d6, 0x40823f6ba7cbba53,
            0xc0005e7ee9a3f4d3, 0xbfe38084aa7a3f76, 0x40816f433c331f3b, 0x0000000000000001,
            0x3fd6f914ce4a7182, 0x3ff4fd6d7e889c20, 0x40896ccb0d3e02c0, 0xc015b2b34989bd03,
            0xbff4b0d3d0ca7411, 0x407d2d2b0aa1bf31, 0x0000000000000001, 0x3fd97606587fe935,
            0x3ffbe89591fa49c3, 0x408f400000000000, 0xc023c8fc87d114c3, 0xc000ef5236d5478e,
            0x4078e5323e39f549, 0x0000000000000001, 0x3fdb7ff52a7beb57, 0x1a36ec8e303acd3f,
            0x0000000000000015, 0x78bf622636fd312b,
        ],
    ),
    (
        "flat/rk45/1500m-25C-830hPa/rh100/no-zones",
        &[
            0x408f400000000000, 0x3fc0c1e4a5ad63ac, 0x3ffbd68b6ed88a51, 0x407913e614620aed,
            0x408c86cceeef0774, 0x4075c3b860e9c586, 0x0000000000000001, 0x408ddecd420fe31f,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000b2,
            0x3f50624dd2f1a9fc, 0x3fea545020a2f84a, 0x3f5ae186f0df1521, 0xbeb0d13449429eb1,
            0x4089b4b34ecfc51c, 0x0000000000000001, 0x3fd342304752aac0, 0x3fda5c1c6088d6d8,
            0x4073358006dfa4b3, 0xbfc52e2c8379c6c9, 0xbfc3acd75fe382dd, 0x4085273802d4b20e,
            0x0000000000000001, 0x3fd4ef9e66b30b0c, 0x3feb947496aad1d6, 0x408240f4cbb276f7,
            0xc00060e5c6e0b36b, 0xbfe378e009761f20, 0x40817762c4b8bf88, 0x0000000000000001,
            0x3fd6fce516bb3087, 0x3ff4d477bbf93ff7, 0x40894e155b0f0cec, 0xc0155e62f745b926,
            0xbff459684e613df4, 0x407d672fba862d4a, 0x0000000000000001, 0x3fd964c8c9f6cd0b,
            0x3ffbd68b6ed88a51, 0x408f400000000000, 0xc023b8f14535e0a0, 0xc000cb3df091c8af,
            0x407913e614620aed, 0x0000000000000001, 0x3fdb73e1923f1a68, 0x01c76564a332faf7,
            0x0000000000000015, 0x08bb0bff88af7800,
        ],
    ),
    (
        "flat/rk45/1500m-25C-830hPa/rh100/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc0c1d3bcaa6f05, 0x3ffbe89591fa49c3, 0x4078e5323e39f549,
            0x408c1cefb9472b9b, 0x4075c3b860e9c586, 0x0000000000000001, 0x408d902aacc6c785,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000b3,
            0x3f50624dd2f1a9fc, 0x3fea544fff9bee16, 0x3f5ae186cf1e4d8d, 0xbeb0d27d096fc3c6,
            0x4089b4b30e53aec3, 0x0000000000000001, 0x3fd342305da292f2, 0x3fda5c1c6088d6d8,
            0x4073355a4c8b6e10, 0xbfc52e415be39120, 0xbfc3ae4ee96bc6c4, 0x408526e07a6c91a0,
            0x0000000000000001, 0x3fd4efc118b9d3c8, 0x3feb947496aad1d6, 0x40823f6ba7cbba53,
            0xc0005e7ee9a3f4d3, 0xbfe38084aa7a3f76, 0x40816f433c331f3b, 0x0000000000000001,
            0x3fd702da33dcab5e, 0x3ff4fd6d7e889c20, 0x40896ccb0d3e02c0, 0xc015b2b34989bd03,
            0xbff4b0d3d0ca7411, 0x407d2d2b0aa1bf31, 0x0000000000000001, 0x3fd980be53c27608,
            0x3ffbe89591fa49c3, 0x408f400000000000, 0xc023c8fc87d114c3, 0xc000ef5236d5478e,
            0x4078e5323e39f549, 0x0000000000000001, 0x3fdb8820217ca3de, 0x1a36ec8e303acd3f,
            0x0000000000000015, 0x78bf622636fd312b,
        ],
    ),
];

/// `TrajectorySolver::solve`, the 30-degree lob to ground impact, Euler, RK4 and RK45 at both
/// stations, RH 0, no zones.
///
/// IMMUTABLE: no later step may change it, S8 included.
#[rustfmt::skip]
const TS_LOB_DRY: Pinned = &[
    (
        "lob30/euler/icao-sea-level/rh0/no-zones",
        &[
            0x40b2543ad8e1d6ab, 0x40919f5e9f3bf77b, 0x403ddcfac17feddb, 0x406190ebf688682a,
            0x405bfe4783a223d0, 0x407544b4bbbd89be, 0x0000000000000001, 0x408512249afa2f10,
            0x0000000000000001, 0x408993e29fbf5f96, 0x0000000000000001, 0x408e21138c4fa19d,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x00000000000074a9,
            0x3f50624dd2f1a9fc, 0x3fe6cb0ba004dbf9, 0x3fda5191c48293ad, 0xbec554a6c13e5dde,
            0x4089b3c67a315a4a, 0x0000000000000001, 0x3fd30dc9f145ee7f, 0x401ddd2f1a9fc21b,
            0x40a06651912f522d, 0x408f987602991214, 0xc032154c5840ce82, 0x4067c1de8d4823d1,
            0x0000000000000001, 0x3fc09339306582ca, 0x402ddd2f1a9fb83a, 0x40a94ce7ef9a7cd6,
            0x40910df7d7571a2c, 0xc044c6d2b5c18d88, 0x40607222bddebfd3, 0x0000000000000001,
            0x3fc093e8d9c4e12b, 0x403665e353f7d3ad, 0x40afd4e676b6c80f, 0x40851f6ecdd5d7c7,
            0xc050ae6f98def0f6, 0x405fcaaf9dd77217, 0x0000000000000001, 0x3fc09503b192c3d2,
            0x403ddcfac17feddb, 0x40b2543ad8e1d6ab, 0xc059000000000000, 0xc05747a279607d1b,
            0x406190ebf688682a, 0x0000000000000001, 0x3fc092c75f23d8de, 0x48d465da64d7819f,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/euler/1500m-25C-830hPa/rh0/no-zones",
        &[
            0x40b58a3ae51188ec, 0x40943b24a38eed21, 0x403fe358dbc3b8c5, 0x40636d14043edb58,
            0x40611e1d24a0f6e1, 0x4075a2608de6e31f, 0x0000000000000001, 0x4089e418db53b02c,
            0x0000000000000001, 0x408f981ac20e6ba8, 0x0000000000000001, 0x4092a428a9178510,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000007c92,
            0x3f50624dd2f1a9fc, 0x3fe6cbcba98febd3, 0x3fda526f835ed7b1, 0xbec1050afce375ad,
            0x4089b49f0694436f, 0x0000000000000001, 0x3fd334941ed79ca5, 0x401fe353f7cedd79,
            0x40a2b59c3b547d67, 0x4091ff6a76cf7dda, 0xc03272af3db3f664, 0x4069f73a9c009fbb,
            0x0000000000000001, 0x3fc096b12e98345b, 0x402fe3d70a3d692b, 0x40ad3c0d6681533c,
            0x4093a92da6213cb0, 0xc04563bc0605c25f, 0x40625e09f7d0e064, 0x0000000000000001,
            0x3fc092ad083f6a79, 0x4037eac083127576, 0x40b290287a082b9a, 0x4088bf4ff4809640,
            0xc0513fe1db7d9616, 0x4061a9cdc267c1a2, 0x0000000000000001, 0x3fc092e03b62e250,
            0x403fe358dbc3b8c5, 0x40b58a3ae51188ec, 0xc059000000000000, 0xc05826c20dfe6fb7,
            0x40636d14043edb58, 0x0000000000000001, 0x3fc0931592d51d9e, 0x0f0e25ac42d0da49,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk4/icao-sea-level/rh0/no-zones",
        &[
            0x40b255a415daa86f, 0x4091a159a9befc41, 0x403dde968cd612b0, 0x40619141cddb0935,
            0x405bff591dd87aa9, 0x407544b4bbbd89be, 0x0000000000000001, 0x40851722b81b66da,
            0x0000000000000001, 0x4089968c25825764, 0x0000000000000001, 0x408e2643e581560c,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x00000000000074af,
            0x3f50624dd2f1a9fc, 0x3fe6cce6e2c14197, 0x3fda53cb1f2d0932, 0xbeb552cdc5c1a9cd,
            0x4089b3c706356418, 0x0000000000000001, 0x3fd30dc9c0300e37, 0x401dde353f7cf136,
            0x40a067da1f522728, 0x408f9bfe623fcb05, 0xc03215387e96f59a, 0x4067c2575cd25631,
            0x0000000000000001, 0x3fc09339755c3fe7, 0x402ddeb851eb7ee1, 0x40a94f302d3b2cc2,
            0x40910fdbac53e01b, 0xc044c78a1465bdc0, 0x406072376d46626b, 0x0000000000000001,
            0x3fc093e8b424f728, 0x4036672b020c4e90, 0x40afd79c851385b6, 0x4085218736c85056,
            0xc050af3aa1e36908, 0x405fcb39660371f9, 0x0000000000000001, 0x3fc0950321cbfc35,
            0x403dde968cd612b0, 0x40b255a415daa86f, 0xc059000000000000, 0xc05748bc85098108,
            0x40619141cddb0935, 0x0000000000000001, 0x3fc092c73ca00942, 0xdc34b14797e3e04d,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk4/1500m-25C-830hPa/rh0/no-zones",
        &[
            0x40b58ba56a7ea7e0, 0x40943d12a723ab8f, 0x403fe4d194675d91, 0x40636d6725de2761,
            0x40611eafa73131fe, 0x4075a2608de6e31f, 0x0000000000000001, 0x4089e62db83a1696,
            0x0000000000000001, 0x408f9d2b38307381, 0x0000000000000001, 0x4092a6c14ccc961c,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000007c97,
            0x3f50624dd2f1a9fc, 0x3fe6cd46d9025355, 0x3fda5439ed9e1503, 0xbeb103d8dd425a35,
            0x4089b49f6152e32b, 0x0000000000000001, 0x3fd33493ff88b5b1, 0x401fe45a1cac0c94,
            0x40a2b71c781a5dd8, 0x4092012284988d11, 0xc03272a15dadf1b6, 0x4069f7a3493c6b04,
            0x0000000000000001, 0x3fc096b287a148aa, 0x402fe4dd2f1a9845, 0x40ad3e1252d1d2f0,
            0x4093ab1bbcbb3cbf, 0xc0456408800e39d8, 0x40625e367b266c81, 0x0000000000000001,
            0x3fc092ad07858da8, 0x4037ebc6a7efa492, 0x40b291713c9cf959, 0x4088c1cda736b9a0,
            0xc0514071eb05dee6, 0x4061aa0962fad32a, 0x0000000000000001, 0x3fc092e01a8498df,
            0x403fe4d194675d91, 0x40b58ba56a7ea7e0, 0xc059000000000000, 0xc05827b890fbbce8,
            0x40636d6725de2761, 0x0000000000000001, 0x3fc09315d460db5c, 0xb60e7fd79a5e2a02,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk45/icao-sea-level/rh0/no-zones",
        &[
            0x40b255a412c5f552, 0x4091a159a1462b16, 0x403dde9688243792, 0x40619141cb38f451,
            0x405bff591573e3fe, 0x407544b4bbbd89be, 0x0000000000000001, 0x40851c7d6b3bc133,
            0x0000000000000001, 0x40899b0913674280, 0x0000000000000001, 0x408e33074a8751c5,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000000baf,
            0x3f50624dd2f1a9fc, 0x3fe6cce6e2c14100, 0x3fda53cb1f2d0885, 0xbeb552cdc5c86cf6,
            0x4089b3c706356418, 0x0000000000000001, 0x3fd30dc9c0300e37, 0x401dc47a17f4120b,
            0x40a05ea90e83d73e, 0x408f917bf2242c1a, 0xc0320217ad90ec58, 0x4067cd5efb98abb0,
            0x0000000000000001, 0x3fc0934024f301cd, 0x402dd7ff9b5631a2, 0x40a94bccf7b699ac,
            0x40911134e057d6b6, 0xc044c208d945185d, 0x4060739be56e5439, 0x0000000000000001,
            0x3fc093e62ce29e0d, 0x403666e115592e5f, 0x40afd764cd3b2630, 0x4085224676f9fa94,
            0xc050aefab6d1281e, 0x405fcb216e0d3945, 0x0000000000000001, 0x3fc095033acf11e2,
            0x403dde9688243792, 0x40b255a412c5f552, 0xc059000000000000, 0xc05748bc8172fb64,
            0x40619141cb38f451, 0x0000000000000001, 0x3fc092c73ca117f2, 0xce891f3a0bb04800,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk45/1500m-25C-830hPa/rh0/no-zones",
        &[
            0x40b58ba566d68b59, 0x40943d12a2db344c, 0x403fe4d18fa6ba52, 0x40636d6723b28851,
            0x40611eafa35df0ce, 0x4075a2608de6e31f, 0x0000000000000001, 0x4089eba3a395ec09,
            0x0000000000000001, 0x408fa6a42052ff69, 0x0000000000000001, 0x4092ab0ea7d637d6,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000000c7a,
            0x3f50624dd2f1a9fc, 0x3fe6cd46d902531f, 0x3fda5439ed9e14c6, 0xbeb103d8dd44cad0,
            0x4089b49f6152e32b, 0x0000000000000001, 0x3fd33493ff88b5b1, 0x401fceb78897e909,
            0x40a2aeabde3ef102, 0x4091fc4577f702db, 0xc032630ede0c8554, 0x406a00742d03acaf,
            0x0000000000000001, 0x3fc096cfd497f04d, 0x402fe23d0bfa08a0, 0x40ad3c969cfb3007,
            0x4093aba39f6a37e9, 0xc04561f0bf1bc410, 0x40625ec41ef789bd, 0x0000000000000001,
            0x3fc092ad055d462c, 0x4037ebffcdab1a27, 0x40b29189dc7a99ab, 0x4088c130095656ff,
            0xc05140a2388bed8c, 0x4061aa11efa77f31, 0x0000000000000001, 0x3fc092e015cefff7,
            0x403fe4d18fa6ba52, 0x40b58ba566d68b59, 0xc059000000000000, 0xc05827b88d8b63d5,
            0x40636d6723b28851, 0x0000000000000001, 0x3fc09315d45f2501, 0x204f17f770e31887,
            0xffffffffffffffff,
        ],
    ),
];

/// `TrajectorySolver::solve`, the 30-degree lob: RH 50 and 100 without zones, and RH 50 under humid
/// zones.
///
/// MAY CHANGE ONLY IN S8 (the unified aloft vapour model, plan C7): byte-identical through every
/// other step.
#[rustfmt::skip]
const TS_LOB_HUMID: Pinned = &[
    (
        "lob30/euler/icao-sea-level/rh50/no-zones",
        &[
            0x40b2602b7ae19e36, 0x4091a943136f9e81, 0x403de545b4022f05, 0x40619767ba00af77,
            0x405c12f5509aba95, 0x40754bf73b88949e, 0x0000000000000001, 0x4085191a72aaeedc,
            0x0000000000000001, 0x40899c4b90990320, 0x0000000000000001, 0x408e2b16edd6a00e,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x00000000000074c9,
            0x3f50624dd2f1a9fc, 0x3fe6cb0e0b537889, 0x3fda51948f9fd0d9, 0xbec546bfb16f49c3,
            0x4089b3c9348be756, 0x0000000000000001, 0x3fd310db58b42b79, 0x401de56041893af3,
            0x40a06f2692b24765, 0x408fa8fd23ac4883, 0xc03217254a5a3d5d, 0x4067cb35c63e5831,
            0x0000000000000001, 0x3fc09339e4d696f6, 0x402de5604189310a, 0x40a95bef63e195ef,
            0x409117dee5c57b4a, 0xc044c99df11a3c0e, 0x406079467f4a775f, 0x0000000000000001,
            0x3fc093e615568cbb, 0x40366c083126ee55, 0x40afe8d4bf94829b, 0x40852d3e41507155,
            0xc050b10025beb01d, 0x405fd730c77f34e1, 0x0000000000000001, 0x3fc09501f8a07e3a,
            0x403de545b4022f05, 0x40b2602b7ae19e36, 0xc059000000000000, 0xc0574ba38bf7c66a,
            0x40619767ba00af77, 0x0000000000000001, 0x3fc092c72d105441, 0x80ee6152bc0cdd91,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/euler/icao-sea-level/rh100/no-zones",
        &[
            0x40b26c2380a3ba1e, 0x4091b32c57db7779, 0x403ded923cd753cd, 0x40619de6a8cb424e,
            0x405c27b4df594a50, 0x40755334ef1fc968, 0x0000000000000001, 0x40851d34976a860d,
            0x0000000000000001, 0x4089a4b7e2fc1d59, 0x0000000000000001, 0x408e374c115b40b3,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x00000000000074ea,
            0x3f50624dd2f1a9fc, 0x3fe6cb107662f65b, 0x3fda51975a742b82, 0xbec538da0c5e0036,
            0x4089b3cbee9f472f, 0x0000000000000001, 0x3fd313e8f7ee0281, 0x401ded916872b3cb,
            0x40a077ff115c4f0d, 0x408fb98b158abc72, 0xc03218fbbad9c33c, 0x4067d48faa719444,
            0x0000000000000001, 0x3fc0933a9fd9b86f, 0x402dee147ae14167, 0x40a96b40300a1e0e,
            0x409121b0c8471230, 0xc044ccd1e7004f2d, 0x40608052a9019011, 0x0000000000000001,
            0x3fc093e3792b3a33, 0x4036726e978d54c4, 0x40affcff3617c48f, 0x40853a6e357bd855,
            0xc050b3c76d26d47c, 0x405fe3cd3b5c3488, 0x0000000000000001, 0x3fc095001ce504fc,
            0x403ded923cd753cd, 0x40b26c2380a3ba1e, 0xc059000000000000, 0xc0574fa38d679c27,
            0x40619de6a8cb424e, 0x0000000000000001, 0x3fc092c6f876c43e, 0xeef116ee21c8e7c3,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/euler/icao-sea-level/rh50/zones-rh20-rh80",
        &[
            0x40b23ed0eacc0b7a, 0x409192c9a5eb8e07, 0x403dd06bc16d98fe, 0x406182de1b1bbfb0,
            0x405bd18e5df5ca75, 0x40754bf73b88949e, 0x0000000000000001, 0x40851c4a90b96932,
            0x0000000000000001, 0x408992a20d737b6e, 0x0000000000000001, 0x408de852a04d8e36,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000007478,
            0x3f50624dd2f1a9fc, 0x3fe6cb148bee4258, 0x3fda519c11be810c, 0xbec5216177e29fc7,
            0x4089b3d0899ddd90, 0x0000000000000001, 0x3fd310d8c768a2b1, 0x401dd0e560418cd7,
            0x40a05b9103766e2e, 0x408f87905f7ecfb1, 0xc03210bbedd8d51f, 0x4067ab1596c8957e,
            0x0000000000000001, 0x3fc0932a8b9184b0, 0x402dd0e560418302, 0x40a9355c9c739904,
            0x40910082eeb292ed, 0xc044c2fddd020867, 0x4060611d239fa609, 0x0000000000000001,
            0x3fc09413006b480d, 0x40365cac08312bb1, 0x40afb2f76cebe8d5, 0x40850a34ff03aaa0,
            0xc050ab506597518d, 0x405fafc24012f0e1, 0x0000000000000001, 0x3fc0952af0782bc1,
            0x403dd06bc16d98fe, 0x40b23ed0eacc0b7a, 0xc059000000000000, 0xc05742c47622da3b,
            0x406182de1b1bbfb0, 0x0000000000000001, 0x3fc092d00c763699, 0x70e6ed4c5fb268bf,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/euler/1500m-25C-830hPa/rh50/no-zones",
        &[
            0x40b5a9772cce6a80, 0x409454207ec714c7, 0x403ff70694f7121a, 0x40637d54b9e77a56,
            0x40613acdaffe0bfb, 0x4075b318f42cf8e0, 0x0000000000000001, 0x4089f23b9073e799,
            0x0000000000000001, 0x408fb28c3190db93, 0x0000000000000001, 0x4092b37d3ba7eeb9,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000007cde,
            0x3f50624dd2f1a9fc, 0x3fe6cbd01341ba63, 0x3fda52749bd4531e, 0xbec0ebaecef0e832,
            0x4089b4a40076fcf5, 0x0000000000000001, 0x3fd33b6a724ce681, 0x401ff6c8b4395c7a,
            0x40a2cbcafbd54d68, 0x4092140f0d027edd, 0xc03276349d173842, 0x406a0da6224b6eb5,
            0x0000000000000001, 0x3fc096b8ecd5153d, 0x402ff74bc6a7e819, 0x40ad627d34d6eb90,
            0x4093c23a1c0360cd, 0xc0456962850195ed, 0x40626fdc4eb802ce, 0x0000000000000001,
            0x3fc092acfde12265, 0x4037f958106254c5, 0x40b2a9f89c318a8d, 0x4088e2871829dcce,
            0xc0514533ea1bafa5, 0x4061b999ca86dcb8, 0x0000000000000001, 0x3fc092df0f0e1a52,
            0x403ff70694f7121a, 0x40b5a9772cce6a80, 0xc059000000000000, 0xc0582f37b1120ff0,
            0x40637d54b9e77a56, 0x0000000000000001, 0x3fc093168c9d0af5, 0x9f1ffe5bd1cabda2,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/euler/1500m-25C-830hPa/rh100/no-zones",
        &[
            0x40b5c8ead43c6de1, 0x40946d416b66a0d5, 0x40400561e3ab612d, 0x40638daef109448d,
            0x406157c37090f5be, 0x4075c3b860e9c586, 0x0000000000000001, 0x408a006d5cd6281f,
            0x0000000000000001, 0x408fcaab0a566a85, 0x0000000000000001, 0x4092c1c753f1c292,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000007d2c,
            0x3f50624dd2f1a9fc, 0x3fe6cbd47d77bbe3, 0x3fda5279b4e274ff, 0xbec0d24fa9443d64,
            0x4089b4a8faeec90c, 0x0000000000000001, 0x3fd34233da20e3ed, 0x402005a1cac08545,
            0x40a2e279fe2d3082, 0x40922908b8b1c821, 0xc0327a63fef804a5, 0x406a23be77072b87,
            0x0000000000000001, 0x3fc096bff7c8bd17, 0x403005a1cac07f55, 0x40ad896a2f14b7e2,
            0x4093db547a94da6d, 0xc0456f5f20f626dc, 0x406281aee282be9d, 0x0000000000000001,
            0x3fc092acfd85af79, 0x40380872b020cba2, 0x40b2c427c15518bf, 0x408904aa79c903ef,
            0xc0514ae6eae24e63, 0x4061c990d60643a6, 0x0000000000000001, 0x3fc092ddc6055915,
            0x40400561e3ab612d, 0x40b5c8ead43c6de1, 0xc059000000000000, 0xc05837a8e9b9f803,
            0x40638daef109448d, 0x0000000000000001, 0x3fc09317abf53344, 0x39447ee6565bed6a,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/euler/1500m-25C-830hPa/rh50/zones-rh20-rh80",
        &[
            0x40b58664b123a0d6, 0x40943c5b3f9e66ea, 0x403fe292425ec2f5, 0x4063692fa9fba50b,
            0x40611741ec6e3548, 0x4075b318f42cf8e0, 0x0000000000000001, 0x4089e965d1afabdc,
            0x0000000000000001, 0x408f9c8f21e682a8, 0x0000000000000001, 0x409288eb0514ff76,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000007c8f,
            0x3f50624dd2f1a9fc, 0x3fe6cbd444393069, 0x3fda527972c8d8a1, 0xbec0d398a2a9eea0,
            0x4089b4a8ba61c664, 0x0000000000000001, 0x3fd33b68d0b81e39, 0x401fe24dd2f1ae5e,
            0x40a2b6f046249a69, 0x409202272ff700d0, 0xc03270b4d365f2ac, 0x4069efa120dbfffe,
            0x0000000000000001, 0x3fc096588d130102, 0x402fe2d0e5603a11, 0x40ad39dc18de9cfb,
            0x4093a9cb3598a0d5, 0xc045637f7cb94fd7, 0x406258b37c6da49f, 0x0000000000000001,
            0x3fc092adeaf29f66, 0x4037ea3d70a3dde8, 0x40b28db6c75a20e7, 0x4088be06b53fb420,
            0xc051403f85afa0cf, 0x4061a634bd64e71b, 0x0000000000000001, 0x3fc092ea280a75e4,
            0x403fe292425ec2f5, 0x40b58664b123a0d6, 0xc059000000000000, 0xc058275c3dadf077,
            0x4063692fa9fba50b, 0x0000000000000001, 0x3fc093073f97e2c0, 0x193a31cd252b33a9,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk4/icao-sea-level/rh50/no-zones",
        &[
            0x40b26194c916e698, 0x4091ab3e00251ac5, 0x403de6e100590812, 0x406197bd83ac16f7,
            0x405c140724370dbf, 0x40754bf73b88949e, 0x0000000000000001, 0x40851b3b1b15134c,
            0x0000000000000001, 0x4089a15ace6ab383, 0x0000000000000001, 0x408e304711d5ea48,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x00000000000074cf,
            0x3f50624dd2f1a9fc, 0x3fe6cce8183b044f, 0x3fda53cc84865671, 0xbeb544e8c16dea0e,
            0x4089b3c9bff527f9, 0x0000000000000001, 0x3fd310db27e2d469, 0x401de66666666a0e,
            0x40a070af0cd7a4a6, 0x408fac8543031496, 0xc03217117f0c8de1, 0x4067cbae5975a938,
            0x0000000000000001, 0x3fc0933a2a7c2547, 0x402de6e978d4f7b1, 0x40a95e37cb033a8e,
            0x409119c291b0034b, 0xc044ca5540d62cd0, 0x4060795b1cd8e3e3, 0x0000000000000001,
            0x3fc093e5f0068445, 0x40366d4fdf3b6938, 0x40afeb8b269b247e, 0x40852f55af6b3684,
            0xc050b1cb1ed6de24, 0x405fd7ba8db4a23b, 0x0000000000000001, 0x3fc0950169060ffd,
            0x403de6e100590812, 0x40b26194c916e698, 0xc059000000000000, 0xc0574cbd1270a28b,
            0x406197bd83ac16f7, 0x0000000000000001, 0x3fc092c70abe88fa, 0xb91ae500d17bcb7b,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk4/icao-sea-level/rh100/no-zones",
        &[
            0x40b26d8cdfd99046, 0x4091b52726be48d7, 0x403def2d0a01b019, 0x40619e3c64ba4afd,
            0x405c28c6ec292f0f, 0x40755334ef1fc968, 0x0000000000000001, 0x40852232add96579,
            0x0000000000000001, 0x4089a9c701f26efb, 0x0000000000000001, 0x408e3c7c0f5339c4,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x00000000000074f0,
            0x3f50624dd2f1a9fc, 0x3fe6cce94d93d6d8, 0x3fda53cde9b99b1c, 0xbeb5370537af3450,
            0x4089b3cc79691703, 0x0000000000000001, 0x3fd313e8c75cf25d, 0x401def9db22d1201,
            0x40a079e552c5f742, 0x408fbd7e0cb69f82, 0xc03219aab6a7fa32, 0x4067d497e8cbf4cf,
            0x0000000000000001, 0x3fc0933aa4a8843d, 0x402def9db22d080e, 0x40a96d88c0f2c40a,
            0x4091239448eba390, 0xc044cd8927bad730, 0x4060806736cf54b0, 0x0000000000000001,
            0x3fc093e354255b89, 0x403673b645a1cfa7, 0x40afffb5f526a391, 0x40853c849ecde346,
            0xc050b492566fc05f, 0x405fe45701f8b1bb, 0x0000000000000001, 0x3fc094ff8d7428a0,
            0x403def2d0a01b019, 0x40b26d8cdfd99046, 0xc059000000000000, 0xc05750bc8e98a476,
            0x40619e3c64ba4afd, 0x0000000000000001, 0x3fc092c6d6589993, 0x63b258aec3d0a211,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk4/icao-sea-level/rh50/zones-rh20-rh80",
        &[
            0x40b24036d554964e, 0x409194c063a535e2, 0x403dd2051d472319, 0x4061833338cbdf6a,
            0x405bd29ccf6e2ded, 0x40754bf73b88949e, 0x0000000000000001, 0x40851e668d30261a,
            0x0000000000000001, 0x408997aa67dbe457, 0x0000000000000001, 0x408defa6a933877c,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x000000000000747e,
            0x3f50624dd2f1a9fc, 0x3fe6cceb5800893e, 0x3fda53d044f6d576, 0xbeb51f90a1956fc4,
            0x4089b3d11339b66d, 0x0000000000000001, 0x3fd310d89738e219, 0x401dd1eb851ebbf2,
            0x40a05d16bfc69304, 0x408f8b11dc3205f9, 0xc03210aa1ab58ab1, 0x4067ab8c89df6ba3,
            0x0000000000000001, 0x3fc0932ab8dcd740, 0x402dd26e978d49a9, 0x40a937a0d43fc354,
            0x4091026219e2cefa, 0xc044c3b6ad8a9c24, 0x406061315786756a, 0x0000000000000001,
            0x3fc09412d9f0b36b, 0x40365db22d0e5acd, 0x40afb5777bf7f2a7, 0x40850cecbf89713f,
            0xc050abe397e81688, 0x405fb036d7eb0406, 0x0000000000000001, 0x3fc0952a77d7040c,
            0x403dd2051d472319, 0x40b24036d554964e, 0xc059000000000000, 0xc05743dd43c4ae9f,
            0x4061833338cbdf6a, 0x0000000000000001, 0x3fc092cfe50e515f, 0x12630b1c03c4b35d,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk4/1500m-25C-830hPa/rh50/no-zones",
        &[
            0x40b5aae1ce48b709, 0x4094560e32f9f020, 0x403ff87e40b2286a, 0x40637da7bc4fb7fc,
            0x40613b6075ef235c, 0x4075b318f42cf8e0, 0x0000000000000001, 0x4089f44e29281390,
            0x0000000000000001, 0x408fb527f9a41236, 0x0000000000000001, 0x4092b615abc4bc4c,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000007ce4,
            0x3f50624dd2f1a9fc, 0x3fe6cd490d920a32, 0x3fda543c7982fc65, 0xbeb0ea7ff8865789,
            0x4089b4a45a3cdf74, 0x0000000000000001, 0x3fd33b6a535119a0, 0x401ff8d4fdf3bab0,
            0x40a2cdb196580799, 0x40921601d03b7cba, 0xc03276e358e46910, 0x406a0da3af5349a1,
            0x0000000000000001, 0x3fc096b8e4c0c2a0, 0x402ff8d4fdf3aec0, 0x40ad64cc8e0a751b,
            0x4093c40d432e20e0, 0xc0456a175618da9a, 0x40626fed0c07d327, 0x0000000000000001,
            0x3fc092acfdc89c3b, 0x4037fa9fbe76cfa8, 0x40b2ab5e08a768a3, 0x4088e44d9835872f,
            0xc05145fb3740cddf, 0x4061b9df319bb3e2, 0x0000000000000001, 0x3fc092dee95b4ddb,
            0x403ff87e40b2286a, 0x40b5aae1ce48b709, 0xc059000000000000, 0xc058302d1d137c95,
            0x40637da7bc4fb7fc, 0x0000000000000001, 0x3fc09316ce2b62f1, 0x93921dc44bd678b0,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk4/1500m-25C-830hPa/rh100/no-zones",
        &[
            0x40b5ca55913ef0fc, 0x40946f2ecf0ad052, 0x4040061d32e27994, 0x40638e01d4481e07,
            0x406158567a5ee571, 0x4075c3b860e9c586, 0x0000000000000001, 0x408a027daec5d4de,
            0x0000000000000001, 0x408fcfbb0e05e52f, 0x0000000000000001, 0x4092c57dba65c4dd,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000007d31,
            0x3f50624dd2f1a9fc, 0x3fe6cd4b425ff05a, 0x3fda543f05afb1a4, 0xbeb0d124490efad6,
            0x4089b4a953aeb44b, 0x0000000000000001, 0x3fd34233bb6b8b89, 0x40200624dd2f1cd2,
            0x40a2e3f9df83b11e, 0x40922ac0309311a9, 0xc0327a5660a0d185, 0x406a24260ea304da,
            0x0000000000000001, 0x3fc096c14dee874b, 0x40300624dd2f16e3, 0x40ad8b6f4d98aa75,
            0x4093dd41bdfe30fa, 0xc0456fab88f6a9a6, 0x406281db238f41c7, 0x0000000000000001,
            0x3fc092acfdc1bd94, 0x403809374bc6aef7, 0x40b2c554a022a46c, 0x408907d97d60f156,
            0xc0514b3f75a67b6c, 0x4061c9c2c5ac09fe, 0x0000000000000001, 0x3fc092ddab523e95,
            0x4040061d32e27994, 0x40b5ca55913ef0fc, 0xc059000000000000, 0xc058389d3e8aef11,
            0x40638e01d4481e07, 0x0000000000000001, 0x3fc09317ed916e51, 0x5d66f55d07c5e72c,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk4/1500m-25C-830hPa/rh50/zones-rh20-rh80",
        &[
            0x40b587ce01af7a70, 0x40943e47a86f5e7e, 0x403fe409f7d8465e, 0x4063698283bd8bc6,
            0x406117d3d3152a6b, 0x4075b318f42cf8e0, 0x0000000000000001, 0x4089ee63b256657a,
            0x0000000000000001, 0x408fa19cc78a4346, 0x0000000000000001, 0x40928c9e9b2c168f,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000007c94,
            0x3f50624dd2f1a9fc, 0x3fe6cd4b25c19e69, 0x3fda543ee4a40e30, 0xbeb0d26d3782806f,
            0x4089b4a913250c99, 0x0000000000000001, 0x3fd33b68b2159afd, 0x401fe45a1cac0c94,
            0x40a2b8d591115fbc, 0x40920418ab7b29a9, 0xc0327164804509fc, 0x4069ef9e7e04ae8f,
            0x0000000000000001, 0x3fc0965884d84dd7, 0x402fe45a1cac00b8, 0x40ad3c29a1e02044,
            0x4093ab9d002c6dc7, 0xc0456434f48b1faf, 0x406258c4552a5ca0, 0x0000000000000001,
            0x3fc092ade9b33e10, 0x4037eb4395810d04, 0x40b28efebf9fefec, 0x4088c0804da7bc56,
            0xc05140cfdde2165c, 0x4061a6707d847854, 0x0000000000000001, 0x3fc092ea0429ad4b,
            0x403fe409f7d8465e, 0x40b587ce01af7a70, 0xc059000000000000, 0xc05828522c7278fa,
            0x4063698283bd8bc6, 0x0000000000000001, 0x3fc093077c197144, 0x3dfb487f26947b9c,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk45/icao-sea-level/rh50/no-zones",
        &[
            0x40b26194ca6bd49c, 0x4091ab3df171ea4a, 0x403de6e101b1457d, 0x406197bd83f4d3ad,
            0x405c1407251f3b8a, 0x40754bf73b88949e, 0x0000000000000001, 0x40852097df5477af,
            0x0000000000000001, 0x4089a10c22d15c5f, 0x0000000000000001, 0x408e3ae1ead7f5bf,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000000bb2,
            0x3f50624dd2f1a9fc, 0x3fe6cce8183b03bf, 0x3fda53cc848655c8, 0xbeb544e8c17477bd,
            0x4089b3c9bff527f9, 0x0000000000000001, 0x3fd310db27e2d469, 0x401dceb78897e915,
            0x40a068361c0c6a5f, 0x408fa2d54e1e05b8, 0xc0320578b3559328, 0x4067d5d4486d0763,
            0x0000000000000001, 0x3fc093405b9d353d, 0x402de23d0bfa08ac, 0x40a95bdbce3e1575,
            0x40911ab2c4e8463b, 0xc044c681c454b8c7, 0x40607a530d7301c7, 0x0000000000000001,
            0x3fc093e42fd5aa0e, 0x40366bffcdab19e5, 0x40afea8d556ef191, 0x408532bbd37c9428,
            0xc050b0a8a274791d, 0x405fd74dd54e3591, 0x0000000000000001, 0x3fc09501da582e23,
            0x403de6e101b1457d, 0x40b26194ca6bd49c, 0xc059000000000000, 0xc0574cbd134e7897,
            0x406197bd83f4d3ad, 0x0000000000000001, 0x3fc092c70abe6beb, 0x6a5ba4ed5a76163c,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk45/icao-sea-level/rh100/no-zones",
        &[
            0x40b26d8cdc05c62c, 0x4091b52726479632, 0x403def2d048add54, 0x40619e3c61ede684,
            0x405c28c6e3372473, 0x40755334ef1fc968, 0x0000000000000001, 0x408524b2f4ff17b3,
            0x0000000000000001, 0x4089bf12eb36ced7, 0x0000000000000001, 0x408e42bea13075b0,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000000bb6,
            0x3f50624dd2f1a9fc, 0x3fe6cce94d93d649, 0x3fda53cde9b99a78, 0xbeb5370537b58e05,
            0x4089b3cc79691701, 0x0000000000000001, 0x3fd313e8c75cf25d, 0x401dd8f4f93bc01f,
            0x40a071c72fe97dfe, 0x408fb43616f8d1d5, 0xc03208d6d4f58398, 0x4067de4c6595fb3f,
            0x0000000000000001, 0x3fc093409790c211, 0x402dec7a7c9ddfb6, 0x40a96bf297813c93,
            0x40912435ad2b2811, 0xc044caf7d48f9e79, 0x4060810db1d7e3b3, 0x0000000000000001,
            0x3fc093e2288982b2, 0x403673ade225fb2e, 0x40afffaf9b9b2efd, 0x40853c9a5c10f47c,
            0xc050b48b16e1f770, 0x405fe4544bd99e98, 0x0000000000000001, 0x3fc094ff9046d6fa,
            0x403def2d048add54, 0x40b26d8cdc05c62c, 0xc059000000000000, 0xc05750bc8a81c9cc,
            0x40619e3c61ede684, 0x0000000000000001, 0x3fc092c6d659b64b, 0x15c1cb4785b1b225,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk45/icao-sea-level/rh50/zones-rh20-rh80",
        &[
            0x40b24031d93ed27b, 0x409194b92c7a28bb, 0x403dd1ffd6d2f6a3, 0x406183320c6e4c72,
            0x405bd29915076db5, 0x40754bf73b88949e, 0x0000000000000001, 0x40852698a239ba9b,
            0x0000000000000001, 0x4089a5b4c7fc4fe5, 0x0000000000000001, 0x408dfa333257d71a,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000000baa,
            0x3f50624dd2f1a9fc, 0x3fe6cceb580088b6, 0x3fda53d044f6d4da, 0xbeb51f90a19b8e79,
            0x4089b3d11339b66d, 0x0000000000000001, 0x3fd310d89738e219, 0x401dba3ca7503b01,
            0x40a054a457146200, 0x408f816592bd9c08, 0xc031ff0c9dd56270, 0x4067b5aed2674e96,
            0x0000000000000001, 0x3fc0932ee8f3df8d, 0x402dcdc22ab25a98, 0x40a93542af65d493,
            0x4091034b3877c582, 0xc044bfe3f754c206, 0x40606226713b1396, 0x0000000000000001,
            0x3fc09411078fa2b1, 0x40365ca3a4b55753, 0x40afb4a60c8a9db0, 0x40850f959134d495,
            0xc050aafaac4b778c, 0x405fafde2f17afa7, 0x0000000000000001, 0x3fc0952ad391ef01,
            0x403dd1ffd6d2f6a3, 0x40b24031d93ed27b, 0xc059000000000000, 0xc05743d9c7a7f689,
            0x406183320c6e4c72, 0x0000000000000001, 0x3fc092cfe59939e7, 0x5013d929178ee788,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk45/1500m-25C-830hPa/rh50/no-zones",
        &[
            0x40b5aae1cbde1e41, 0x4094560e326ea98f, 0x403ff87e3d8faf31, 0x40637da7bae21054,
            0x40613b6073689574, 0x4075b318f42cf8e0, 0x0000000000000001, 0x4089f6db67012373,
            0x0000000000000001, 0x408fb74c04d4a030, 0x0000000000000001, 0x4092b5f125c47eb0,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000000c81,
            0x3f50624dd2f1a9fc, 0x3fe6cd490d920a00, 0x3fda543c7982fc2d, 0xbeb0ea7ff8889b78,
            0x4089b4a45a3cdf74, 0x0000000000000001, 0x3fd33b6a535119a0, 0x401fe33269df971d,
            0x40a2c539e0690c12, 0x40921120644075dd, 0xc0326755f4b9495d, 0x406a1672b94fee0c,
            0x0000000000000001, 0x3fc096d62e7f0a2d, 0x402ff19934efcb2f, 0x40ad60b235216916,
            0x4093c583291faf11, 0xc0456454d003055d, 0x406271743ebd095f, 0x0000000000000001,
            0x3fc092acfc76eb35, 0x4037f8cc9a77e6f6, 0x40b2aa93d4861977, 0x4088e9585cd53473,
            0xc0514470b3c8ce25, 0x4061b999a8814af0, 0x0000000000000001, 0x3fc092df0f2098b6,
            0x403ff87e3d8faf31, 0x40b5aae1cbde1e41, 0xc059000000000000, 0xc058302d1acf4225,
            0x40637da7bae21054, 0x0000000000000001, 0x3fc09316ce2a41ff, 0xe2846ed26211c240,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk45/1500m-25C-830hPa/rh100/no-zones",
        &[
            0x40b5ca558b6fc7a5, 0x40946f2ecdd67053, 0x4040061d2f2a1fa1, 0x40638e01d0ef6642,
            0x40615856746f1358, 0x4075c3b860e9c586, 0x0000000000000001, 0x408a1f6c19d170b7,
            0x0000000000000001, 0x408fe0a2af7f6a37, 0x0000000000000001, 0x4092cc0daae84822,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000000c89,
            0x3f50624dd2f1a9fc, 0x3fe6cd4b425ff02b, 0x3fda543f05afb16e, 0xbeb0d124491114ea,
            0x4089b4a953aeb44b, 0x0000000000000001, 0x3fd34233bb6b8b89, 0x401ff7ad4b274531,
            0x40a2dbe1e884ed40, 0x409226159445eacc, 0xc0326b8a59c0e009, 0x406a2c88b6cf6e41,
            0x0000000000000001, 0x3fc096dd3751d88f, 0x4030030a0b1bbca2, 0x40ad87e622b50bd0,
            0x4093de82c70f5546, 0xc0456abb028839f6, 0x4062832b67ff8f65, 0x0000000000000001,
            0x3fc092ad0045b817, 0x40380828c36da988, 0x40b2c4df0381b43c, 0x40890ac6a6a3afc0,
            0xc0514a5b2ba76ba1, 0x4061c99aa9ddcfec, 0x0000000000000001, 0x3fc092ddc0c39c5e,
            0x4040061d2f2a1fa1, 0x40b5ca558b6fc7a5, 0xc059000000000000, 0xc058389d3930722b,
            0x40638e01d0ef6642, 0x0000000000000001, 0x3fc09317ed8ec7c5, 0xbd9de310ab5a25d3,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk45/1500m-25C-830hPa/rh50/zones-rh20-rh80",
        &[
            0x40b587cac388eed6, 0x40943e43058d7be2, 0x403fe406dfab2808, 0x40636981c83925c4,
            0x406117d288da14f9, 0x4075b318f42cf8e0, 0x0000000000000001, 0x4089f82c8a9afd4b,
            0x0000000000000001, 0x408fb60ffe8efcb5, 0x0000000000000001, 0x40929173837fef17,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000000c79,
            0x3f50624dd2f1a9fc, 0x3fe6cd4b25c19e38, 0x3fda543ee4a40df8, 0xbeb0d26d3784a954,
            0x4089b4a913250c98, 0x0000000000000001, 0x3fd33b68b2159afd, 0x401fce2f9c88755c,
            0x40a2b02f8182c972, 0x4091ff1b54491812, 0xc03261707bc73657, 0x4069f8a35904ed71,
            0x0000000000000001, 0x3fc09674e9ef6b77, 0x402fdcda5da06344, 0x40ad37eaa444720e,
            0x4093ad1c7a39a531, 0xc0455e3ba1c4380e, 0x40625a5768908484, 0x0000000000000001,
            0x3fc092adccd364a5, 0x4037e94e767e4779, 0x40b28e2507409ffe, 0x4088c5dab9acf6c8,
            0xc0513f28f223fcc4, 0x4061a6252febe6ae, 0x0000000000000001, 0x3fc092ea3162ec27,
            0x403fe406dfab2808, 0x40b587cac388eed6, 0xc059000000000000, 0xc05828503e715c19,
            0x40636981c83925c4, 0x0000000000000001, 0x3fc093077b9068c7, 0xfacf3e831b117c2e,
            0xffffffffffffffff,
        ],
    ),
];

/// `TrajectorySolver::solve` (RK4, flat) with Litz spin drift plus aerodynamic jump, and with
/// Magnus, at both stations, RH 0.
///
/// IMMUTABLE: no later step may change it, S8 included.
#[rustfmt::skip]
const TS_EFFECTS_DRY: Pinned = &[
    (
        "spin-drift+aerodynamic-jump/rk4/icao-sea-level/rh0",
        &[
            0x408f400000000000, 0x3fc29531b4c02da5, 0x3fff710188b4d261, 0x4074332862a9c83f,
            0x40828225fa61ccf0, 0x407544b4bbbd89be, 0x0000000000000001, 0x40882a073dd2e191,
            0x0000000000000001, 0x408d3ef52a1a69b4, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007af,
            0x3f50624dd2f1a9fc, 0x3fea53dbc7de79f6, 0x3f5c8f65ba5fe652, 0xbeadbec738d84c10,
            0x4089b3d0ff1ddf48, 0x0000000000000001, 0x3fd30dc641365e64, 0x3fdf6c8b4395810d,
            0x4075dc6a7669ebe7, 0xbfd5557103e41af9, 0xbfcda719c4e9f052, 0x40834a1a499f043b,
            0x0000000000000001, 0x3fd59320db3ce36a, 0x3fef74bc6a7ef9e2, 0x40834e66a4b6f473,
            0xc00574635b00c43e, 0xbfea0d9e61980b91, 0x407e08fd4faf4fb9, 0x0000000000000001,
            0x3fd8bf390195d7fb, 0x3ff79999999998b1, 0x4089e85de04a1ce6, 0xc01a9dc0bfc09660,
            0xbffa0eca931a63cc, 0x4077fc3accdae759, 0x0000000000000001, 0x3fdbac6e2ee21434,
            0x3fff710188b4d261, 0x408f400000000000, 0xc028086312555087, 0xc004ba5e3e4cf062,
            0x4074332862a9c83f, 0x0000000000000001, 0x3fcc72b45786cd9f, 0x9ac952a5b59216c1,
            0x0000000000000015, 0xe59fce7cd55f319a,
        ],
    ),
    (
        "magnus/rk4/icao-sea-level/rh0",
        &[
            0x408f400000000000, 0x3fc087c6d1b42c1f, 0x3fff7107101ff208, 0x40743327780611d9,
            0x408282244c660303, 0x407544b4bbbd89be, 0x0000000000000001, 0x40882a05d07809a1,
            0x0000000000000001, 0x408d3ef24ff99167, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007af,
            0x3f50624dd2f1a9fc, 0x3fea53dc3977f487, 0x3f5ae11060bd4959, 0xbeb552e2905e82da,
            0x4089b3d0ffc19150, 0x0000000000000001, 0x3fd30dc640fcfb58, 0x3fdf6c8b4395810d,
            0x4075dc6a836ce518, 0xbfd82087ee4d6b58, 0xbfd0e53d85fb1648, 0x40834a1a82c88555,
            0x0000000000000001, 0x3fd59320bd922206, 0x3fef74bc6a7ef9e2, 0x40834e660f7ba5f1,
            0xc006126799fd7dd6, 0xbfedbd00c394fb89, 0x407e08fcc99622dc, 0x0000000000000001,
            0x3fd8bf3940641fb5, 0x3ff79999999998b1, 0x4089e85c0701dec4, 0xc01b07e71bfde9c9,
            0xbffdee06c920eb2c, 0x4077fc3974602873, 0x0000000000000001, 0x3fdbac6eb64144eb,
            0x3fff7107101ff208, 0x408f400000000000, 0xc02848902d8fb2cd, 0xc008003533209bd1,
            0x40743327780611d9, 0x0000000000000001, 0x3fcc72a495bcaf36, 0x71987e71af93d83a,
            0x0000000000000015, 0xaa762a73a693ef92,
        ],
    ),
    (
        "spin-drift+aerodynamic-jump/rk4/1500m-25C-830hPa/rh0",
        &[
            0x408f400000000000, 0x3fc32443069d22fb, 0x3ffbfac9ab1b5ad9, 0x4078d36683a1f16a,
            0x408bf4ccf7fb1ac3, 0x4075a2608de6e31f, 0x0000000000000001, 0x408da7a198bca561,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d6,
            0x3f50624dd2f1a9fc, 0x3fea544a8688bb03, 0x3f5ccbf777e6924b, 0xbea2d255c8f3678e,
            0x4089b4a95bb1c741, 0x0000000000000001, 0x3fd334908e3db85d, 0x3fdbf7ced9168731,
            0x40743e14f1864420, 0xbfc71f3207396f50, 0xbfc24fd3083f0ccc, 0x4084de4990ba5f46,
            0x0000000000000001, 0x3fd4ffe0fff80b42, 0x3fec000000000006, 0x408271853b056bc0,
            0xc00052686a2289a0, 0xbfe0b2a632bf70b6, 0x40814d174568b8bf, 0x0000000000000001,
            0x3fd708664e9219d6, 0x3ff4fdf3b645a133, 0x4089615b6c351753, 0xc0153dfae6a4f8fb,
            0xbff12057fc0e0656, 0x407d1fce988385d5, 0x0000000000000001, 0x3fd9719b87a3153b,
            0x3ffbfac9ab1b5ad9, 0x408f400000000000, 0xc0239c1c6d04bc89, 0xbffbebf9a83bd07e,
            0x4078d36683a1f16a, 0x0000000000000001, 0x3fdb7f5aaf933507, 0x3448d64aeb122fbe,
            0x0000000000000015, 0xd9d82288fab0d9c7,
        ],
    ),
    (
        "magnus/rk4/1500m-25C-830hPa/rh0",
        &[
            0x408f400000000000, 0x3fc0bf443a04e734, 0x3ffbfacd6d3b1ca6, 0x4078d362f54d1ea0,
            0x408bf4c4f5bd48cd, 0x4075a2608de6e31f, 0x0000000000000001, 0x408da79f2fca5b85,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d6,
            0x3f50624dd2f1a9fc, 0x3fea544b088d8034, 0x3f5ae1819f25a179, 0xbeb103e8f01aad5f,
            0x4089b4a95c6c726f, 0x0000000000000001, 0x3fd334908dfd512b, 0x3fdbf7ced9168731,
            0x40743e151d21d3fb, 0xbfcd042e51494aa3, 0xbfc640ce5525ed9b, 0x4084de4a01facba7,
            0x0000000000000001, 0x3fd4ffe0d21efce2, 0x3fec000000000006, 0x40827184e7814f67,
            0xc000fe5dac9b69f1, 0xbfe4355c06e47a3e, 0x40814d176421850d, 0x0000000000000001,
            0x3fd7086637a4cc16, 0x3ff4fdf3b645a133, 0x4089615a14feb703, 0xc015b461cd3e6944,
            0xbff4cf6550dd380b, 0x407d1fcddb019791, 0x0000000000000001, 0x3fd9719be388f8c6,
            0x3ffbfacd6d3b1ca6, 0x408f400000000000, 0xc023e51707b817fd, 0xc00113c1ed56f059,
            0x4078d362f54d1ea0, 0x0000000000000001, 0x3fdb7f5c3beaf030, 0x903b459eb1743253,
            0x0000000000000015, 0xae15ca3192a5ecd4,
        ],
    ),
];

/// As `TS_EFFECTS_DRY`, at RH 50 and 100.
///
/// MAY CHANGE ONLY IN S8 (the unified aloft vapour model, plan C7): byte-identical through every
/// other step.
#[rustfmt::skip]
const TS_EFFECTS_HUMID: Pinned = &[
    (
        "spin-drift+aerodynamic-jump/rk4/icao-sea-level/rh50",
        &[
            0x408f400000000000, 0x3fc296032b994665, 0x3fff64ddc11a0b11, 0x40743e6391effff1,
            0x408296c08f2205b7, 0x40754bf73b88949e, 0x0000000000000001, 0x40882eb14cd345b7,
            0x0000000000000001, 0x408d487cbcbbd670, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ac,
            0x3f50624dd2f1a9fc, 0x3fea53dd2d38b88a, 0x3f5c8f673e64dddb, 0xbeada2fd27c7f794,
            0x4089b3d3b8df679d, 0x0000000000000001, 0x3fd310d7a9ee2237, 0x3fdf6c8b4395810d,
            0x4075de59f130add0, 0xbfd555193489dd97, 0xbfcd93d5d6f3fec2, 0x40834d8233a3cd7d,
            0x0000000000000001, 0x3fd594c719cdb993, 0x3fef6c8b4395810d, 0x40834d8a8dac0216,
            0xc0056858ebdd5828, 0xbfe9f2e28bdce82c, 0x407e15e51fab0ff3, 0x0000000000000001,
            0x3fd8bdfb5ac9b0aa, 0x3ff7916872b01fdd, 0x4089e792712a3b6b, 0xc01a8cb929f6ec1f,
            0xbff9f2757ee81982, 0x40780b5ab28ce727, 0x0000000000000001, 0x3fdba9b21ef5590e,
            0x3fff64ddc11a0b11, 0x408f400000000000, 0xc027f80225a61ff1, 0xc004a466475a724e,
            0x40743e6391effff1, 0x0000000000000001, 0x3fccbf065285c968, 0x8d14ad166e332b7a,
            0x0000000000000015, 0x8733f844b3720fd0,
        ],
    ),
    (
        "magnus/rk4/icao-sea-level/rh50",
        &[
            0x408f400000000000, 0x3fc0887808ed9384, 0x3fff64e343dbf014, 0x40743e62a553acf2,
            0x408296bedc987001, 0x40754bf73b88949e, 0x0000000000000001, 0x40882eafdef287fb,
            0x0000000000000001, 0x408d4879e065cb2e, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ac,
            0x3f50624dd2f1a9fc, 0x3fea53dd9ed23924, 0x3f5ae111cdf479fc, 0xbeb544fd87d61bea,
            0x4089b3d3b98319c2, 0x0000000000000001, 0x3fd310d7a9b4cfee, 0x3fdf6c8b4395810d,
            0x4075de59fe3e3b99, 0xbfd8206f54409a23, 0xbfd0db9b8ef5050c, 0x40834d826cd77608,
            0x0000000000000001, 0x3fd594c6fc151afe, 0x3fef6c8b4395810d, 0x40834d89f8ce073e,
            0xc0060655f69924a9, 0xbfeda0838524bef0, 0x407e15e498f02436, 0x0000000000000001,
            0x3fd8bdfb99c8a585, 0x3ff7916872b01fdd, 0x4089e79098bef283, 0xc01af6dbf8522691,
            0xbffdcf3c683f36ad, 0x40780b5957cddbf0, 0x0000000000000001, 0x3fdba9b2a83c8f79,
            0x3fff64e343dbf014, 0x408f400000000000, 0xc028382ec0b1c869, 0xc007e7ed9a9897fb,
            0x40743e62a553acf2, 0x0000000000000001, 0x3fccbef59baf5347, 0xadf77a9c05cdff00,
            0x0000000000000015, 0x4937c7dad155420a,
        ],
    ),
    (
        "spin-drift+aerodynamic-jump/rk4/icao-sea-level/rh100",
        &[
            0x408f400000000000, 0x3fc296d48c658ee5, 0x3fff58c24f3fd6ad, 0x407449a9acbbfdc8,
            0x4082ab7aae8a7b19, 0x40755334ef1fc968, 0x0000000000000001, 0x408836a26588e83c,
            0x0000000000000001, 0x408d520850cba8a2, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007a9,
            0x3f50624dd2f1a9fc, 0x3fea53de926cf2ef, 0x3f5c8f68c2408e7f, 0xbead87360b88c248,
            0x4089b3d67255280f, 0x0000000000000001, 0x3fd313e54a03357b, 0x3fdf5c28f5c28f63,
            0x4075d664edc8a5ca, 0xbfd52577e255a7e4, 0xbfcd63e97505972a, 0x40835393145dfb7b,
            0x0000000000000001, 0x3fd5950926a2f93a, 0x3fef5c28f5c28f63, 0x408348d04a50aa5b,
            0xc0054ef774525b7b, 0xbfe9cc62621ab311, 0x407e2683e26979a5, 0x0000000000000001,
            0x3fd8baffb5020bdc, 0x3ff7851eb851ea9f, 0x4089e3abfe79d4a9, 0xc01a71f223ac6782,
            0xbff9ceb1caaa2573, 0x40781d26f3f8ffbb, 0x0000000000000001, 0x3fdba5de50edd114,
            0x3fff58c24f3fd6ad, 0x408f400000000000, 0xc027e7b475d2d3af, 0xc0048e7ca8bd0bc2,
            0x407449a9acbbfdc8, 0x0000000000000001, 0x3fcd107243205730, 0xf5ab00c940256341,
            0x0000000000000015, 0x229176cc20e79ecf,
        ],
    ),
    (
        "magnus/rk4/icao-sea-level/rh100",
        &[
            0x408f400000000000, 0x3fc089292cfab686, 0x3fff58c7cd624b67, 0x407449a8bd6e087b,
            0x4082ab78f619a03d, 0x40755334ef1fc968, 0x0000000000000001, 0x408836a0f669d4ae,
            0x0000000000000001, 0x408d52057242a47a, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007a9,
            0x3f50624dd2f1a9fc, 0x3fea53df04067995, 0x3f5ae1133b04ccb9, 0xbeb53719f9b641ae,
            0x4089b3d672f8da51, 0x0000000000000001, 0x3fd313e549c9ed25, 0x3fdf5c28f5c28f63,
            0x4075d664fb25e6b8, 0xbfd7efc95aa96bb4, 0xbfd0c1ac65698d0d, 0x408353934dca6a1d,
            0x0000000000000001, 0x3fd5950908d4070c, 0x3fef5c28f5c28f63, 0x408348cfb650a2cd,
            0xc005eccd8db4a56d, 0xbfed7681aca484c6, 0x407e26835ba07c6e, 0x0000000000000001,
            0x3fd8bafff3e4d819, 0x3ff7851eb851ea9f, 0x4089e3aa27bd5495, 0xc01adc04977acf3b,
            0xbffda7ca14203d21, 0x40781d259743b378, 0x0000000000000001, 0x3fdba5dedc725365,
            0x3fff58c7cd624b67, 0x408f400000000000, 0xc02827e091cafeab, 0xc007cfb6ada54e95,
            0x407449a8bd6e087b, 0x0000000000000001, 0x3fcd10607c40dcf1, 0xde1c9e3bcac8e5d8,
            0x0000000000000015, 0x19ce4ebbe9c97afb,
        ],
    ),
    (
        "spin-drift+aerodynamic-jump/rk4/1500m-25C-830hPa/rh50",
        &[
            0x408f400000000000, 0x3fc325d650cca229, 0x3ffbe89cdf5b016b, 0x4078f39a0fecf529,
            0x408c3d8214e4a275, 0x4075b318f42cf8e0, 0x0000000000000001, 0x408db7b81953d778,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d2,
            0x3f50624dd2f1a9fc, 0x3fea544d126f8026, 0x3f5ccbfa418f18d3, 0xbea29fa3ed1968e2,
            0x4089b4ae549f6f53, 0x0000000000000001, 0x3fd33b66e1b80f22, 0x3fdbe76c8b439587,
            0x40743666426e3754, 0xbfc6ccb9a3e1a255, 0xbfc21d714e835f82, 0x4084e6e9975cc450,
            0x0000000000000001, 0x3fd502ea00556913, 0x3febef9db22d0e5c, 0x40826d9f9082e0bf,
            0xc0003bb81ca1a6bc, 0xbfe08843b6acf9a2, 0x408159f53ea163b9, 0x0000000000000001,
            0x3fd708c6637b00e5, 0x3ff4f1a9fbe76bf5, 0x40895f819109afed, 0xc015261733a02084,
            0xbff0f7bd7bdc4062, 0x407d3d96a5ab87cd, 0x0000000000000001, 0x3fd96e167ebe3a92,
            0x3ffbe89cdf5b016b, 0x408f400000000000, 0xc023860d35d16c13, 0xbffbaaad10803dd0,
            0x4078f39a0fecf529, 0x0000000000000001, 0x3fdb79abad8c69b6, 0xdd987e4d9bde35b1,
            0x0000000000000015, 0x5ae680d038593916,
        ],
    ),
    (
        "magnus/rk4/1500m-25C-830hPa/rh50",
        &[
            0x408f400000000000, 0x3fc0c08f81c847b1, 0x3ffbe8a0998f2888, 0x4078f3967f2f9e36,
            0x408c3d7a02cfe8e4, 0x4075b318f42cf8e0, 0x0000000000000001, 0x408db7b5adb076ee,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d2,
            0x3f50624dd2f1a9fc, 0x3fea544d947451f8, 0x3f5ae1843965d703, 0xbeb0ea90022d1491,
            0x4089b4ae555a1ac7, 0x0000000000000001, 0x3fd33b66e177a245, 0x3fdbe76c8b439587,
            0x407436666e4c7ddf, 0xbfccaf787e774bcb, 0xbfc60a33f662b746, 0x4084e6ea08f028fd,
            0x0000000000000001, 0x3fd502e9d248b316, 0x3febef9db22d0e5c, 0x40826d9f3dda403a,
            0xc000e788cafc344d, 0xbfe407381a215135, 0x408159f55d8e9f94, 0x0000000000000001,
            0x3fd708c64c751557, 0x3ff4f1a9fbe76bf5, 0x40895f803b78957b, 0xc0159c751f4dba4f,
            0xbff4a2d98f8f4920, 0x407d3d95e5e8d338, 0x0000000000000001, 0x3fd96e16dba60528,
            0x3ffbe8a0998f2888, 0x408f400000000000, 0xc023cf075d17bc11, 0xc000ef6845ff07e9,
            0x4078f3967f2f9e36, 0x0000000000000001, 0x3fdb79ad3e8c13c7, 0x07a7b57dee9bb5c3,
            0x0000000000000015, 0x28cd706839aae109,
        ],
    ),
    (
        "spin-drift+aerodynamic-jump/rk4/1500m-25C-830hPa/rh100",
        &[
            0x408f400000000000, 0x3fc32769fd06d77e, 0x3ffbd684d9b7bcac, 0x407913e8b56aa1b5,
            0x408c86d2ea1ef9f2, 0x4075c3b860e9c586, 0x0000000000000001, 0x408dc7e0bd09e494,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006cd,
            0x3f50624dd2f1a9fc, 0x3fea544f9e9e216f, 0x3f5ccbfd0b864b89, 0xbea26cec7aae39b6,
            0x4089b4b34e15196c, 0x0000000000000001, 0x3fd342304793419a, 0x3fdbd70a3d70a3dd,
            0x40742eb14f395a32, 0xbfc67a81d8715968, 0xbfc1eb4be5be723c, 0x4084ef885ddd81a6,
            0x0000000000000001, 0x3fd505e8cd6227ad, 0x3febd70a3d70a3dd, 0x4082653d21c24383,
            0xc00018ee3d045822, 0xbfe0554f3b83b9a8, 0x408168a6c931465a, 0x0000000000000001,
            0x3fd707bc4408cf14, 0x3ff4e147ae147a4d, 0x408959da17d10e77, 0xc0150515814c22ba,
            0xbff0c984eb257389, 0x407d5e4d988ba1be, 0x0000000000000001, 0x3fd96918ec005268,
            0x3ffbd684d9b7bcac, 0x408f400000000000, 0xc023701fed78d433, 0xbffb69a723924f68,
            0x407913e8b56aa1b5, 0x0000000000000001, 0x3fdb73e06a0b0eba, 0x306e394985830c2d,
            0x0000000000000015, 0x586865d7016e9e2c,
        ],
    ),
    (
        "magnus/rk4/1500m-25C-830hPa/rh100",
        &[
            0x408f400000000000, 0x3fc0c1db156ac65a, 0x3ffbd6888b9c4e08, 0x407913e5227c6129,
            0x408c86cac89afd4a, 0x4075c3b860e9c586, 0x0000000000000001, 0x408dc7de4ec5f506,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006cd,
            0x3f50624dd2f1a9fc, 0x3fea545020a2ffe7, 0x3f5ae186d3ef683a, 0xbeb0d13448f6d8f8,
            0x4089b4b34ecfc527, 0x0000000000000001, 0x3fd342304752aabc, 0x3fdbd70a3d70a3dd,
            0x40742eb17b5b7039, 0xbfcc5b01722dba2c, 0xbfc5d3d7f6ac045d, 0x4084ef88cfc844de,
            0x0000000000000001, 0x3fd505e89f1a6839, 0x3febd70a3d70a3dd, 0x4082653cd0634ed4,
            0xc000c470710c1ae0, 0xbfe3cea4e268b9fb, 0x408168a6e8a8f966, 0x0000000000000001,
            0x3fd707bc2cb6eb6f, 0x3ff4e147ae147a4d, 0x408959d8c4b78f3f, 0xc0157b58a700ec45,
            0xbff46f624dacf344, 0x407d5e4cd753207a, 0x0000000000000001, 0x3fd9691949a21807,
            0x3ffbd6888b9c4e08, 0x408f400000000000, 0xc023b919a17593cc, 0xc000cb382a1952fa,
            0x407913e5227c6129, 0x0000000000000001, 0x3fdb73e1fcb4a411, 0x0f19ce03648f5afe,
            0x0000000000000015, 0xb19f6c728a8b96fb,
        ],
    ),
];

/// `calculate_zero_angle_with_station_presence` to 300 m (every zero trial flies
/// `calculate_acceleration`), Euler, RK4 and RK45 at both stations, RH 0. The bisection stops at
/// 1e-7 rad, so a zero pin moves only when a trial height changes sign; the flight groups are the
/// sensitive ones.
///
/// IMMUTABLE: no later step may change it, S8 included.
#[rustfmt::skip]
const ZERO_DRY: Pinned = &[
    (
        "zero-300m/euler/icao-sea-level/rh0",
        &[
            0x3f666cccccccccce,
        ],
    ),
    (
        "zero-300m/euler/1500m-25C-830hPa/rh0",
        &[
            0x3f65b1999999999a,
        ],
    ),
    (
        "zero-300m/rk4/icao-sea-level/rh0",
        &[
            0x3f665b3333333334,
        ],
    ),
    (
        "zero-300m/rk4/1500m-25C-830hPa/rh0",
        &[
            0x3f65a0cccccccccd,
        ],
    ),
    (
        "zero-300m/rk45/icao-sea-level/rh0",
        &[
            0x3f665b3333333334,
        ],
    ),
    (
        "zero-300m/rk45/1500m-25C-830hPa/rh0",
        &[
            0x3f65a1999999999a,
        ],
    ),
];

/// As `ZERO_DRY`, at RH 50 and 100.
///
/// MAY CHANGE ONLY IN S8 (the unified aloft vapour model, plan C7): byte-identical through every
/// other step.
#[rustfmt::skip]
const ZERO_HUMID: Pinned = &[
    (
        "zero-300m/euler/icao-sea-level/rh50",
        &[
            0x3f666b3333333334,
        ],
    ),
    (
        "zero-300m/euler/icao-sea-level/rh100",
        &[
            0x3f66680000000000,
        ],
    ),
    (
        "zero-300m/euler/1500m-25C-830hPa/rh50",
        &[
            0x3f65acccccccccce,
        ],
    ),
    (
        "zero-300m/euler/1500m-25C-830hPa/rh100",
        &[
            0x3f65a9999999999a,
        ],
    ),
    (
        "zero-300m/rk4/icao-sea-level/rh50",
        &[
            0x3f6658cccccccccd,
        ],
    ),
    (
        "zero-300m/rk4/icao-sea-level/rh100",
        &[
            0x3f66566666666667,
        ],
    ),
    (
        "zero-300m/rk4/1500m-25C-830hPa/rh50",
        &[
            0x3f659ccccccccccd,
        ],
    ),
    (
        "zero-300m/rk4/1500m-25C-830hPa/rh100",
        &[
            0x3f65980000000000,
        ],
    ),
    (
        "zero-300m/rk45/icao-sea-level/rh50",
        &[
            0x3f6658cccccccccd,
        ],
    ),
    (
        "zero-300m/rk45/icao-sea-level/rh100",
        &[
            0x3f66566666666667,
        ],
    ),
    (
        "zero-300m/rk45/1500m-25C-830hPa/rh50",
        &[
            0x3f659ccccccccccd,
        ],
    ),
    (
        "zero-300m/rk45/1500m-25C-830hPa/rh100",
        &[
            0x3f6598cccccccccd,
        ],
    ),
];

/// Seeded `TrajectorySolver` Monte Carlo at the `BallisticInputs` default humidity (50 %): the
/// CLI's seeded entry at ICAO sea level, and the presence-aware entry the FFI Monte Carlo runs on,
/// at 1500 m.
///
/// MAY CHANGE ONLY IN S8 (the unified aloft vapour model, plan C7): byte-identical through every
/// other step.
#[rustfmt::skip]
const MONTE_CARLO_DEFAULT_HUMIDITY: Pinned = &[
    (
        "seeded/icao-sea-level/default-rh",
        &[
            0x0000000000000028, 0x408283ec179233e4, 0x407de78c5bd1c144, 0x3fe9a388297250f1,
            0x3fbd9a387ecfc050, 0x3fb999999999999a, 0x1cc4cc9120623eec, 0x26d90fbd1be6c87a,
            0x84d2ca5a6bb6d013,
        ],
    ),
    (
        "seeded/1500m-25C-830hPa/default-rh",
        &[
            0x0000000000000028, 0x40811edd92f414bb, 0x40819cd6c796e5f5, 0x3fc6881e5105c4ac,
            0xbfe46500230625c3, 0x3fb999999999999a, 0x3d5340a556cc4237, 0x9f0551ec01e0b64c,
            0x8e64d9aa9c5de2eb,
        ],
    ),
];

/// Station density and sound speed: `calculate_atmosphere` over altitude, T/P override and RH, plus
/// `calculate_air_density_cimp` and `moist_speed_of_sound`. Station air is not part of the aloft
/// model, so this holds through S8 at every RH.
///
/// IMMUTABLE: no later step may change it, S8 included.
#[rustfmt::skip]
const STATION_AIR: Pinned = &[
    (
        "calculate_atmosphere/0m/None/None",
        &[
            0x3ff39bbc4502831a, 0x407544b4bbbd89be, 0x3ff38fc7f384bdc2, 0x40754a270eebb447,
            0x3ff38bcd293e575b, 0x40754bf73b88949e, 0x3ff37be889cadc0f, 0x40755334ef1fc968,
        ],
    ),
    (
        "calculate_atmosphere/0m/Some(-10.0)/Some(1013.25)",
        &[
            0x3ff57a5c3875328d, 0x40745327d13f55d3, 0x3ff578294d9eb953, 0x40745407d21d0c74,
            0x3ff5776db99ec930, 0x4074545278ea08e8, 0x3ff5747fb79a126c, 0x4074557cff5a700d,
        ],
    ),
    (
        "calculate_atmosphere/0m/Some(15.0)/Some(1013.25)",
        &[
            0x3ff39bbc4502831a, 0x407544b4bbbd89be, 0x3ff38fc7f384bdc2, 0x40754a270eebb447,
            0x3ff38bcd293e575b, 0x40754bf73b88949e, 0x3ff37be889cadc0f, 0x40755334ef1fc968,
        ],
    ),
    (
        "calculate_atmosphere/0m/Some(35.0)/Some(1013.25)",
        &[
            0x3ff255141827cfe6, 0x4075fe7d7dfedebf, 0x3ff23030bf3f1d97, 0x4076110c17cdd4b6,
            0x3ff223edc141c06e, 0x40761734e3b44084, 0x3ff1f30ccbb9f7f0, 0x40762fb6ae909fdf,
        ],
    ),
    (
        "calculate_atmosphere/0m/Some(25.0)/Some(830.0)",
        &[
            0x3fef0abe88bf329d, 0x4075a2608de6e31f, 0x3feedfcb5c869ad7, 0x4075aeed3440992e,
            0x3feed1817ad31cde, 0x4075b318f42cf8e0, 0x3fee987d4e8b1180, 0x4075c3b860e9c586,
        ],
    ),
    (
        "calculate_atmosphere/0m/Some(35.0)/Some(850.0)",
        &[
            0x3feec19ef85a5434, 0x4075fe7d7dfedebf, 0x3fee77d514c3cab2, 0x4076149623784ade,
            0x3fee5f4de33de950, 0x40761bea233a66c6, 0x3fedfd860ea25751, 0x4076390ae07a80f6,
        ],
    ),
    (
        "calculate_atmosphere/0m/Some(30.0)/None",
        &[
            0x3ff2a2ad8bd46b61, 0x4075d09fa46719c8, 0x3ff28663394efd16, 0x4075de84bfbb5ab1,
            0x3ff27cfad6bd3993, 0x4075e322a44ca630, 0x3ff257757348568f, 0x4075f5874a8c8459,
        ],
    ),
    (
        "calculate_atmosphere/0m/None/Some(900.0)",
        &[
            0x3ff16a7b662c3c1d, 0x407544b4bbbd89be, 0x3ff15e869133bbb5, 0x40754ad5c392d181,
            0x3ff15a8b99f8f012, 0x40754ce00376620b, 0x3ff14aa63fb9b08a, 0x407555053851103a,
        ],
    ),
    (
        "calculate_atmosphere/1500m/None/None",
        &[
            0x3ff0eff2f5ec5434, 0x4074e7d205c76b24, 0x3ff0e9817076b595, 0x4074eb28c958b4cc,
            0x3ff0e75c08c70bc0, 0x4074ec457a661d16, 0x3ff0dec895243a4a, 0x4074f0b7191433f8,
        ],
    ),
    (
        "calculate_atmosphere/1500m/Some(-10.0)/Some(1013.25)",
        &[
            0x3ff57a5c3875328d, 0x40745327d13f55d3, 0x3ff578294d9eb953, 0x40745407d21d0c74,
            0x3ff5776db99ec930, 0x4074545278ea08e8, 0x3ff5747fb79a126c, 0x4074557cff5a700d,
        ],
    ),
    (
        "calculate_atmosphere/1500m/Some(15.0)/Some(1013.25)",
        &[
            0x3ff39bbc4502831a, 0x407544b4bbbd89be, 0x3ff38fc7f384bdc2, 0x40754a270eebb447,
            0x3ff38bcd293e575b, 0x40754bf73b88949e, 0x3ff37be889cadc0f, 0x40755334ef1fc968,
        ],
    ),
    (
        "calculate_atmosphere/1500m/Some(35.0)/Some(1013.25)",
        &[
            0x3ff255141827cfe6, 0x4075fe7d7dfedebf, 0x3ff23030bf3f1d97, 0x4076110c17cdd4b6,
            0x3ff223edc141c06e, 0x40761734e3b44084, 0x3ff1f30ccbb9f7f0, 0x40762fb6ae909fdf,
        ],
    ),
    (
        "calculate_atmosphere/1500m/Some(25.0)/Some(830.0)",
        &[
            0x3fef0abe88bf329d, 0x4075a2608de6e31f, 0x3feedfcb5c869ad7, 0x4075aeed3440992e,
            0x3feed1817ad31cde, 0x4075b318f42cf8e0, 0x3fee987d4e8b1180, 0x4075c3b860e9c586,
        ],
    ),
    (
        "calculate_atmosphere/1500m/Some(35.0)/Some(850.0)",
        &[
            0x3feec19ef85a5434, 0x4075fe7d7dfedebf, 0x3fee77d514c3cab2, 0x4076149623784ade,
            0x3fee5f4de33de950, 0x40761bea233a66c6, 0x3fedfd860ea25751, 0x4076390ae07a80f6,
        ],
    ),
    (
        "calculate_atmosphere/1500m/Some(30.0)/None",
        &[
            0x3fef1a498232812d, 0x4075d09fa46719c8, 0x3feee1b219c0575b, 0x4075e141bd3752da,
            0x3feecee05390637a, 0x4075e6c7a5a34716, 0x3fee83d0f06c2510, 0x4075fcc43438e639,
        ],
    ),
    (
        "calculate_atmosphere/1500m/None/Some(900.0)",
        &[
            0x3ff2070a10e82ba0, 0x4074e7d205c76b24, 0x3ff20098b3165965, 0x4074eaf548a3c6b0,
            0x3ff1fe7358c59f82, 0x4074ec00d5a97108, 0x3ff1f5e01b8a4767, 0x4074f02e067f8272,
        ],
    ),
    (
        "calculate_atmosphere/3000m/None/None",
        &[
            0x3fed1c125d30c882, 0x40748956c61b0f3a, 0x3fed157a3b22af8e, 0x40748b4a84573010,
            0x3fed1347c698902c, 0x40748bf1049f34ca, 0x3fed0a7f3b0dc53f, 0x40748e8a9f87b8c1,
        ],
    ),
    (
        "calculate_atmosphere/3000m/Some(-10.0)/Some(1013.25)",
        &[
            0x3ff57a5c3875328d, 0x40745327d13f55d3, 0x3ff578294d9eb953, 0x40745407d21d0c74,
            0x3ff5776db99ec930, 0x4074545278ea08e8, 0x3ff5747fb79a126c, 0x4074557cff5a700d,
        ],
    ),
    (
        "calculate_atmosphere/3000m/Some(15.0)/Some(1013.25)",
        &[
            0x3ff39bbc4502831a, 0x407544b4bbbd89be, 0x3ff38fc7f384bdc2, 0x40754a270eebb447,
            0x3ff38bcd293e575b, 0x40754bf73b88949e, 0x3ff37be889cadc0f, 0x40755334ef1fc968,
        ],
    ),
    (
        "calculate_atmosphere/3000m/Some(35.0)/Some(1013.25)",
        &[
            0x3ff255141827cfe6, 0x4075fe7d7dfedebf, 0x3ff23030bf3f1d97, 0x4076110c17cdd4b6,
            0x3ff223edc141c06e, 0x40761734e3b44084, 0x3ff1f30ccbb9f7f0, 0x40762fb6ae909fdf,
        ],
    ),
    (
        "calculate_atmosphere/3000m/Some(25.0)/Some(830.0)",
        &[
            0x3fef0abe88bf329d, 0x4075a2608de6e31f, 0x3feedfcb5c869ad7, 0x4075aeed3440992e,
            0x3feed1817ad31cde, 0x4075b318f42cf8e0, 0x3fee987d4e8b1180, 0x4075c3b860e9c586,
        ],
    ),
    (
        "calculate_atmosphere/3000m/Some(35.0)/Some(850.0)",
        &[
            0x3feec19ef85a5434, 0x4075fe7d7dfedebf, 0x3fee77d514c3cab2, 0x4076149623784ade,
            0x3fee5f4de33de950, 0x40761bea233a66c6, 0x3fedfd860ea25751, 0x4076390ae07a80f6,
        ],
    ),
    (
        "calculate_atmosphere/3000m/Some(30.0)/None",
        &[
            0x3fe9ca7b897e158e, 0x4075d09fa46719c8, 0x3fe991e1bfcb16cb, 0x4075e4a927cf34c5,
            0x3fe97f0f19278dfc, 0x4075eb4f132fa734, 0x3fe933fb905ffd9f, 0x407605bf88cac107,
        ],
    ),
    (
        "calculate_atmosphere/3000m/None/Some(900.0)",
        &[
            0x3ff2aef16e3a2dcb, 0x40748956c61b0f3a, 0x3ff2aba5b2756a38, 0x40748adc68b32d38,
            0x3ff2aa8c94c30344, 0x40748b5e3d1fea94, 0x3ff2a628c22235fb, 0x40748d6550abfd61,
        ],
    ),
    (
        "primitives/-10C/1013.25hPa",
        &[
            0x3ff57a5c3875328d, 0x40745327d13f55d3, 0x3ff578294d9eb953, 0x40745407d21d0c74,
            0x3ff5776db99ec930, 0x4074545278ea08e8, 0x3ff5747fb79a126c, 0x4074557cff5a700d,
        ],
    ),
    (
        "primitives/15C/1013.25hPa",
        &[
            0x3ff39bbc4502831a, 0x407544b4bbbd89be, 0x3ff38fc7f384bdc2, 0x40754a270eebb447,
            0x3ff38bcd293e575b, 0x40754bf73b88949e, 0x3ff37be889cadc0f, 0x40755334ef1fc968,
        ],
    ),
    (
        "primitives/35C/850hPa",
        &[
            0x3feec19ef85a5434, 0x4075fe7d7dfedebf, 0x3fee77d514c3cab2, 0x4076149623784ade,
            0x3fee5f4de33de950, 0x40761bea233a66c6, 0x3fedfd860ea25751, 0x4076390ae07a80f6,
        ],
    ),
];

/// `get_local_atmosphere_humid` from 200 m below each station base to 25 km above it, RH 0.
///
/// IMMUTABLE: no later step may change it, S8 included.
#[rustfmt::skip]
const LOCAL_AIR_DRY: Pinned = &[
    (
        "get_local_atmosphere_humid/0m/15C/1013.25hPa/rh0",
        &[
            0x3ff001be42c894d1, 0x3ff3fcd49940ab38, 0x407550f97c7d0dde, 0x3ff39bbc4502831a,
            0x407544b4bbbd89be, 0x3ff396eb09fbdf67, 0x4075441783de95ea, 0x3ff2af2fb3680ee8,
            0x407525ea95679bf5, 0x3fed1bc8a1e25138, 0x40748956c61b0f3a, 0x3fd3f8f95cba1384,
            0x4072711d212dfca7, 0x3fa4881ddf6be4e0, 0x4072a639fe0c27ee,
        ],
    ),
    (
        "get_local_atmosphere_humid/0m/35C/1013.25hPa/rh0",
        &[
            0x3fedee2b418f5ddf, 0x3ff2a9ec477cbd1f, 0x40760a5ad58891ee, 0x3ff255141827cfe6,
            0x4075fe7d7dfedebf, 0x3ff250ddec1b5b6a, 0x4075fde5764a36fd, 0x3ff18607936558ca,
            0x4075e0b8d42c9e64, 0x3febc38c5057fb06, 0x40754983bdf6a58b, 0x3fd4b40bd3b64cf0,
            0x40734637ee5a76ed, 0x3fa94bb28b2eb94f, 0x4073790fadafc534,
        ],
    ),
    (
        "get_local_atmosphere_humid/1500m/25C/830hPa/rh0",
        &[
            0x3fe957235ff4e566, 0x3fef9f315c288890, 0x4075ae6ed93deaf1, 0x3fef0abe88bf329e,
            0x4075a2608de6e31f, 0x3fef03608df8e1b1, 0x4075a1c6117d757d, 0x3feda0c4a0631172,
            0x407584201ecddd6a, 0x3fe7489f5cbf0513, 0x4074ea609b63c1f3, 0x3fd02d73554d0683,
            0x40734396676250bd, 0x3fa3b8b984c9e59e, 0x407385cb8d148e12,
        ],
    ),
];

/// As `LOCAL_AIR_DRY`, at RH 37.5, 50 and 100. The legacy wrapper shares `air_at`'s composition, so
/// the aloft model reaches it.
///
/// MAY CHANGE ONLY IN S8 (the unified aloft vapour model, plan C7): byte-identical through every
/// other step.
#[rustfmt::skip]
const LOCAL_AIR_HUMID: Pinned = &[
    (
        "get_local_atmosphere_humid/0m/15C/1013.25hPa/rh37.5",
        &[
            0x3fefeff8156b0190, 0x3ff3f0a515eb40f8, 0x407556c4e83fe7c3, 0x3ff38fc7f384bdc2,
            0x40754a270eebb447, 0x3ff38af9a8573184, 0x407549858048fff7, 0x3ff2a3cb98ab4bb8,
            0x40752a90ecf6fa10, 0x3fed0a099919cabc, 0x40748b4a84573010, 0x3fd3eccc3344e3ae,
            0x40727128438c0e5e, 0x3fa47b99715eda28, 0x4072a6d2d5880b92,
        ],
    ),
    (
        "get_local_atmosphere_humid/0m/15C/1013.25hPa/rh50",
        &[
            0x3fefe978c0b95f99, 0x3ff3ec9696e126cb, 0x407558b2b3d661b2, 0x3ff38bcd293e575b,
            0x40754bf73b88949e, 0x3ff386ffd85bb999, 0x40754b543b94cc68, 0x3ff2a000d0c51b9b,
            0x40752c1d4418e6f0, 0x3fed042134f33a87, 0x40748bf1049f34ca, 0x3fd3e8be7c989b5f,
            0x4072712bf9a8842b, 0x3fa4776ead709a35, 0x4072a705c5eb047f,
        ],
    ),
    (
        "get_local_atmosphere_humid/0m/15C/1013.25hPa/rh100",
        &[
            0x3fefcf861a6fc551, 0x3ff3dc6344be4da5, 0x40756066804ea148, 0x3ff37be889cadc0e,
            0x40755334ef1fc968, 0x3ff3771f2078cb2c, 0x4075528c2e93d0ba, 0x3ff290dbebf414ae,
            0x4075324c6e4f3b92, 0x3fecec8958dff5c3, 0x40748e8a9f87b8c1, 0x3fd3d88e4aa3e131,
            0x4072713ad20be048, 0x3fa466ca762dca3c, 0x4072a7d17ced8500,
        ],
    ),
    (
        "get_local_atmosphere_humid/0m/35C/1013.25hPa/rh37.5",
        &[
            0x3fedb1f18bd4c298, 0x3ff2845e358632dc, 0x40761de47ea09667, 0x3ff23030bf3f1d97,
            0x4076110c17cdd4b6, 0x3ff22c030cae4489, 0x40761067c5639cf7, 0x3ff162c4d9713dab,
            0x4075f0fe67895834, 0x3feb8bae9ab79ead, 0x40755171d9d085ff, 0x3fd48a6325d105a8,
            0x40734694b1e295e5, 0x3fa918cc3c3e6277, 0x40737cd8ca4f2c26,
        ],
    ),
    (
        "get_local_atmosphere_humid/0m/35C/1013.25hPa/rh50",
        &[
            0x3fed9deca93c544a, 0x3ff277e277dda207, 0x407624604351cb19, 0x3ff223edc141c06d,
            0x40761734e3b44084, 0x3ff21fc2dfd52074, 0x4076168c811d023d, 0x3ff1570c57d34684,
            0x4075f665b8abccde, 0x3feb791cae37a824, 0x40755415481d3eb1, 0x3fd47c8a2e895f7d,
            0x407346b39cfa4779, 0x3fa907e0edab71cd, 0x40737e1b82b44d35,
        ],
    ),
    (
        "get_local_atmosphere_humid/0m/35C/1013.25hPa/rh100",
        &[
            0x3fed4e1f5710466d, 0x3ff2461f4ab6c110, 0x40763e2a6a8fd863, 0x3ff1f30ccbb9f7f0,
            0x40762fb6ae909fdf, 0x3ff1eeed24fcbc97, 0x40762efe365895f3, 0x3ff128536e07382f,
            0x40760be925a65a29, 0x3feb2f161f3b62a2, 0x40755e9cab57291e, 0x3fd44556e3da1f97,
            0x4073472f45978236, 0x3fa8c46f0c00d980, 0x40738324cf5fcdd8,
        ],
    ),
    (
        "get_local_atmosphere_humid/1500m/25C/830hPa/rh37.5",
        &[
            0x3fe93413b93edc70, 0x3fef7370caf40144, 0x4075bbb5d570479a, 0x3feedfcb5c869ad7,
            0x4075aeed3440992e, 0x3feed8779334c2f6, 0x4075ae499e288ad2, 0x3fed77c649627711,
            0x40758efe1fa1c767, 0x3fe7286833b3c171, 0x4074ef5127d348c9, 0x3fd01711313f5a9e,
            0x4073440a34ee82b2, 0x3fa39d7011195105, 0x40738b3d4ba22864,
        ],
    ),
    (
        "get_local_atmosphere_humid/1500m/25C/830hPa/rh50",
        &[
            0x3fe92869b2a1e351, 0x3fef64e29492ef49, 0x4075c01f5437a8c7, 0x3feed1817ad31cdd,
            0x4075b318f42cf8e0, 0x3feeca31159ec7e4, 0x4075b27259f86c6a, 0x3fed6a2305bed675,
            0x4075929b18532010, 0x3fe71db087ee4228, 0x4074f0f62b589c66, 0x3fd00f9ee29a8fdf,
            0x40734430cd9bfef8, 0x3fa3945c244f526c, 0x40738d0d39190a7a,
        ],
    ),
    (
        "get_local_atmosphere_humid/1500m/25C/830hPa/rh100",
        &[
            0x3fe8f9de745ca092, 0x3fef2acdbea5e1b4, 0x4075d1b3ebdce4fc, 0x3fee987d4e8b1180,
            0x4075c3b860e9c586, 0x3fee913a7168d867, 0x4075c305cc4fe920, 0x3fed33b7b575542d,
            0x4075a1033a4c34b7, 0x3fe6f2ec5d613d4f, 0x4074f787b7e22b68, 0x3fcfe3d0296dc83e,
            0x407344cb2a76bd5b, 0x3fa37022e7218b82, 0x40739449ac7142e7,
        ],
    ),
];

/// NWS pressure altitude, the FAA-rule and density-matched density altitudes, the density-altitude
/// entry inverse, and the `atmosphere.density_altitude` service.
///
/// IMMUTABLE: no later step may change it, S8 included.
#[rustfmt::skip]
const DENSITY_ALTITUDE: Pinned = &[
    (
        "nws+faa_rule",
        &[
            0x40c34b31491ee685, 0x40c206b86858f15a, 0x40c7e2b86858f15a, 0x40cc92b86858f15a,
            0x40b526b9369782b4, 0x40aceedaf8d55920, 0x40ba2f6d7c6aac90, 0x40c1c7b6be355648,
            0x40a9538c7d7d8c7f, 0x408faa889a7c0dc0, 0x40af5aa2269f0370, 0x40b90d51134f81b8,
            0x0000000000000000, 0xc0a7700000000000, 0x0000000000000000, 0x40a2c00000000000,
            0xc08ee69fff9eb3cf, 0xc0b0802b51dc7694, 0xc09320ad4771da4f, 0x40925f52b88e25b2,
        ],
    ),
    (
        "density_matched",
        &[
            0x0000000000000000, 0x4038d4354b1ba398, 0x40408d64a6af37c4, 0x40508d18421d1bfb,
            0x4080775c583b7191, 0x408260c130bbfe1c, 0x408303fe7d390079, 0x408591860b8a29b5,
            0x40a2829345b4bba8, 0x40a2ed5fa618853c, 0x40a31100757b4010, 0x40a39fa9dc2ee9dc,
        ],
    ),
    (
        "resolve_atmosphere_for_density_altitude",
        &[
            0xc07f400000000000, 0x403240e558eab901, 0x4090cb342f4c10b1, 0xc08a785954c82b03,
            0x403e000000000000, 0x40917d1a7ead66d4, 0x0000000000000000, 0x402e000000000000,
            0x408faa0000000000, 0xc07bb2a9905605cc, 0x403e000000000000, 0x4090ae9fa2909cb6,
            0x4097700000000000, 0x4014f53fd4ff53fe, 0x408a6bf72e919069, 0x408803b8ddac7574,
            0x403e000000000000, 0x408ce1eafaa5cbd2, 0x40a7700000000000, 0xc012158056015804,
            0x4085e7cb77736cee, 0x409ef06341c1f6e7, 0x403e000000000000, 0x4088e6b33ebef6bd,
        ],
    ),
    (
        "service/icao-sea-level/rh0.5",
        &[
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x40408d64a6af37c4, 0x405b271ab20539e2, 0x3ff38bcd293e575b,
        ],
    ),
    (
        "service/1500m/rh0",
        &[
            0x4099c9a1813a27c0, 0x40b526b9369782b4, 0x40a2d1edf499808b, 0x40bedf6d7c6aac90,
            0x40a2829345b4bba8, 0x40be5d40dfde913b, 0x3fef0abe88bf329d,
        ],
    ),
    (
        "service/1500m/rh0.9",
        &[
            0x4099c9a1813a27c0, 0x40b526b9369782b4, 0x40a2d1edf499808b, 0x40bedf6d7c6aac90,
            0x40a3831c99adb95a, 0x40c0010a5b1bedb2, 0x3feea3e00d26d766,
        ],
    ),
];

/// Muzzle Miller Sg (the published dry T/P correction, kept deliberately: plan C8 was rejected),
/// Litz spin drift, the fast path's aerodynamic-jump offset, and `TrajectorySolver`'s
/// aerodynamic-jump components, which depend on station T/P only and so are pinned at every RH.
///
/// IMMUTABLE: no later step may change it, S8 included.
#[rustfmt::skip]
const STABILITY_SPIN_JUMP: Pinned = &[
    (
        "miller-sg+litz",
        &[
            0x40045740e6bb5acd, 0x40045740e6bb5acd, 0x3fa11cb106e70547, 0x3fbe6b7dbdc6f437,
            0x3fdb09dac2f5964e, 0x4005c0ae1b64e63c, 0x4005c0ae1b64e63c, 0x3fa1eb41ee5424aa,
            0x3fbfdab359ae9a75, 0x3fdc503ef5884f5b, 0x4009b189fc2b839c, 0x4009b189fc2b839c,
            0x3fa42bd9b24b539e, 0x3fc1edd9c1d7d891, 0x3fdfdf4f6cab5d10, 0x40043714bfd9ca79,
        ],
    ),
    (
        "fast-path-aerodynamic-jump",
        &[
            0x3f205867a5a4382e, 0x3f22a0c071f3acbf, 0x3f2154bea136b61a,
        ],
    ),
    (
        "trajectory-solver-jump/icao-sea-level/rh0",
        &[
            0x3fdb6fe38aae2234, 0x0000000000000000, 0x3f205867a5a4382e, 0x0000000000000000,
            0x0000000000000000, 0x3ff0000000000000,
        ],
    ),
    (
        "trajectory-solver-jump/icao-sea-level/rh50",
        &[
            0x3fdb6fe38aae2234, 0x0000000000000000, 0x3f205867a5a4382e, 0x0000000000000000,
            0x0000000000000000, 0x3ff0000000000000,
        ],
    ),
    (
        "trajectory-solver-jump/icao-sea-level/rh100",
        &[
            0x3fdb6fe38aae2234, 0x0000000000000000, 0x3f205867a5a4382e, 0x0000000000000000,
            0x0000000000000000, 0x3ff0000000000000,
        ],
    ),
    (
        "trajectory-solver-jump/1500m-25C-830hPa/rh0",
        &[
            0x3fdf44c41c7ffa25, 0x0000000000000000, 0x3f22a0c071f3acbf, 0x0000000000000000,
            0x0000000000000000, 0x3ff0000000000000,
        ],
    ),
    (
        "trajectory-solver-jump/1500m-25C-830hPa/rh50",
        &[
            0x3fdf44c41c7ffa25, 0x0000000000000000, 0x3f22a0c071f3acbf, 0x0000000000000000,
            0x0000000000000000, 0x3ff0000000000000,
        ],
    ),
    (
        "trajectory-solver-jump/1500m-25C-830hPa/rh100",
        &[
            0x3fdf44c41c7ffa25, 0x0000000000000000, 0x3f22a0c071f3acbf, 0x0000000000000000,
            0x0000000000000000, 0x3ff0000000000000,
        ],
    ),
];

/// Kernel direct (uniform) mode, density 1.1 kg/m3 and sound speed 335 m/s: `compute_derivatives`,
/// `try_integrate_trajectory` RK4 and RK45, `try_solve_trajectory_rust`, `fast_integrate` and
/// `fast_integrate_with_segments`. S3 migrates the tuple to the uniform `KernelAtmosphere`; the
/// values do not move.
///
/// IMMUTABLE: no later step may change it, S8 included.
#[rustfmt::skip]
const KERNEL_DIRECT: Pinned = &[
    (
        "direct(1.1,335)/compute_derivatives",
        &[
            0x4085e00000000000, 0x4000000000000000, 0x3fe0000000000000, 0xc076ac6e3f7438fc,
            0xc025afb2130696ad, 0xc002a8379c1390e7, 0x4085e00000000000, 0x4000000000000000,
            0x3fe0000000000000, 0xc076ac6e3f7438fc, 0xc025afb2130696ad, 0xc002a8379c1390e7,
            0x4085e00000000000, 0x4000000000000000, 0x3fe0000000000000, 0xc076ac6e3f7438fc,
            0xc025afb2130696ad, 0xc002a8379c1390e7,
        ],
    ),
    (
        "direct(1.1,335)/try_integrate_trajectory/RK4",
        &[
            0x000000000000053a, 0x3f50624dd2f1a9fc, 0x3fea54181aca9560, 0x3f5ae14daff80ea8,
            0xbeb2fed8eca8f5f0, 0x4089b44295d3f77e, 0x3ffa2a09e555ab17, 0xbf628c34051cb06f,
            0x3fc5604189374bcb, 0x406067aeaaee9976, 0x3fc09eae2efadc71, 0xbf9edc34d4689311,
            0x4087739e526510eb, 0xbfb06961ee44a839, 0xbfd690e6182acdb6, 0x3fd5604189374bcb,
            0x406f65aa5e236415, 0xbf8b1b7fbeb5fc06, 0xbfbd7e66a9f794ac, 0x40857af117b0fe51,
            0xbffa08daae0792ef, 0xbfe5182fd95fa7ce, 0x3fe56872b020c4a0, 0x407cf52390bd4e4c,
            0xbff09720d107704e, 0xbfdb24f535fe4ec5, 0x40823cd14916c2ea, 0xc011a9d801be921f,
            0xbff29dd37d5b09f4, 0x3ff00c49ba5e3541, 0x408420dc2238e699, 0xc007671c1cb3701e,
            0xbfec34ea03686fcb, 0x407f5ec19431d7af, 0xc01b5d0a39f7bc40, 0xbff8f7e2a782ba52,
            0x3ff56170bb1bfe32, 0x4089000000000000, 0xc01647a18a633f5c, 0xbff74f9b29ede477,
            0x407b3c74d8066844, 0xc021faa63476d1f3, 0xbffe1ca79d293b0a, 0x22b5035cbded18b2,
        ],
    ),
    (
        "direct(1.1,335)/try_integrate_trajectory/RK45",
        &[
            0x000000000000002a, 0x3f9db22d0e56041a, 0x4037ac2ea7be9f12, 0x3fa62411703499fc,
            0xbf4ef4767df78f84, 0x40894d312c3aee86, 0x3ff5653d7f108c6e, 0xbfb09c132f35a428,
            0x3fc1cac083126e97, 0x405b8453920f005c, 0x3fc05b57646898f2, 0xbf958bc38f1abbe8,
            0x4087cf4b8a3cd3e7, 0x3fcb157db51b948c, 0xbfd3005d0f7f6409, 0x3fd27ef9db22d0e7,
            0x406b7c4bd6d6527d, 0x3faa070517279370, 0xbfb65965c7b5d659, 0x4085fc2dd89948ed,
            0xbff380b9624b0d66, 0xbfe294fb8addaa95, 0x3fe3ced916872b05, 0x407b1cbcd87224e2,
            0xbfea6d58f1cdcaf9, 0xbfd78743437dc959, 0x4082abdca4216d8e, 0xc010195521b557e9,
            0xbff1897f1820e25a, 0x3fee0c49ba5e3545, 0x40831c3c265a6284, 0xc00403483de37f27,
            0xbfe914df5a52a091, 0x408021eadffa48b3, 0xc019990793234b72, 0xbff7dadbc1f74bf0,
            0x3ff56174f2369494, 0x4089000000000000, 0xc01647baab32e716, 0xbff74fac065851db,
            0x407b3c7638edfc56, 0xc021faa71d5845f0, 0xbffe1ca5e610f280, 0x79eef4abc7699d33,
        ],
    ),
    (
        "direct(1.1,335)/try_solve_trajectory_rust/RK45",
        &[
            0x000000000000002a, 0x3f9db22d0e56041a, 0x4037ac2ea7be9f12, 0x3fa62411703499fc,
            0xbf4ef4767df78f84, 0x40894d312c3aee86, 0x3ff5653d7f108c6e, 0xbfb09c132f35a428,
            0x3fc1cac083126e97, 0x405b8453920f005c, 0x3fc05b57646898f2, 0xbf958bc38f1abbe8,
            0x4087cf4b8a3cd3e7, 0x3fcb157db51b948c, 0xbfd3005d0f7f6409, 0x3fd27ef9db22d0e7,
            0x406b7c4bd6d6527d, 0x3faa070517279370, 0xbfb65965c7b5d659, 0x4085fc2dd89948ed,
            0xbff380b9624b0d66, 0xbfe294fb8addaa95, 0x3fe3ced916872b05, 0x407b1cbcd87224e2,
            0xbfea6d58f1cdcaf9, 0xbfd78743437dc959, 0x4082abdca4216d8e, 0xc010195521b557e9,
            0xbff1897f1820e25a, 0x3fee0c49ba5e3545, 0x40831c3c265a6284, 0xc00403483de37f27,
            0xbfe914df5a52a091, 0x408021eadffa48b3, 0xc019990793234b72, 0xbff7dadbc1f74bf0,
            0x3ff56174f2369494, 0x4089000000000000, 0xc01647baab32e716, 0xbff74fac065851db,
            0x407b3c7638edfc56, 0xc021faa71d5845f0, 0xbffe1ca5e610f280, 0x79eef4abc7699d33,
        ],
    ),
    (
        "direct(1.1,335)/fast_integrate",
        &[
            0x0000000000000001, 0x000000000000053a, 0x3f50624dd2f1a9fc, 0x3fea54181aca9564,
            0x3f5ae14daff80eab, 0xbeb2fed8d32a14cd, 0x4089b44295d3f782, 0x3ffa2a09e555ab1a,
            0xbf628c33ec37af30, 0x3fc5604189374bc7, 0x406067aeaaee9a46, 0x3fc09eae2efadd8d,
            0xbf9edc34aafcf8b3, 0x4087739e5265133f, 0xbfb06961ee44958a, 0xbfd690e5f9e11733,
            0x3fd5604189374bc7, 0x406f65aa5e236714, 0xbf8b1b7fbeb5dc53, 0xbfbd7e6682617cb6,
            0x40857af117b10289, 0xbffa08daae079323, 0xbfe5182fbd0f90d6, 0x3fe56872b020c49c,
            0x407cf52390bd538d, 0xbff09720d107714b, 0xbfdb24f5118f7904, 0x40823cd14916c985,
            0xc011a9d801be9403, 0xbff29dd3645e6441, 0x3ff00c49ba5e353f, 0x408420dc2238eb9c,
            0xc007671c1cb37297, 0xbfec34e9dd8c97dd, 0x407f5ec19431e67f, 0xc01b5d0a39f7c03a,
            0xbff8f7e285ff9cdf, 0x3ff56170bb1bf5f6, 0x4089000000000000, 0xc01647a18a632e5f,
            0xbff74f9b0aa42b4d, 0x407b3c74d8067de9, 0xc021faa63476ce43, 0xbffe1ca774bebe9a,
            0x2f74a56952e3bc31, 0x0000000000000001, 0x3ff56170bb1bf5f6, 0x0000000000000001,
            0x3fc47ae147ae147b, 0x0000000000000000,
        ],
    ),
    (
        "direct(1.1,335)/fast_integrate_with_segments",
        &[
            0x0000000000000001, 0x000000000000002c, 0x3f947ae147ae147b, 0x40305e0670ff89ce,
            0x3f9f84f64a414d9d, 0xbf3d85692ffcd02f, 0x40896e0d5a46cc76, 0x3ff6eb74d24c6fcd,
            0xbfa6ff79e1320dc2, 0x3fbeb851eb851eb7, 0x4057e1084047d2f4, 0x3fbf396f8546be4b,
            0xbf9024bda2ca489f, 0x40880ebee39d4342, 0x3fd9ab7b395e543e, 0xbfd088c93ac97ef6,
            0x3fd147ae147ae149, 0x4069ce6e7d104d5a, 0x3fb284713698f383, 0xbfb39b77f956111e,
            0x4086343173736edd, 0xbff0b36d15396809, 0xbfe17e346ae36007, 0x3fe3333333333336,
            0x407a664878fa8886, 0xbfe8069f25a07b92, 0xbfd6362429003682, 0x4082d71c1d833219,
            0xc00efdec3f21c2d3, 0xbff11de026ae2956, 0x3fee66666666666c, 0x4083498dc83277b2,
            0xc004944fec2947d4, 0xbfe99bc59b01f5fb, 0x40800de6bc0fa09b, 0xc019e776bb02e7d7,
            0xbff80caafc7eadf1, 0x3ff56175b79e5436, 0x4089000000000000, 0xc01647bf4544f048,
            0xbff74faf1bf75493, 0x407b3c7679667d3d, 0xc021faa747f27ecd, 0xbffe1ca595d98829,
            0xadfbd6fa0a7cce53, 0x0000000000000001, 0x3ff56175b79e5436, 0x0000000000000001,
            0x3fc48b78fc905657, 0x0000000000000000,
        ],
    ),
];

/// Kernel standard mode at RH 0, station and zones dry, the same families plus the lean Monte
/// Carlo, at both stations. Compared to a relative 1e-15, not bits: the plan allows last-ulp
/// movement here when S3 moves the kernels onto `air_at`.
///
/// IMMUTABLE: no later step may change it, S8 included.
///
/// Floats compared to a relative 1e-15; counts exactly.
#[rustfmt::skip]
const KERNEL_STANDARD_RH0: Pinned = &[
    (
        "standard/icao-sea-level/rh0/no-zones/compute_derivatives",
        &[
            0x4085e00000000000, 0x4000000000000000, 0x3fe0000000000000, 0xc0796c4b1919ea8e,
            0xc025f00c7bc547d3, 0xc004eb654ac7cb3f, 0x4085e00000000000, 0x4000000000000000,
            0x3fe0000000000000, 0xc0796cf2a2d1aec5, 0xc025f01bcd1dc317, 0xc004ebef26e420a2,
            0x4085e00000000000, 0x4000000000000000, 0x3fe0000000000000, 0xc078aa3fc70dd055,
            0xc025de4ebd697a17, 0xc0044bb9998d8fa2,
        ],
    ),
    (
        "standard/icao-sea-level/rh0/no-zones/try_integrate_trajectory/RK4",
        &[
            0x000000000000057a, 0x3f50624dd2f1a9fc, 0x3fea53dc3977ed11, 0x3f5ae1107de50452,
            0xbeb552e2ad478984, 0x4089b3cdab92e112, 0x3ffa299288053224, 0xbf64d2149f8527f2,
            0x3fc666666666666b, 0x40610e3ca59022cd, 0x3fc061ffe1f5d1b9, 0xbfa2dcb062de0734,
            0x408715b4ccd304b6, 0xbfc3379a50ff431f, 0xbfda37afc06bb48e, 0x3fd666666666666b,
            0x407034b26c3dedb2, 0xbfa5978562e1c1b8, 0xbfc1e3b439e53acf, 0x4084d99864b5a375,
            0xbffc4a33226c1eae, 0xbfe83b323375b99c, 0x3fe66e978d4fdf40, 0x407d871d64b1ecc4,
            0xbff2ccc1e8b29683, 0xbfe03fccaf414fd8, 0x4081462cfbdab2e2, 0xc0126dedd88b917b,
            0xbff50395bbd47315, 0x3ff0d0e560418921, 0x40844f93a2279a94, 0xc009c6ab668cdeb4,
            0xbff0b8a2da880892, 0x407d11cad15def76, 0xc01c212189a99279, 0xbffbd4b14a05b9bd,
            0x3ff6687ae452635f, 0x4089000000000000, 0xc0183c6842b024ac, 0xbffb6bc3d44ba868,
            0x4078c487226e33c3, 0xc02253058c797605, 0xc0009775bd9e8292,
        ],
    ),
    (
        "standard/icao-sea-level/rh0/no-zones/try_integrate_trajectory/RK45",
        &[
            0x000000000000002b, 0x3f9db22d0e56041a, 0x4037a624fdb3b837, 0x3fa61e3dc784bc1d,
            0xbf515b05b3d0f96a, 0x4089404c9235b779, 0x3ff5592e57cc6b88, 0xbfb29d69ef4f001e,
            0x3fc1cac083126e97, 0x405b63e7edabd079, 0x3fc0434aa0167c6c, 0xbf98112e50782fae,
            0x4087976aded25135, 0x3fca0b9219067bc8, 0xbfd52c8f91d9217c, 0x3fd27ef9db22d0e7,
            0x406b3c2d3619edbf, 0x3fa925ca720f4b0c, 0xbfb8d7a27e607947, 0x40859699f80a49dc,
            0xbff38d1d8e8ad712, 0xbfe48e86dd7230e2, 0x3fe3ced916872b05, 0x407aa2561baee546,
            0xbfea522d40ed679e, 0xbfd9e871df61e61d, 0x408205b1b142ef72, 0xc00fe21837a53e60,
            0xbff326ffa743e4f0, 0x3ff024dd2f1a9fc1, 0x4083b1b8f4b2d2b9, 0xc0077581604e9ccc,
            0xbfef22a504769fe6, 0x407da8a6b53f1005, 0xc01b091496e7933b, 0xbffb18fd394b965d,
            0x3ff6687dc8c6c7af, 0x4089000000000000, 0xc0183c78676e6178, 0xbffb6bcf661d462f,
            0x4078c48847802e76, 0xc022530665865545, 0xc0009775074baf58,
        ],
    ),
    (
        "standard/icao-sea-level/rh0/no-zones/try_solve_trajectory_rust/RK45",
        &[
            0x000000000000002b, 0x3f9db22d0e56041a, 0x4037a624fdb3b837, 0x3fa61e3dc784bc1d,
            0xbf515b05b3d0f96a, 0x4089404c9235b779, 0x3ff5592e57cc6b88, 0xbfb29d69ef4f001e,
            0x3fc1cac083126e97, 0x405b63e7edabd079, 0x3fc0434aa0167c6c, 0xbf98112e50782fae,
            0x4087976aded25135, 0x3fca0b9219067bc8, 0xbfd52c8f91d9217c, 0x3fd27ef9db22d0e7,
            0x406b3c2d3619edbf, 0x3fa925ca720f4b0c, 0xbfb8d7a27e607947, 0x40859699f80a49dc,
            0xbff38d1d8e8ad712, 0xbfe48e86dd7230e2, 0x3fe3ced916872b05, 0x407aa2561baee546,
            0xbfea522d40ed679e, 0xbfd9e871df61e61d, 0x408205b1b142ef72, 0xc00fe21837a53e60,
            0xbff326ffa743e4f0, 0x3ff024dd2f1a9fc1, 0x4083b1b8f4b2d2b9, 0xc0077581604e9ccc,
            0xbfef22a504769fe6, 0x407da8a6b53f1005, 0xc01b091496e7933b, 0xbffb18fd394b965d,
            0x3ff6687dc8c6c7af, 0x4089000000000000, 0xc0183c78676e6178, 0xbffb6bcf661d462f,
            0x4078c48847802e76, 0xc022530665865545, 0xc0009775074baf58,
        ],
    ),
    (
        "standard/icao-sea-level/rh0/no-zones/fast_integrate",
        &[
            0x0000000000000001, 0x000000000000057a, 0x3f50624dd2f1a9fc, 0x3fea53dc3977ed14,
            0x3f5ae1107de50454, 0xbeb552e290a8a792, 0x4089b3cdab92e117, 0x3ffa299288053228,
            0xbf64d21483932793, 0x3fc6666666666667, 0x40610e3ca59023cc, 0x3fc061ffe1f5d30c,
            0xbfa2dcb0498cfe7a, 0x408715b4ccd30777, 0xbfc3379a50ff38b8, 0xbfda37af9d3b5384,
            0x3fd6666666666667, 0x407034b26c3def8a, 0xbfa5978562e1b902, 0xbfc1e3b421e26242,
            0x4084d99864b5a84f, 0xbffc4a33226c1f36, 0xbfe83b3212efd7b0, 0x3fe66e978d4fdf3c,
            0x407d871d64b1f2f2, 0xbff2ccc1e8b297b8, 0xbfe03fcc99721069, 0x4081462cfbdab9ee,
            0xc0126dedd88b938a, 0xbff503959fa0068c, 0x3ff0d0e560418937, 0x40844f93a227a043,
            0xc009c6ab668ce17a, 0xbff0b8a2c4169b41, 0x407d11cad15dfed6, 0xc01c212189a996bf,
            0xbffbd4b124ab0d5d, 0x3ff6687ae452593e, 0x4089000000000000, 0xc0183c6842b00f0b,
            0xbffb6bc3af7db857, 0x4078c487226e4b40, 0xc02253058c7971c2, 0xc0009775a7599c16,
            0x0000000000000001, 0x3ff6687ae452593e, 0x0000000000000001, 0x3fc47ae147ae147b,
            0x0000000000000000,
        ],
    ),
    (
        "standard/icao-sea-level/rh0/no-zones/fast_integrate_with_segments",
        &[
            0x0000000000000001, 0x000000000000002c, 0x3f947ae147ae147b, 0x40305b2312b14a3b,
            0x3f9f7f482e0513c2, 0xbf408ead903d4827, 0x408965158cab9959, 0x3ff6e2d2342e52a9,
            0xbfa9c99a7cc57dea, 0x3fc0a3d70a3d70a3, 0x4059afe0929af724, 0x3fbff4296517063b,
            0xbf951c52e2d5dd8b, 0x4087b889cefc5d6d, 0x3fd2c214c18bbc37, 0xbfd3e2e181c1b057,
            0x3fd1eb851eb851ed, 0x406a74b38c8697a5, 0x3fae96ddaeccbed3, 0xbfb761e284f31f9c,
            0x4085b365735ad9af, 0xbff23cbd4c9b4622, 0xbfe3ff378f6bd0a9, 0x3fe3851eb851eb88,
            0x407a4f1820c2afdc, 0xbfe92ee8712f57bf, 0xbfd938e911eb0326, 0x40821b742cc1421b,
            0xc00f536850676a4b, 0xbff2f0d9b829ab4b, 0x3ff0000000000003, 0x40838f7b4c0ad1ae,
            0xc006f97891bfffdf, 0xbfeea6258eceae28, 0x407dc989fcb4acc7, 0xc01acc7c4e25b9af,
            0xbffaf011b20e92ec, 0x3ff6687bf1fbafc4, 0x4089000000000000, 0xc0183c6e23ba9b6a,
            0xbffb6bc80af0dcaa, 0x4078c4878d12003c, 0xc0225305dbbf9262, 0xc00097757b46e1d2,
            0x0000000000000001, 0x3ff6687bf1fbafc4, 0x0000000000000001, 0x3fc46e8825b30b45,
            0x0000000000000000,
        ],
    ),
    (
        "standard/icao-sea-level/rh0/zones-rh0-rh0/compute_derivatives",
        &[
            0x4085e00000000000, 0x4000000000000000, 0x3fe0000000000000, 0xc07935246c5fc8d6,
            0xc025eb01a21aa9ae, 0xc004be03a3c83bef, 0x4085e00000000000, 0x4000000000000000,
            0x3fe0000000000000, 0xc079b0b56bd8b005, 0xc025f64dcb041233, 0xc00523b113fce89d,
            0x4085e00000000000, 0x4000000000000000, 0x3fe0000000000000, 0xc078e9e0000fc7c4,
            0xc025e41ff2e5ffe4, 0xc00480147aee43d5,
        ],
    ),
    (
        "standard/icao-sea-level/rh0/zones-rh0-rh0/try_integrate_trajectory/RK4",
        &[
            0x0000000000000578, 0x3f50624dd2f1a9fc, 0x3fea53e0bd8822aa, 0x3f5ae1151b5a9df1,
            0xbeb525ef65ec6524, 0x4089b3d67cd5eb5f, 0x3ffa299b8892b41f, 0xbf64a63266f2592f,
            0x3fc666666666666b, 0x40611022d66bf589, 0x3fc06489674e79d9, 0xbfa2b6e17241c946,
            0x40871ada74b994d6, 0xbfc323d8eeda746b, 0xbfda04740158cfd4, 0x3fd666666666666b,
            0x4070382169798170, 0xbfa58667650051a9, 0xbfc1c186e6b60b7d, 0x4084e26f9e429e01,
            0xbffc4aefc0f8fe17, 0xbfe80f324d8f4799, 0x3fe666666666666b, 0x407d88ec6c4a32c6,
            0xbff2baf09032d051, 0xbfe01a87d81b6ae3, 0x4081507da5141924, 0xc01265b090ce080b,
            0xbff4e9ea8c893530, 0x3ff0ccccccccccb7, 0x4084528ce02cd7ad, 0xc009b52ba6886d36,
            0xbff0a0da3a109e12, 0x407d161bb1b318ef, 0xc01c1099fbbdcf28, 0xbffbcf52915a4e07,
            0x3ff660d5cfc75575, 0x4089000000000000, 0xc0182423fb4bca91, 0xbffb4d2f81f665e8,
            0x4078c0bbc2279f69, 0xc022436e2d2d30c2, 0xc00099d20c7fb61b,
        ],
    ),
    (
        "standard/icao-sea-level/rh0/zones-rh0-rh0/try_integrate_trajectory/RK45",
        &[
            0x000000000000002b, 0x3f9db22d0e56041a, 0x4037a699b5c4afdc, 0x3fa61eae695cd45a,
            0xbf5136b746ac5915, 0x408941460187d8b5, 0x3ff55a179db2c330, 0xbfb2769e96a308e4,
            0x3fc1cac083126e97, 0x405b665aba5e84c6, 0x3fc0451b4000ab8a, 0xbf97e0702d9ea3bb,
            0x40879ba2d74392f0, 0x3fca1f9c005ce502, 0xbfd50291b5227f5b, 0x3fd27ef9db22d0e7,
            0x406b4106148aa9e1, 0x3fa936b0537e83af, 0xbfb8a76312ad869b, 0x40859e4b2f2b1137,
            0xbff38c3b56c1f202, 0xbfe4683e3ad13ff1, 0x3fe3ced916872b05, 0x407aab802b86f0fc,
            0xbfea54042c6774d7, 0xbfd9bad5b35eea79, 0x4082103745842a48, 0xc00fe4eeb28bcb87,
            0xbff30cd0c8b32d61, 0x3ff024dd2f1a9fc1, 0x4083b86879c11a40, 0xc007733343709a49,
            0xbfef015ee9e04531, 0x407daaefc7c47ad0, 0xc01b005406a608e2, 0xbffb16254236887f,
            0x3ff660cb3732aa49, 0x4089000000000000, 0xc0182414bc23271e, 0xbffb4d051fa37d35,
            0x4078c0d39fcc44f7, 0xc022436886cf2d38, 0xc00099c3339858c7,
        ],
    ),
    (
        "standard/icao-sea-level/rh0/zones-rh0-rh0/fast_integrate",
        &[
            0x0000000000000001, 0x0000000000000578, 0x3f50624dd2f1a9fc, 0x3fea53e0bd8822ac,
            0x3f5ae1151b5a9df3, 0xbeb525ef4989d81b, 0x4089b3d67cd5eb63, 0x3ffa299b8892b423,
            0xbf64a6324b3b3f3b, 0x3fc6666666666667, 0x40611022d66bf688, 0x3fc06489674e7b24,
            0xbfa2b6e159237f80, 0x40871ada74b9978c, 0xbfc323d8eeda6a3b, 0xbfda0473de6d328d,
            0x3fd6666666666667, 0x407038216979833e, 0xbfa5866765004939, 0xbfc1c186cee11229,
            0x4084e26f9e42a2c8, 0xbffc4aefc0f8fe9f, 0xbfe80f322d4473e1, 0x3fe6666666666667,
            0x407d88ec6c4a38c6, 0xbff2baf09032d184, 0xbfe01a87c27e3109, 0x4081507da5142014,
            0xc01265b090ce0a17, 0xbff4e9ea70773c62, 0x3ff0cccccccccccc, 0x4084528ce02cdd3f,
            0xc009b52ba6887008, 0xbff0a0da23bf1ccd, 0x407d161bb1b32803, 0xc01c1099fbbdd373,
            0xbffbcf526c06d6b0, 0x3ff660d5cfc74b89, 0x4089000000000000, 0xc0182423fb4bb589,
            0xbffb4d2f5d51814e, 0x4078c0bbc227b67d, 0xc022436e2d2d2caf, 0xc00099d1f637a489,
            0x0000000000000001, 0x3ff660d5cfc74b89, 0x0000000000000001, 0x3fc47ae147ae147b,
            0x0000000000000000,
        ],
    ),
    (
        "standard/icao-sea-level/rh0/zones-rh0-rh0/fast_integrate_with_segments",
        &[
            0x0000000000000001, 0x000000000000002c, 0x3f947ae147ae147b, 0x40305b5ae082575c,
            0x3f9f7fb5ece794d8, 0xbf406bf622a61af5, 0x408965c2ecaea0cc, 0x3ff6e379229de08a,
            0xbfa993ac7f0dc334, 0x3fc0a3d70a3d70a3, 0x4059b207d2fbc864, 0x3fbff76ed1eb7a47,
            0xbf94f174bd33548c, 0x4087bc843d6f48f5, 0x3fd2cbfa6359af90, 0xbfd3bb4831e709f9,
            0x3fd1eb851eb851ed, 0x406a79465cae5906, 0x3faea77b0b23de4c, 0xbfb7345c6c66b9ad,
            0x4085bae81d87570b, 0xbff23ba862ec17f4, 0xbfe3d9d69a91afe1, 0x3fe3851eb851eb88,
            0x407a5811d6de51e6, 0xbfe930a5bc35d4f0, 0xbfd90c3d88f56c41, 0x4082264e17b876c3,
            0xc00f56bd6953f05d, 0xbff2d5d8f9eedf43, 0x3ff0000000000003, 0x4083962c27887ead,
            0xc006f758a59fc228, 0xbfee84d8cbf241d9, 0x407dcc381ea28f76, 0xc01ac40232bf0c6f,
            0xbffaecbbfcf7596e, 0x3ff660c400c39c27, 0x4089000000000000, 0xc018240fc0b96b2d,
            0xbffb4ce845e71dee, 0x4078c0e764538686, 0xc022436462b509a7, 0xc00099b6e75ca78d,
            0x0000000000000001, 0x3ff660c400c39c27, 0x0000000000000001, 0x3fc470907ed70a3d,
            0x0000000000000000,
        ],
    ),
    (
        "standard/icao-sea-level/rh0/solve_trajectory_for_monte_carlo",
        &[
            0x40186f9b75e3423e, 0xbffb6bc3af7db857, 0x3ff6687ae452593e, 0x4078c64f21272d6d,
            0x408bd758c54b58d3, 0x3ff2a342190f650a, 0x0000000000000000, 0x4089000000000000,
        ],
    ),
    (
        "standard/1500m-25C-830hPa/rh0/no-zones/compute_derivatives",
        &[
            0x4085e00000000000, 0x4000000000000000, 0x3fe0000000000000, 0xc07442c83e1e18ac,
            0xc0257739922c1010, 0xc000abfb1464d566, 0x4085e00000000000, 0x4000000000000000,
            0x3fe0000000000000, 0xc074434916b27050, 0xc025774559e7b49a, 0xc000ac6519fd9e3f,
            0x4085e00000000000, 0x4000000000000000, 0x3fe0000000000000, 0xc073ad7dce3fa010,
            0xc0256993505b0393, 0xc0003122c40b64fd,
        ],
    ),
    (
        "standard/1500m-25C-830hPa/rh0/no-zones/try_integrate_trajectory/RK4",
        &[
            0x0000000000000508, 0x3f50624dd2f1a9fc, 0x3fea544b088d78ca, 0x3f5ae181bc23bf86,
            0xbeb103e9073ae28b, 0x4089b4a60821c294, 0x3ffa2a6f6d01224c, 0xbf609d43f96b5fb7,
            0x3fc49ba5e353f7d3, 0x405fd3d2b7d4161b, 0x3fc0bf5b6609227a, 0xbf99d85c73ab3a13,
            0x4087be4d421310c6, 0x3f5bff49860546c8, 0xbfd3a9831b1d4f35, 0x3fd49ba5e353f7d3,
            0x406ea1517c64f952, 0x3f7d6e021e317074, 0xbfb8d7153baa87bc, 0x4085feb4b1a17030,
            0xbff84d2b635ff729, 0xbfe288683aeaf653, 0x3fe49ba5e353f7d3, 0x407c7f844490979a,
            0xbfeddaa0e60acb06, 0xbfd707f46b586b9f, 0x40830e702c29461e, 0xc0110acf10b04db9,
            0xbff094316eae86cb, 0x3feee978d4fdf3bd, 0x4083fd2a3836b58a, 0xc0059f2e41ac94a2,
            0xbfe82a26327c3409, 0x4080ac1a5a338351, 0xc01abdc9f0a66591, 0xbff682fd3992a91c,
            0x3ff4937ef8389c4b, 0x4089000000000000, 0xc014c67fd7367cd4, 0xbff417d41a0eb491,
            0x407d68a7e9c21d7f, 0xc021b1ac426e44f3, 0xbffb689d3494e594,
        ],
    ),
    (
        "standard/1500m-25C-830hPa/rh0/no-zones/try_integrate_trajectory/RK45",
        &[
            0x000000000000002a, 0x3f9db22d0e56041a, 0x4037b156c360126a, 0x3fa6290b576b1bee,
            0xbf4bbf3407840676, 0x408958394dcfb131, 0x3ff56f8ed711d5ec, 0xbfadc9ae95dedc3a,
            0x3fc1cac083126e97, 0x405ba04b76654587, 0x3fc070131cc91103, 0xbf935efa9abf7976,
            0x4087ffcb0a779252, 0x3fcbfc3a6c826d96, 0xbfd11d9f0cbb08f1, 0x3fd27ef9db22d0e7,
            0x406bb451e88e01ed, 0x3faaca977ac0cd26, 0xbfb42bbf003779cc, 0x4086561cd699f288,
            0xbff375eb4f02e25c, 0xbfe0d563dbae7a0c, 0x3fe3ced916872b05, 0x407b8a4827e45375,
            0xbfea8690b824ee31, 0xbfd56611d64aea8d, 0x40834318c77f5989, 0xc0103d7f68ad928a,
            0xbff01127743eb46f, 0x3fed6872b020c4a1, 0x408332a59c6c733c, 0xc0032bd9c1096a3f,
            0xbfe615f84540b7f8, 0x4080fe0dc6678a2b, 0xc01964bc685eda17, 0xbff5b70eba5a234e,
            0x3ff493832b7e427b, 0x4089000000000000, 0xc014c69ac39c7bc0, 0xbff417e4e7256198,
            0x407d68a92348a836, 0xc021b1acfec8e619, 0xbffb689bae7bc0ae,
        ],
    ),
    (
        "standard/1500m-25C-830hPa/rh0/no-zones/try_solve_trajectory_rust/RK45",
        &[
            0x000000000000002a, 0x3f9db22d0e56041a, 0x4037b156c360126a, 0x3fa6290b576b1bee,
            0xbf4bbf3407840676, 0x408958394dcfb131, 0x3ff56f8ed711d5ec, 0xbfadc9ae95dedc3a,
            0x3fc1cac083126e97, 0x405ba04b76654587, 0x3fc070131cc91103, 0xbf935efa9abf7976,
            0x4087ffcb0a779252, 0x3fcbfc3a6c826d96, 0xbfd11d9f0cbb08f1, 0x3fd27ef9db22d0e7,
            0x406bb451e88e01ed, 0x3faaca977ac0cd26, 0xbfb42bbf003779cc, 0x4086561cd699f288,
            0xbff375eb4f02e25c, 0xbfe0d563dbae7a0c, 0x3fe3ced916872b05, 0x407b8a4827e45375,
            0xbfea8690b824ee31, 0xbfd56611d64aea8d, 0x40834318c77f5989, 0xc0103d7f68ad928a,
            0xbff01127743eb46f, 0x3fed6872b020c4a1, 0x408332a59c6c733c, 0xc0032bd9c1096a3f,
            0xbfe615f84540b7f8, 0x4080fe0dc6678a2b, 0xc01964bc685eda17, 0xbff5b70eba5a234e,
            0x3ff493832b7e427b, 0x4089000000000000, 0xc014c69ac39c7bc0, 0xbff417e4e7256198,
            0x407d68a92348a836, 0xc021b1acfec8e619, 0xbffb689bae7bc0ae,
        ],
    ),
    (
        "standard/1500m-25C-830hPa/rh0/no-zones/fast_integrate",
        &[
            0x0000000000000001, 0x0000000000000508, 0x3f50624dd2f1a9fc, 0x3fea544b088d78cc,
            0x3f5ae181bc23bf88, 0xbeb103e8f0646847, 0x4089b4a60821c298, 0x3ffa2a6f6d012250,
            0xbf609d43e31eaa5a, 0x3fc49ba5e353f7cf, 0x405fd3d2b7d4178d, 0x3fc0bf5b66092374,
            0xbf99d85c50faca12, 0x4087be4d421312e5, 0x3f5bff498609b7f0, 0xbfd3a98300b95ce2,
            0x3fd49ba5e353f7cf, 0x406ea1517c64fbfa, 0x3f7d6e021e31ad80, 0xbfb8d7151a53697a,
            0x4085feb4b1a17411, 0xbff84d2b635ff723, 0xbfe28868220b0d7c, 0x3fe49ba5e353f7cf,
            0x407c7f8444909c3e, 0xbfeddaa0e60acc63, 0xbfd707f44c6eea10, 0x40830e702c294c56,
            0xc0110acf10b04f4a, 0xbff09431586e01ee, 0x3feee978d4fdf3b6, 0x4083fd2a3836ba08,
            0xc0059f2e41ac969a, 0xbfe82a26120d3535, 0x4080ac1a5a338a80, 0xc01abdc9f0a66913,
            0xbff682fd1b5ba38a, 0x3ff4937ef8389561, 0x4089000000000000, 0xc014c67fd7366ec4,
            0xbff417d3ff16a8df, 0x407d68a7e9c231f0, 0xc021b1ac426e41c6, 0xbffb689d0fcb4247,
            0x0000000000000001, 0x3ff4937ef8389561, 0x0000000000000001, 0x3fc49ba5e353f7cf,
            0x0000000000000000,
        ],
    ),
    (
        "standard/1500m-25C-830hPa/rh0/no-zones/fast_integrate_with_segments",
        &[
            0x0000000000000001, 0x000000000000002c, 0x3f947ae147ae147b, 0x4030607d2adb373c,
            0x3f9f89ceb12401d2, 0xbf3a74a45ba75f73, 0x408975b688aabc52, 0x3ff6f2d52b2cd838,
            0xbfa49d6f7ac5a24c, 0x3fbeb851eb851eb7, 0x4057f6153a1ded8f, 0x3fbf5a45e17d9cca,
            0xbf8d03571452f22b, 0x4088393cb8c10198, 0x3fda1aa398736d21, 0xbfcdc3a9df228283,
            0x3fd147ae147ae149, 0x4069ffcc55338fd4, 0x3fb2e23a14d4d2ae, 0xbfb1b014ca71dcb0,
            0x4086896c9a08d8a0, 0xbff0a40274e35094, 0xbfdfac08525dde52, 0x3fe3333333333336,
            0x407ace21616227d7, 0xbfe81a91e7fce853, 0xbfd4314d17972723, 0x40836bb24b2ec493,
            0xc00f404a2de892a8, 0xbfef583f2d12e90c, 0x3fedc28f5c28f5c8, 0x408362640c481d9f,
            0xc003bbc2dc11f158, 0xbfe690ccd33b39fa, 0x4080eaacb681d8f5, 0xc019b5fddc93afa0,
            0xbff5e7482484f6d0, 0x3ff49383c673d15c, 0x4089000000000000, 0xc014c69ea5e28963,
            0xbff417e752fba05f, 0x407d68a950eb4191, 0xc021b1ad1a48bec0, 0xbffb689b75b3e689,
            0x0000000000000001, 0x3ff49383c673d15c, 0x0000000000000001, 0x3fc4a307b3bb46f1,
            0x0000000000000000,
        ],
    ),
    (
        "standard/1500m-25C-830hPa/rh0/zones-rh0-rh0/compute_derivatives",
        &[
            0x4085e00000000000, 0x4000000000000000, 0x3fe0000000000000, 0xc074144820ee96dc,
            0xc02572f931e1808a, 0xc00085b7b1c5c9ab, 0x4085e00000000000, 0x4000000000000000,
            0x3fe0000000000000, 0xc0747a4a709b9a8b, 0xc0257c4cc9fd05ff, 0xc000d9a80abd7ac8,
            0x4085e00000000000, 0x4000000000000000, 0x3fe0000000000000, 0xc073e158e36cfa79,
            0xc0256e5109261bee, 0xc0005bce432f4035,
        ],
    ),
    (
        "standard/1500m-25C-830hPa/rh0/zones-rh0-rh0/try_integrate_trajectory/RK4",
        &[
            0x0000000000000506, 0x3f50624dd2f1a9fc, 0x3fea544edb8ca76b, 0x3f5ae185a4a48e5d,
            0xbeb0ddd83e20de09, 0x4089b4ad7fc01c14, 0x3ffa2a770cac3529, 0xbf607819fe6f04de,
            0x3fc47ae147ae147f, 0x405fa5f122f5e508, 0x3fc0c110fc04b248, 0xbf99512c71c264cc,
            0x4087c56594ddf185, 0x3f88afb0196620e5, 0xbfd362e3b2d076c8, 0x3fd48b4395810629,
            0x406e8fef362225c3, 0x3f820df7811e882d, 0xbfb87df9b84df297, 0x408608bdae5f2b07,
            0xbff827c02b6093a2, 0xbfe2567693565b1b, 0x3fe49374bc6a7efe, 0x407c7e5994bdb849,
            0xbfedb93ef176c5a0, 0xbfd6cc39be59fac5, 0x408318da051bcd55, 0xc01102d8ca8d51f6,
            0xbff07a4792793711, 0x3feed916872b0213, 0x4083fa67e2766a7a, 0xc0058212d119f9a4,
            0xbfe7f6585af93349, 0x4080b0f1c8d01720, 0xc01aa82ff1614e4a, 0xbff676f0f1f9dad6,
            0x3ff48d96911cf141, 0x4089000000000000, 0xc014b505b91ba8a5, 0xbff400327d805086,
            0x407d663a879ac172, 0xc021a557b1cd1abe, 0xbffb6ba259940aeb,
        ],
    ),
    (
        "standard/1500m-25C-830hPa/rh0/zones-rh0-rh0/try_integrate_trajectory/RK45",
        &[
            0x000000000000002a, 0x3f9db22d0e56041a, 0x4037b1b9ca3ced02, 0x3fa6296ae8841c90,
            0xbf4b8198dd4db814, 0x4089590d09138371, 0x3ff57054e3c27011, 0xbfad87d23d5156b2,
            0x3fc1cac083126e97, 0x405ba2659c46f833, 0x3fc071a2088e4b39, 0xbf933521516972ff,
            0x4088037221e30254, 0x3fcc0d9dd64c2130, 0xbfd0f9434a849dd5, 0x3fd27ef9db22d0e7,
            0x406bb88d53445652, 0x3faad959c54e3d32, 0xbfb4019ed5860440, 0x40865cf0ec74dfad,
            0xbff375185d049254, 0xbfe0b367c51b3383, 0x3fe3ced916872b05, 0x407b927740b38ee3,
            0xbfea87fbc9efca99, 0xbfd53d56b2e4d321, 0x40834c1eed5fd860, 0xc0103e3cb231f207,
            0xbfeff5657a276a8c, 0x3fed6872b020c4a1, 0x4083385d14a1582f, 0xc0032ab49cd9c63e,
            0xbfe5f984b1922782, 0x4081004d272af615, 0xc0195ecf7fd348ed, 0xbff5b176ebf29b94,
            0x3ff48d92e8f87266, 0x4089000000000000, 0xc014b51782f770df, 0xbff40023dcee3f86,
            0x407d66536b80104a, 0xc021a55448747602, 0xbffb6b836177cc4d,
        ],
    ),
    (
        "standard/1500m-25C-830hPa/rh0/zones-rh0-rh0/fast_integrate",
        &[
            0x0000000000000001, 0x0000000000000506, 0x3f50624dd2f1a9fc, 0x3fea544edb8ca76d,
            0x3f5ae185a4a48e5f, 0xbeb0ddd8277d7b22, 0x4089b4ad7fc01c18, 0x3ffa2a770cac352d,
            0xbf607819e8543117, 0x3fc47ae147ae147b, 0x405fa5f122f5e65e, 0x3fc0c110fc04b33b,
            0xbf99512c4fc76767, 0x4087c56594ddf391, 0x3f88afb01966ab51, 0xbfd362e398cb4e74,
            0x3fd48b4395810625, 0x406e8fef3622284f, 0x3f820df7811ea5da, 0xbfb87df9976e6db9,
            0x408608bdae5f2ece, 0xbff827c02b6093a3, 0xbfe256767ab97b06, 0x3fe49374bc6a7efa,
            0x407c7e5994bdbce0, 0xbfedb93ef176c713, 0xbfd6cc399fc0a431, 0x408318da051bd37f,
            0xc01102d8ca8d538f, 0xbff07a477c5b7a09, 0x3feed916872b020c, 0x4083fa67e2766eee,
            0xc0058212d119fba8, 0xbfe7f6583acfbc55, 0x4080b0f1c8d01e53, 0xc01aa82ff16151d4,
            0xbff676f0d3d300e5, 0x3ff48d96911cea5e, 0x4089000000000000, 0xc014b505b91b9abf,
            0xbff4003262a7fc9d, 0x407d663a879ad5dd, 0xc021a557b1cd1798, 0xbffb6ba234c659d2,
            0x0000000000000001, 0x3ff48d96911cea5e, 0x0000000000000001, 0x3fc49ba5e353f7cf,
            0x0000000000000000,
        ],
    ),
    (
        "standard/1500m-25C-830hPa/rh0/zones-rh0-rh0/fast_integrate_with_segments",
        &[
            0x0000000000000001, 0x000000000000002c, 0x3f947ae147ae147b, 0x403060ac7f60ba61,
            0x3f9f8a2bc67a491e, 0xbf3a39c09b27754b, 0x40897649a52694f7, 0x3ff6f362d2e6f289,
            0xbfa46face4dfa87a, 0x3fbeb851eb851eb7, 0x4057f7aa0a0344b4, 0x3fbf5cbd537a4fc0,
            0xbf8cc46153b2cb15, 0x40883c6f830a42f8, 0x3fda230213ad293a, 0xbfcd83fda535eb02,
            0x3fd147ae147ae149, 0x406a03866b457b37, 0x3fb2e94d35de407c, 0xbfb18afbf626345a,
            0x40868fe381d3dc52, 0xbff0a2d62493e538, 0xbfdf6bafa4edcb9c, 0x3fe3333333333336,
            0x407ad5f714f492a6, 0xbfe81bdb7fba13cf, 0xbfd40a4ede481587, 0x4083754b022ebd53,
            0xc00f42975be69beb, 0xbfef287c4af520e5, 0x3fedc28f5c28f5c8, 0x408368252d90399b,
            0xc003ba80b350c270, 0xbfe674292b3e1def, 0x4080ecc12e51e3b3, 0xc019afd3f46e2b98,
            0xbff5e21b1dbf711b, 0x3ff48d8c00fcc464, 0x4089000000000000, 0xc014b50cb67172b5,
            0xbff400083cff6258, 0x407d66665f64cf36, 0xc021a550dfdd7ed8, 0xbffb6b6bcc976d48,
            0x0000000000000001, 0x3ff48d8c00fcc464, 0x0000000000000001, 0x3fc4a4ce753fc5cb,
            0x0000000000000000,
        ],
    ),
    (
        "standard/1500m-25C-830hPa/rh0/solve_trajectory_for_monte_carlo",
        &[
            0x4014f9b30a69a1f7, 0xbff417d3ff16a8df, 0x3ff4937ef8389561, 0x407d6a0951259e5d,
            0x40939f64d7cb8015, 0x3ff5c0f70ff05668, 0x0000000000000000, 0x4089000000000000,
        ],
    ),
];

/// One in-process call through each JSON front end, RH 0: solve-json decoded and run by `solve_v1`
/// (Euler at sea level, RK45 at 1500 m), and `card_service::come_ups_v1` at both stations.
///
/// IMMUTABLE: no later step may change it, S8 included.
#[rustfmt::skip]
const FRONT_ENDS_DRY: Pinned = &[
    (
        "solve_v1/euler/icao-sea-level/rh0",
        &[
            0x4082c00000000000, 0x3ff436ef709c1d9e, 0x3fee4b889ee731a0, 0x407e90198a455bdc,
            0x40952f671003d2f8, 0x40045740e6bb5acd, 0x3fbb8509cf43f5f8, 0x3fdb6fe38aae2234,
            0x3f54ccccccccccce, 0x000000000000000d, 0x4049000000000000, 0x3fafba77b3a2b9ca,
            0x4088bc2927236a2e, 0xbf4ae168fd534800, 0xbf70ed01b044652e, 0x40029b9faa4040e9,
            0x4072c00000000000, 0x3fda62134285c8d4, 0x408422d7573ed223, 0x3fd9cb4e1bc474d2,
            0xbfc56ce710c032fa, 0x3ffe4bdffe0c7288, 0x4082c00000000000, 0x3fee4b889ee731a0,
            0x407e90198a455bdc, 0x40078aafe32b2ab9, 0xbfe86c4bb8892959, 0x3ff6fdfebd41b9ce,
            0x302e890e74b9b523,
        ],
    ),
    (
        "solve_v1/rk45/1500m-25C-830hPa/rh0",
        &[
            0x4082c00000000000, 0x3ff43d1c6529737f, 0x3fec9187d9e44b66, 0x40812d41801d439e,
            0x409ac42523670cde, 0x4009b189fc2b839c, 0x3fbd22f26ba06303, 0x3fdf44c41c7ffa25,
            0x3f54933333333334, 0x000000000000000d, 0x4049000000000000, 0x3faf97acb8d0b8be,
            0x4088eecd8b46edd3, 0xbf571d8b8a979000, 0xbf68225767add3c3, 0x400270825efdf73e,
            0x4072c00000000000, 0x3fd9b28549cf850d, 0x408534d2f4eb9ae1, 0x3fd7d9f0a65bdeb6,
            0xbfbf232ef28015bb, 0x3fff5df4d02347f0, 0x4082c00000000000, 0x3fec9187d9e44b66,
            0x40812d41801d439e, 0x400526318b97f1cc, 0xbfe1508b88e5cab0, 0x3ff968210bd57583,
            0x7acbebacbcb90cf3,
        ],
    ),
    (
        "card.come_ups/icao-sea-level/rh0",
        &[
            0x0000000000000008, 0x0000000000000001, 0x4059000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x3f18948ea9ea5c2a, 0x0000000000000001,
            0x3f18948ea9ea5c2a, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a3a1cdfc64981f, 0x0000000000000001,
            0x40a32aa762ec026e, 0x0000000000000001, 0x3fbd7c91dcc992f7, 0x0000000000000001,
            0x4069000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x3fe1d35f10f796b6, 0x0000000000000001, 0x3fe1d29a6c824763, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40a23a57948b651d, 0x0000000000000001, 0x40a085df965baef3, 0x0000000000000001,
            0x3fce9a3e43bef706, 0x0000000000000001, 0x4072c00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x3ff55bca427d9128, 0x0000000000000001,
            0x3fe8e43574038b9a, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a0e2584666720f, 0x0000000000000001,
            0x409c5a6f89ac6be1, 0x0000000000000001, 0x3fd7da11cc019af2, 0x0000000000000001,
            0x4079000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4001efec184a9293, 0x0000000000000001, 0x3fed081bdc2f27fc, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x409f33672d656f29, 0x0000000000000001, 0x409834a0bee6e2b7, 0x0000000000000001,
            0x3fe08c263d3c0f29, 0x0000000000000001, 0x407f400000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x400a2c58f2f11743, 0x0000000000000001,
            0x3ff078d9b54d0960, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x409cbcf08b7f5834, 0x0000000000000001,
            0x40948901f2694137, 0x0000000000000001, 0x3fe58e89c503c7da, 0x0000000000000001,
            0x4082c00000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4011bd907259b53e, 0x0000000000000001, 0x3ff29d8fe384a672, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x409a5e5390598680, 0x0000000000000001, 0x409149d1fa743d8e, 0x0000000000000001,
            0x3feb0185431ac291, 0x0000000000000001, 0x4085e00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x4017021207dd6314, 0x0000000000000001,
            0x3ff51206560eb758, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x409816099cc7171d, 0x0000000000000001,
            0x408cd99cecec92c8, 0x0000000000000001, 0x3ff07aa83c1d35d4, 0x0000000000000001,
            0x4089000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x401cff788a29e696, 0x0000000000000001, 0x3ff7f59a09320e08, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4095e4cacf499f7f, 0x0000000000000001, 0x4087d64bb5e407e5, 0x0000000000000001,
            0x3ff3bedfca51ca17,
        ],
    ),
    (
        "card.come_ups/1500m-77F-24.51inHg/rh0",
        &[
            0x0000000000000008, 0x0000000000000001, 0x4059000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x3f4738c20aa8c9dc, 0x0000000000000001,
            0x3f4738c20aa8c9dc, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a3ec470b899053, 0x0000000000000001,
            0x40a3bd2555740d4e, 0x0000000000000001, 0x3fbd451206cbf144, 0x0000000000000001,
            0x4069000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x3fe11fc4f82e741f, 0x0000000000000001, 0x3fe119f6c7abc9ed, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40a2c9f68fa9ed3b, 0x0000000000000001, 0x40a18e41485376aa, 0x0000000000000001,
            0x3fce2420d7c8ecdb, 0x0000000000000001, 0x4072c00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x3ff4527fb00bfc51, 0x0000000000000001,
            0x3fe7853a67e98483, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a1b16871e5acd9, 0x0000000000000001,
            0x409f22872f876cb6, 0x0000000000000001, 0x3fd74c369c08e4dd, 0x0000000000000001,
            0x4079000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4000dee57a05fa8d, 0x0000000000000001, 0x3fead69687fff192, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40a0a2e9f7f737f0, 0x0000000000000001, 0x409b86f651a26401, 0x0000000000000001,
            0x3fe0053f5175ae87, 0x0000000000000001, 0x407f400000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40084f3155edf096, 0x0000000000000001,
            0x3fedc12f6f9fd824, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x409f3b8bdb0402b0, 0x0000000000000001,
            0x40984144fb49c0e4, 0x0000000000000001, 0x3fe4ac827961d3e1, 0x0000000000000001,
            0x4082c00000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40103f1ad73eabb1, 0x0000000000000001, 0x3ff05e08b11ecd98, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x409d41d127604ee8, 0x0000000000000001, 0x4095489f16ad10dc, 0x0000000000000001,
            0x3fe9a2d6b78d8cab, 0x0000000000000001, 0x4085e00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x4014bf21146811cc, 0x0000000000000001,
            0x3ff20018f4a5986c, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x409b5752c152621a, 0x0000000000000001,
            0x409296596f48a2d4, 0x0000000000000001, 0x3feef0a6a74fffe7, 0x0000000000000001,
            0x4089000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4019b36412076420, 0x0000000000000001, 0x3ff3d10bf67d4950, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40997b176c82709a, 0x0000000000000001, 0x409024dd4cbb0fae, 0x0000000000000001,
            0x3ff24fce3c88b7f9,
        ],
    ),
];

/// As `FRONT_ENDS_DRY` with humidity: solve-json at RH 0.5 (RK4, sea level) and 0.9 (RK45, 1500 m),
/// and the come-ups card at 50 and 100 %.
///
/// MAY CHANGE ONLY IN S8 (the unified aloft vapour model, plan C7): byte-identical through every
/// other step.
#[rustfmt::skip]
const FRONT_ENDS_HUMID: Pinned = &[
    (
        "solve_v1/rk4/icao-sea-level/rh0.5",
        &[
            0x4082c00000000000, 0x3ff43727c9bba23e, 0x3fee41cef59a48ed, 0x407e9f363f6945dc,
            0x4095445f8ed12820, 0x40045740e6bb5acd, 0x3fbb74e145b1de64, 0x3fdb6fe38aae2234,
            0x3f54b33333333334, 0x000000000000000d, 0x4049000000000000, 0x3fafb7608ac85fa7,
            0x4088bcfee5a9c210, 0xbf4b462d19c85000, 0xbf708aa51ebaefd6, 0x400295e8654c19b8,
            0x4072c00000000000, 0x3fda5d20c9d62b88, 0x4084274f4f5af2d0, 0x3fd9b909df52d576,
            0xbfc5465b660a619d, 0x3ffe48430655bf6a, 0x4082c00000000000, 0x3fee41cef59a48ed,
            0x407e9f363f6945dc, 0x400779d2f09c1656, 0xbfe8476a31d2fb44, 0x3ff70182cd387405,
            0xc73ada8bf3458070,
        ],
    ),
    (
        "solve_v1/rk45/1500m-25C-830hPa/rh0.9",
        &[
            0x4082c00000000000, 0x3ff43d363a70ba00, 0x3fec8053c33c40bd, 0x408141d3b02ab805,
            0x409b0467e466a43f, 0x4009b189fc2b839c, 0x3fbd02deb725f8b5, 0x3fdf44c41c7ffa25,
            0x3f54933333333334, 0x000000000000000d, 0x4049000000000000, 0x3faf965ac2a1e71c,
            0x4088f0e3ef9e0f75, 0xbf572c7ead4cd400, 0xbf67ce5a2c57d452, 0x400258972e416488,
            0x4072c00000000000, 0x3fd9ab83fcb39f3b, 0x40854070c07f767c, 0x3fd7c8af9312334c,
            0xbfbeb66fff324261, 0x3fff43bf091de1c8, 0x4082c00000000000, 0x3fec8053c33c40bd,
            0x408141d3b02ab805, 0x40050f53d587d7aa, 0xbfe10fbda4d4ed5e, 0x3ff96352050f72ac,
            0xa8aa06f7356f134b,
        ],
    ),
    (
        "card.come_ups/icao-sea-level/rh50",
        &[
            0x0000000000000008, 0x0000000000000001, 0x4059000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x3ed40ef2e12ce85f, 0x0000000000000001,
            0x3ed40ef2e12ce85f, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a3a2be8e4171fc, 0x0000000000000001,
            0x40a32c7d29c24146, 0x0000000000000001, 0x3fbd7bdce0e72d8d, 0x0000000000000001,
            0x4069000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x3fe1d048482e2d8c, 0x0000000000000001, 0x3fe1d03e40b4bcf6, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40a23c28139ffdfe, 0x0000000000000001, 0x40a08929d99cbffc, 0x0000000000000001,
            0x3fce98b93eba93e9, 0x0000000000000001, 0x4072c00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x3ff557d4d00962c0, 0x0000000000000001,
            0x3fe8df6157e497f4, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a0e4f752ab5c9c, 0x0000000000000001,
            0x409c633dfa817a92, 0x0000000000000001, 0x3fd7d83712ec3710, 0x0000000000000001,
            0x4079000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4001ec19612cb27c, 0x0000000000000001, 0x3fed00bbe4a00470, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x409f3a159874318a, 0x0000000000000001, 0x40983effbda63151, 0x0000000000000001,
            0x3fe08a5eee8d5a50, 0x0000000000000001, 0x407f400000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x400a25a9cc1c03ac, 0x0000000000000001,
            0x3ff07320d5dea260, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x409cc4e5c9d54020, 0x0000000000000001,
            0x409494630154a77b, 0x0000000000000001, 0x3fe58b8674affc74, 0x0000000000000001,
            0x4082c00000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4011b8377cfbab86, 0x0000000000000001, 0x3ff2958a5bb6a6c0, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x409a676f8e1f5beb, 0x0000000000000001, 0x409155c5f5cee53e, 0x0000000000000001,
            0x3feafccfbe2409e0, 0x0000000000000001, 0x4085e00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x4016fa036158e270, 0x0000000000000001,
            0x3ff5072f9174dba8, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x4098202fff46103d, 0x0000000000000001,
            0x408cf1f2aed2bc00, 0x0000000000000001, 0x3ff0772bdbc30521, 0x0000000000000001,
            0x4089000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x401cf3827e00e1e9, 0x0000000000000001, 0x3ff7e5fc729ffde4, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4095efda2e96ecd5, 0x0000000000000001, 0x4087ee672f907fd2, 0x0000000000000001,
            0x3ff3b9dd686a5aae,
        ],
    ),
    (
        "card.come_ups/icao-sea-level/rh100",
        &[
            0x0000000000000008, 0x0000000000000001, 0x4059000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0xbf160e4b708c60f4, 0x0000000000000001,
            0xbf160e4b708c60f4, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a3a3aefa4fb39e, 0x0000000000000001,
            0x40a32e52bd3b3c0c, 0x0000000000000001, 0x3fbd7b280c465958, 0x0000000000000001,
            0x4069000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x3fe1cd322ed321d0, 0x0000000000000001, 0x3fe1cde2a12ea633, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40a23df8538ff85b, 0x0000000000000001, 0x40a08c73fe20b1c8, 0x0000000000000001,
            0x3fce9734a304edc4, 0x0000000000000001, 0x4072c00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x3ff553e06156a24b, 0x0000000000000001,
            0x3fe8da8e93da22c6, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a0e796082c9ad3, 0x0000000000000001,
            0x409c6c0ca57a9813, 0x0000000000000001, 0x3fd7d65cf6adad6e, 0x0000000000000001,
            0x4079000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4001e847c70111be, 0x0000000000000001, 0x3fecf95e59570262, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x409f40c33ef75919, 0x0000000000000001, 0x4098495fc33a6663, 0x0000000000000001,
            0x3fe0889850082692, 0x0000000000000001, 0x407f400000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x400a1efcd69a3a69, 0x0000000000000001,
            0x3ff06d6a1f325156, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x409cccda686bc2a7, 0x0000000000000001,
            0x40949fc651b0410c, 0x0000000000000001, 0x3fe588847796dab9, 0x0000000000000001,
            0x4082c00000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4011b2e074cedfd6, 0x0000000000000001, 0x3ff28d8826070a86, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x409a708b33d9c285, 0x0000000000000001, 0x409161bd9dcc3bac, 0x0000000000000001,
            0x3feaf81c801daabc, 0x0000000000000001, 0x4085e00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x4016f1f7dbe430dc, 0x0000000000000001,
            0x3ff4fc5d9c554418, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40982a5664ce0c2b, 0x0000000000000001,
            0x408d0a52b71843f6, 0x0000000000000001, 0x3ff073b14d34e6cf, 0x0000000000000001,
            0x4089000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x401ce798fc57f2ac, 0x0000000000000001, 0x3ff7d68481cf0740, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4095faea4df95266, 0x0000000000000001, 0x4088069077b9697e, 0x0000000000000001,
            0x3ff3b4dee301a74c,
        ],
    ),
    (
        "card.come_ups/1500m-77F-24.51inHg/rh50",
        &[
            0x0000000000000008, 0x0000000000000001, 0x4059000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x3f41f7109cbc2087, 0x0000000000000001,
            0x3f41f7109cbc2087, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a3edff41424c9b, 0x0000000000000001,
            0x40a3c08dc3bbf116, 0x0000000000000001, 0x3fbd43ce2238295b, 0x0000000000000001,
            0x4069000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x3fe11a3db7caf9e0, 0x0000000000000001, 0x3fe115bff3a3cad8, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40a2cd4f4255b20c, 0x0000000000000001, 0x40a19482c63cb611, 0x0000000000000001,
            0x3fce217211460212, 0x0000000000000001, 0x4072c00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x3ff44bc2dbfb0285, 0x0000000000000001,
            0x3fe77d48002b0b2a, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a1b6498031a5b8, 0x0000000000000001,
            0x409f33b58d33a818, 0x0000000000000001, 0x3fd74900dbae4893, 0x0000000000000001,
            0x4079000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4000d891d93804bf, 0x0000000000000001, 0x3feacac1acea0df2, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40a0a9357606352d, 0x0000000000000001, 0x409b9bcf225b4421, 0x0000000000000001,
            0x3fe002386ca13674, 0x0000000000000001, 0x407f400000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x400844538e2baf0c, 0x0000000000000001,
            0x3fedaf06d3cea934, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x409f4ab3b8cab63c, 0x0000000000000001,
            0x409858d4d4e3cca8, 0x0000000000000001, 0x3fe4a77800f3b10d, 0x0000000000000001,
            0x4082c00000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x401036b3e757a2a6, 0x0000000000000001, 0x3ff0522881072c80, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x409d5351483e2069, 0x0000000000000001, 0x4095621d06e54b6c, 0x0000000000000001,
            0x3fe99b21424e3986, 0x0000000000000001, 0x4085e00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x4014b294694178ca, 0x0000000000000001,
            0x3ff1ef8207a75890, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x409b6afb1c37cd39, 0x0000000000000001,
            0x4092b11d5584f58e, 0x0000000000000001, 0x3feee57152eb41e8, 0x0000000000000001,
            0x4089000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4019a17182605127, 0x0000000000000001, 0x3ff3bb74647b6174, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x409990bdaba39aff, 0x0000000000000001, 0x40904057ce574e52, 0x0000000000000001,
            0x3ff247f86843770e,
        ],
    ),
    (
        "card.come_ups/1500m-77F-24.51inHg/rh100",
        &[
            0x0000000000000008, 0x0000000000000001, 0x4059000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x3f396a51c5003162, 0x0000000000000001,
            0x3f396a51c5003162, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a3efb7c0df3233, 0x0000000000000001,
            0x40a3c3f70fd17650, 0x0000000000000001, 0x3fbd428a3173ef91, 0x0000000000000001,
            0x4069000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x3fe114b625a75de2, 0x0000000000000001, 0x3fe11188db6ebddc, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40a2d0a8b6e070a0, 0x0000000000000001, 0x40a19ac6cc15e699, 0x0000000000000001,
            0x3fce1ec353b95211, 0x0000000000000001, 0x4072c00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x3ff44505e2ad24b0, 0x0000000000000001,
            0x3fe775559fb2eb7e, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a1bb2be824f72c, 0x0000000000000001,
            0x409f44ed6a468216, 0x0000000000000001, 0x3fd745cb6d156a06, 0x0000000000000001,
            0x4079000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4000d23e64ff3e1a, 0x0000000000000001, 0x3feabeedcea2af08, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40a0af82a97e804b, 0x0000000000000001, 0x409bb0b580bb8bbc, 0x0000000000000001,
            0x3fdffe6442f5c5f4, 0x0000000000000001, 0x407f400000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40083976ca9f6b9c, 0x0000000000000001,
            0x3fed9ce19680b608, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x409f59e0786a3c41, 0x0000000000000001,
            0x40987077b74c2d1e, 0x0000000000000001, 0x3fe4a26f14bee790, 0x0000000000000001,
            0x4082c00000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40102e4e40a650e3, 0x0000000000000001, 0x3ff0464b6d5a6c54, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x409d64d85d49659d, 0x0000000000000001, 0x40957bb45c4b58ef, 0x0000000000000001,
            0x3fe9936eceda3dd9, 0x0000000000000001, 0x4085e00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x4014a60a64768bee, 0x0000000000000001,
            0x3ff1def08f40ec2c, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x409b7eacdf781915, 0x0000000000000001,
            0x4092cc0150271b2a, 0x0000000000000001, 0x3feeda4145c3e533, 0x0000000000000001,
            0x4089000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40198fae1b05e043, 0x0000000000000001, 0x3ff3a68eda3d5154, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4099a6714ae3a27a, 0x0000000000000001, 0x40905bfaae803f80, 0x0000000000000001,
            0x3ff2402bc108abc0,
        ],
    ),
];

/// The bridge `solve` command (`bridge_call`) on the `FRONT_ENDS_DRY` requests.
///
/// IMMUTABLE: no later step may change it, S8 included.
#[cfg(feature = "bridge")]
#[rustfmt::skip]
const BRIDGE_SOLVE_DRY: Pinned = &[
    (
        "bridge-solve/euler/icao-sea-level/rh0",
        &[
            0x4082c00000000000, 0x3ff436ef709c1d9e, 0x3fee4b889ee731a0, 0x407e90198a455bdc,
            0x40952f671003d2f8, 0x40045740e6bb5acd, 0x3fbb8509cf43f5f8, 0x3fdb6fe38aae2234,
            0x3f54ccccccccccce, 0x000000000000000d, 0x4049000000000000, 0x3fafba77b3a2b9ca,
            0x4088bc2927236a2e, 0xbf4ae168fd534800, 0xbf70ed01b044652e, 0x40029b9faa4040e9,
            0x4072c00000000000, 0x3fda62134285c8d4, 0x408422d7573ed223, 0x3fd9cb4e1bc474d2,
            0xbfc56ce710c032fa, 0x3ffe4bdffe0c7288, 0x4082c00000000000, 0x3fee4b889ee731a0,
            0x407e90198a455bdc, 0x40078aafe32b2ab9, 0xbfe86c4bb8892959, 0x3ff6fdfebd41b9ce,
            0x302e890e74b9b523,
        ],
    ),
    (
        "bridge-solve/rk45/1500m-25C-830hPa/rh0",
        &[
            0x4082c00000000000, 0x3ff43d1c6529737f, 0x3fec9187d9e44b66, 0x40812d41801d439e,
            0x409ac42523670cde, 0x4009b189fc2b839c, 0x3fbd22f26ba06303, 0x3fdf44c41c7ffa25,
            0x3f54933333333334, 0x000000000000000d, 0x4049000000000000, 0x3faf97acb8d0b8be,
            0x4088eecd8b46edd3, 0xbf571d8b8a979000, 0xbf68225767add3c3, 0x400270825efdf73e,
            0x4072c00000000000, 0x3fd9b28549cf850d, 0x408534d2f4eb9ae1, 0x3fd7d9f0a65bdeb6,
            0xbfbf232ef28015bb, 0x3fff5df4d02347f0, 0x4082c00000000000, 0x3fec9187d9e44b66,
            0x40812d41801d439e, 0x400526318b97f1cc, 0xbfe1508b88e5cab0, 0x3ff968210bd57583,
            0x7acbebacbcb90cf3,
        ],
    ),
];

/// The bridge `solve` command on the `FRONT_ENDS_HUMID` requests.
///
/// MAY CHANGE ONLY IN S8 (the unified aloft vapour model, plan C7): byte-identical through every
/// other step.
#[cfg(feature = "bridge")]
#[rustfmt::skip]
const BRIDGE_SOLVE_HUMID: Pinned = &[
    (
        "bridge-solve/rk4/icao-sea-level/rh0.5",
        &[
            0x4082c00000000000, 0x3ff43727c9bba23e, 0x3fee41cef59a48ed, 0x407e9f363f6945dc,
            0x4095445f8ed12820, 0x40045740e6bb5acd, 0x3fbb74e145b1de64, 0x3fdb6fe38aae2234,
            0x3f54b33333333334, 0x000000000000000d, 0x4049000000000000, 0x3fafb7608ac85fa7,
            0x4088bcfee5a9c210, 0xbf4b462d19c85000, 0xbf708aa51ebaefd6, 0x400295e8654c19b8,
            0x4072c00000000000, 0x3fda5d20c9d62b88, 0x4084274f4f5af2d0, 0x3fd9b909df52d576,
            0xbfc5465b660a619d, 0x3ffe48430655bf6a, 0x4082c00000000000, 0x3fee41cef59a48ed,
            0x407e9f363f6945dc, 0x400779d2f09c1656, 0xbfe8476a31d2fb44, 0x3ff70182cd387405,
            0xc73ada8bf3458070,
        ],
    ),
    (
        "bridge-solve/rk45/1500m-25C-830hPa/rh0.9",
        &[
            0x4082c00000000000, 0x3ff43d363a70ba00, 0x3fec8053c33c40bd, 0x408141d3b02ab805,
            0x409b0467e466a43f, 0x4009b189fc2b839c, 0x3fbd02deb725f8b5, 0x3fdf44c41c7ffa25,
            0x3f54933333333334, 0x000000000000000d, 0x4049000000000000, 0x3faf965ac2a1e71c,
            0x4088f0e3ef9e0f75, 0xbf572c7ead4cd400, 0xbf67ce5a2c57d452, 0x400258972e416488,
            0x4072c00000000000, 0x3fd9ab83fcb39f3b, 0x40854070c07f767c, 0x3fd7c8af9312334c,
            0xbfbeb66fff324261, 0x3fff43bf091de1c8, 0x4082c00000000000, 0x3fec8053c33c40bd,
            0x408141d3b02ab805, 0x40050f53d587d7aa, 0xbfe10fbda4d4ed5e, 0x3ff96352050f72ac,
            0xa8aa06f7356f134b,
        ],
    ),
];

/// The FFI single solve `ballistics_calculate_trajectory` (RK4, sampling and spin drift on) at both
/// stations, RH 0.
///
/// IMMUTABLE: no later step may change it, S8 included.
#[cfg(feature = "ffi")]
#[rustfmt::skip]
const FFI_TRAJECTORY_DRY: Pinned = &[
    (
        "ffi-trajectory/rk4/icao-sea-level/rh0",
        &[
            0x4082c00000000000, 0x3fc087dde518e09a, 0x3fee47f6a570ee2b, 0x407e92d3d073c95e,
            0x4095332f45852df1, 0x7ff8000000000000, 0x00000000000003b4, 0x3fee47f6a570ee2b,
            0x4082c00000000000, 0xc0042c09ce29ae8c, 0xbfebce9ec312358b, 0x407e92d3d073c95e,
            0xb1decbb229999215, 0x000000000000000d, 0xacd0a11caa18e600, 0x408423b098f4d8f8,
            0x3ffe4d26df634962,
        ],
    ),
    (
        "ffi-trajectory/rk4/1500m-25C-830hPa/rh0",
        &[
            0x4082c00000000000, 0x3fc0bf5b66092374, 0x3fec91b1f79d890c, 0x40812d06af2af2e5,
            0x409ac36dd6d33524, 0x7ff8000000000000, 0x000000000000037e, 0x3fec91b1f79d890c,
            0x4082c00000000000, 0xc001dc1ed9a630e2, 0xbfe4f58c0bc4a119, 0x40812d06af2af2e5,
            0xfc54a9bf20059056, 0x000000000000000d, 0x19bb55dd301e8e81, 0x408534b005de9e8b,
            0x3fff5dc1242bacd2,
        ],
    ),
];

/// As `FFI_TRAJECTORY_DRY`, at RH 50 and 100.
///
/// MAY CHANGE ONLY IN S8 (the unified aloft vapour model, plan C7): byte-identical through every
/// other step.
#[cfg(feature = "ffi")]
#[rustfmt::skip]
const FFI_TRAJECTORY_HUMID: Pinned = &[
    (
        "ffi-trajectory/rk4/icao-sea-level/rh50",
        &[
            0x4082c00000000000, 0x3fc0888f199305b8, 0x3fee4210509fc9ad, 0x407e9ea9d3099490,
            0x4095439c83e67e19, 0x7ff8000000000000, 0x00000000000003b3, 0x3fee4210509fc9ad,
            0x4082c00000000000, 0xc00423ec3a6a9459, 0xbfebb7056fcda396, 0x407e9ea9d3099490,
            0xd5e5814389e460e4, 0x000000000000000d, 0x1c66cab1ccaee473, 0x40842724cd5814f4,
            0x3ffe480327459bb3,
        ],
    ),
    (
        "ffi-trajectory/rk4/icao-sea-level/rh100",
        &[
            0x4082c00000000000, 0x3fc089403ae45153, 0x3fee3c2d68741bb5, 0x407eaa7fbc63d0ac,
            0x4095540ff9eee471, 0x7ff8000000000000, 0x00000000000003b2, 0x3fee3c2d68741bb5,
            0x4082c00000000000, 0xc0041bd3f4825fb2, 0xbfeb9f79cf1eebeb, 0x407eaa7fbc63d0ac,
            0x6ab5fa9cd0acefd8, 0x000000000000000d, 0xecb0189c1107e69b, 0x40842a9892fcc96f,
            0x3ffe42e917aee944,
        ],
    ),
    (
        "ffi-trajectory/rk4/1500m-25C-830hPa/rh50",
        &[
            0x4082c00000000000, 0x3fc0c0a6a7b622a1, 0x3fec88205b5d26fb, 0x4081387242281158,
            0x409ae710acc5b6bc, 0x7ff8000000000000, 0x000000000000037d, 0x3fec88205b5d26fb,
            0x4082c00000000000, 0xc001cf634d67a3d7, 0xbfe4cf459ac318e4, 0x4081387242281158,
            0x1998c1cfeba5b32b, 0x000000000000000d, 0x7cd0f2328eee3885, 0x40853b234c7a909b,
            0x3fff4f1909738067,
        ],
    ),
    (
        "ffi-trajectory/rk4/1500m-25C-830hPa/rh100",
        &[
            0x4082c00000000000, 0x3fc0c1f23553721f, 0x3fec7e941121f567, 0x408143e2fdbbddae,
            0x409b0adb5577c88e, 0x7ff8000000000000, 0x000000000000037c, 0x3fec7e941121f567,
            0x4082c00000000000, 0xc001c2b0fcb1b697, 0xbfe4a91471d6524c, 0x408143e2fdbbddae,
            0x91b3d0bff2d008ef, 0x000000000000000d, 0x8d3a0f45ae36c0c3, 0x408541985640201d,
            0x3fff40ade7905a82,
        ],
    ),
];
