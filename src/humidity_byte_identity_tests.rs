//! In-crate byte-identity pins for the humidity-consistency train (0.50.0), step S1 (plan step
//! C0a): the `pub(crate)` entry points `tests/humidity_byte_identity.rs` cannot reach.
//!
//! These are the gravity-aware variants every public entry point calls with Earth gravity,
//! pinned here at the Moon's 1.62 m/s² and at 0 m/s² so a change to how gravity reaches the
//! hydrostatic lapse cannot hide behind the Earth-only public surface, plus the uncapped
//! CIPM vapour mole fraction the density-altitude refusal reads.
//!
//! The rules are those of `tests/humidity_byte_identity.rs`: pins are only ever migrated to
//! renamed APIs, never edited; each group's header says which later step may change it; kernel
//! outputs at RH > 0 are not pinned. Values were generated on 261c1eb (ballistics-engine 0.49.0)
//! on macOS (aarch64) with
//!
//! ```text
//! cargo test --lib humidity_byte_identity_tests::print_current_pin_values -- --ignored --nocapture
//! ```
//!
//! and confirmed identical in Linux (glibc) containers on aarch64 and x86_64.
//!
//! S8 (the unified aloft vapour model, plan C7) re-pinned `LOCAL_AIR_HUMID_GRAVITY`, the one
//! group here it may change; its note gives the old -> new shifts. The new bits were confirmed
//! the same way, on macOS (aarch64) and in Linux glibc containers on aarch64 and x86_64.
//!
//! # How the pins are checked, and on which platforms
//!
//! No group here can be compared with the released 0.49.0 crate at runtime, as most of
//! `tests/humidity_byte_identity.rs`'s groups are: 0.49.0 kept every one of these entry points
//! `pub(crate)` (`get_local_atmosphere_humid_with_gravity`, `water_vapor_mole_fraction_uncapped`,
//! `compute_derivatives_with_gravity`, `integrate_trajectory_with_gravity`,
//! `fast_integrate_with_gravity`, and its `Gravity` type is private). So every group is held to
//! its literal pins, as that file's `Check::Pinned` groups are:
//!
//! - bit for bit on the platforms the pins were captured and checked on ([`exact_bit_platform`]:
//!   macOS aarch64, and Linux glibc on aarch64 and x86_64), kernel standard mode at RH 0 included (it was
//!   compared to a relative 1e-15 until 0.50.0's R5; the bits never moved);
//! - on any other platform, each float within [`UNVERIFIED_LIBM_ULPS`] of its pin and each count
//!   exactly. A libm may round a `pow` or `exp` a last bit apart: FreeBSD's moved
//!   `VAPOR_MOLE_FRACTION` (three values at -10 °C) and `LOCAL_AIR_HUMID_GRAVITY` by one or two
//!   ulps even on the unchanged 0.49.0 code.

use crate::atmosphere::{
    calculate_atmosphere, get_local_atmosphere_humid_with_gravity,
    water_vapor_mole_fraction_uncapped, KernelAtmosphere, RelativeHumidity,
};
use crate::derivatives::compute_derivatives_with_gravity;
use crate::fast_trajectory::{fast_integrate_with_gravity, FastIntegrationParams};
use crate::gravity::Gravity;
use crate::trajectory_integration::{integrate_trajectory_with_gravity, TrajectoryParams};
use crate::wind::{WindSegment, WindSock};
use crate::{BallisticInputs, DragModel};
use nalgebra::Vector3;

/// One case: a name and its labelled values (each an `f64`'s bits or an exact count).
type Case = (String, Vec<(String, Pinned)>);
type Pins = &'static [(&'static str, &'static [u64])];

#[derive(Clone, Copy)]
enum Pinned {
    F(f64),
    N(u64),
}

impl Pinned {
    fn bits(self) -> u64 {
        match self {
            Pinned::F(value) => value.to_bits(),
            Pinned::N(value) => value,
        }
    }
}

/// Whether this build links a libm the pins were captured or checked on, on an architecture they
/// were checked on: Apple's on macOS aarch64, or glibc on Linux aarch64 and x86_64. The same rule
/// as `tests/humidity_byte_identity.rs`'s `Libm`, where the reasons are given.
fn exact_bit_platform() -> bool {
    cfg!(any(
        all(target_os = "macos", target_arch = "aarch64"),
        all(
            target_os = "linux",
            target_env = "gnu",
            any(target_arch = "aarch64", target_arch = "x86_64")
        )
    ))
}

/// How far, in units in the last place, a float may lie from its pin on any other platform:
/// the largest distance measured in 0.50.0's R5 on FreeBSD 14.3, NetBSD 10.1 and OpenBSD 7.8
/// (amd64), Alpine musl on aarch64 and x86_64, and Windows 11 (x86_64, UCRT). It is the
/// bound `tests/humidity_byte_identity.rs` uses, where the measurement is described.
const UNVERIFIED_LIBM_ULPS: u64 = 2;

/// The distance between two finite `f64`s in units in the last place (0 for equal values, 0.0
/// and -0.0 included).
fn ulps_apart(a: f64, b: f64) -> u64 {
    let line = |x: f64| -> i128 {
        let bits = x.to_bits() as i64;
        i128::from(if bits < 0 { i64::MIN - bits } else { bits })
    };
    (line(a) - line(b)).unsigned_abs() as u64
}

/// Every value must equal its pin bit for bit on an exact-bit platform; elsewhere each float
/// must lie within [`UNVERIFIED_LIBM_ULPS`] of it, and each count must match exactly.
fn check(group: &str, actual: &[Case], expected: Pins) {
    let exact = exact_bit_platform();
    let names: Vec<&str> = actual.iter().map(|(name, _)| name.as_str()).collect();
    let pinned: Vec<&str> = expected.iter().map(|(name, _)| *name).collect();
    assert_eq!(names, pinned, "{group}: the probe's cases no longer match");
    let mut moved = Vec::new();
    for ((name, values), (_, bits)) in actual.iter().zip(expected) {
        assert_eq!(values.len(), bits.len(), "{group}/{name}: quantity count");
        for ((label, value), &want) in values.iter().zip(bits.iter()) {
            let ok = match *value {
                Pinned::F(got) if !exact => {
                    let pin = f64::from_bits(want);
                    if got.is_nan() || pin.is_nan() {
                        got.to_bits() == want
                    } else {
                        ulps_apart(got, pin) <= UNVERIFIED_LIBM_ULPS
                    }
                }
                value => value.bits() == want,
            };
            if !ok {
                moved.push(format!(
                    "  {name} :: {label}: pinned {want:#018x}, now {:#018x}",
                    value.bits()
                ));
            }
        }
    }
    let rule = if exact {
        "bit for bit".to_string()
    } else {
        format!("within {UNVERIFIED_LIBM_ULPS} ulps (a libm the pins were not captured on)")
    };
    assert!(
        moved.is_empty(),
        "{group}: {} pins moved ({rule}):\n{}",
        moved.len(),
        moved.join("\n")
    );
}

fn print(name: &str, cases: &[Case]) {
    println!("#[rustfmt::skip]");
    println!("const {name}: Pins = &[");
    for (case, values) in cases {
        println!("    (");
        println!("        {case:?},");
        println!("        &[");
        let bits: Vec<u64> = values.iter().map(|(_, value)| value.bits()).collect();
        for chunk in bits.chunks(4) {
            let line: Vec<String> = chunk.iter().map(|b| format!("{b:#018x}")).collect();
            println!("            {},", line.join(", "));
        }
        println!("        ],");
        println!("    ),");
    }
    println!("];");
    println!();
}

fn gravities() -> [(&'static str, Gravity); 2] {
    [
        ("g1.62", Gravity::new(1.62).expect("lunar gravity")),
        ("g0", Gravity::new(0.0).expect("zero gravity")),
    ]
}

const BASES: [(f64, f64, f64); 2] = [(0.0, 15.0, 1013.25), (1500.0, 25.0, 830.0)];

/// `get_local_atmosphere_humid_with_gravity` at non-Earth gravity, base ratio from the station's
/// own CIPM density at that RH.
fn local_air(rhs: &[f64]) -> Vec<Case> {
    let mut cases = Vec::new();
    for (gravity_name, gravity) in gravities() {
        for (base_alt, base_t, base_p) in BASES {
            for &rh in rhs {
                let (density, _) = calculate_atmosphere(base_alt, Some(base_t), Some(base_p), rh);
                let ratio = density / 1.225;
                let mut values = Vec::new();
                for offset in [-200.0, 10.0, 500.0, 3000.0, 12_000.0] {
                    let (rho, c) = get_local_atmosphere_humid_with_gravity(
                        base_alt + offset,
                        base_alt,
                        base_t,
                        base_p,
                        ratio,
                        rh,
                        gravity,
                    );
                    values.push((format!("{offset:+}m.density"), Pinned::F(rho)));
                    values.push((format!("{offset:+}m.speed_of_sound"), Pinned::F(c)));
                }
                cases.push((
                    format!("{gravity_name}/{base_alt}m/{base_t}C/{base_p}hPa/rh{rh}"),
                    values,
                ));
            }
        }
    }
    cases
}

fn probe_vapor_mole_fraction() -> Vec<Case> {
    let mut values = Vec::new();
    for (t, p) in [
        (-10.0, 1013.25),
        (15.0, 1013.25),
        (35.0, 850.0),
        (95.0, 500.0),
    ] {
        for rh in [0.0, 37.5, 50.0, 100.0] {
            values.push((
                format!("x_v({t}, {p}, {rh})"),
                Pinned::F(water_vapor_mole_fraction_uncapped(t, p, rh)),
            ));
        }
    }
    vec![("water_vapor_mole_fraction_uncapped".to_string(), values)]
}

const MASS_KG: f64 = 0.01134;
const DIAMETER_M: f64 = 0.00782;
const BC_G7: f64 = 0.243;
const MUZZLE_VELOCITY_MPS: f64 = 823.0;
const ANGLE_RAD: f64 = 0.002;

fn kernel_inputs(base_alt: f64, base_t: f64, base_p: f64) -> BallisticInputs {
    BallisticInputs {
        bc_value: BC_G7,
        bc_type: DragModel::G7,
        bullet_mass: MASS_KG,
        bullet_diameter: DIAMETER_M,
        bullet_length: 0.031,
        muzzle_velocity: MUZZLE_VELOCITY_MPS,
        twist_rate: 10.0,
        muzzle_angle: ANGLE_RAD,
        altitude: base_alt,
        temperature: base_t,
        pressure: base_p,
        wind_speed: 4.0,
        wind_angle: std::f64::consts::FRAC_PI_2,
        ..BallisticInputs::default()
    }
}

fn launch_state() -> [f64; 6] {
    [
        0.0,
        0.0,
        0.0,
        MUZZLE_VELOCITY_MPS * ANGLE_RAD.cos(),
        MUZZLE_VELOCITY_MPS * ANGLE_RAD.sin(),
        0.0,
    ]
}

fn wind_segments() -> Vec<WindSegment> {
    vec![WindSegment::new(14.4, 90.0, 2000.0)]
}

/// The pinned kernel atmospheres were written as legacy tuples; S3 moved the kernels onto
/// `KernelAtmosphere`, and `from_legacy_params` at 0 % RH is the one-line migration (a dry
/// station keeps its tuple's ratio verbatim; the direct sentinel becomes uniform air).
fn kernel_atmosphere(atmos_params: (f64, f64, f64, f64)) -> KernelAtmosphere {
    KernelAtmosphere::from_legacy_params(atmos_params, RelativeHumidity::DRY)
        .expect("pinned atmospheres are valid")
}

fn trajectory_params(atmos_params: (f64, f64, f64, f64)) -> TrajectoryParams {
    TrajectoryParams {
        mass_kg: MASS_KG,
        bc: BC_G7,
        drag_model: DragModel::G7,
        wind_segments: wind_segments(),
        atmosphere: kernel_atmosphere(atmos_params),
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
        bullet_length: 0.031,
        twist_rate: 10.0,
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

fn record_states(values: &mut Vec<(String, Pinned)>, times: &[f64], states: &[[f64; 6]]) {
    let n = times.len();
    values.push(("len".to_string(), Pinned::N(n as u64)));
    for index in [1, n / 4, n / 2, (3 * n) / 4, n - 1] {
        values.push((format!("[{index}].t"), Pinned::F(times[index])));
        for (component, value) in ["x", "y", "z", "vx", "vy", "vz"].iter().zip(states[index]) {
            values.push((format!("[{index}].{component}"), Pinned::F(value)));
        }
    }
}

/// Every gravity-aware kernel once for an atmosphere tuple, at each non-Earth gravity.
fn kernels(
    prefix: &str,
    atmos_params: (f64, f64, f64, f64),
    inputs: &BallisticInputs,
) -> Vec<Case> {
    let mut cases = Vec::new();
    for (gravity_name, gravity) in gravities() {
        let mut values = Vec::new();
        for (label, position) in [
            ("muzzle", Vector3::new(0.0, 0.0, 0.0)),
            ("300m-above", Vector3::new(600.0, 300.0, 0.0)),
        ] {
            let derivative = compute_derivatives_with_gravity(
                position,
                Vector3::new(700.0, 2.0, 0.5),
                inputs,
                Vector3::new(0.0, 0.0, -4.0),
                &kernel_atmosphere(atmos_params),
                BC_G7,
                None,
                0.0,
                None,
                gravity,
            );
            for (component, value) in ["vx", "vy", "vz", "ax", "ay", "az"].iter().zip(derivative) {
                values.push((format!("{label}.{component}"), Pinned::F(value)));
            }
        }
        cases.push((
            format!("{prefix}/{gravity_name}/compute_derivatives_with_gravity"),
            values,
        ));

        let trajectory = integrate_trajectory_with_gravity(
            launch_state(),
            (0.0, 3.0),
            trajectory_params(atmos_params),
            "RK45",
            1e-6,
            0.01,
            gravity,
        );
        let times: Vec<f64> = trajectory.iter().map(|(t, _)| *t).collect();
        let states: Vec<[f64; 6]> = trajectory
            .iter()
            .map(|(_, state)| std::array::from_fn(|component| state[component]))
            .collect();
        let mut values = Vec::new();
        record_states(&mut values, &times, &states);
        cases.push((
            format!("{prefix}/{gravity_name}/integrate_trajectory_with_gravity/RK45"),
            values,
        ));

        let solution = fast_integrate_with_gravity(
            inputs,
            &WindSock::new(wind_segments()),
            FastIntegrationParams {
                horiz: 800.0,
                vert: 0.0,
                initial_state: launch_state(),
                t_span: (0.0, 30.0),
                atmosphere: kernel_atmosphere(atmos_params),
                atmo_sock: None,
            },
            gravity,
        );
        let states: Vec<[f64; 6]> = (0..solution.t.len())
            .map(|i| std::array::from_fn(|component| solution.y[component][i]))
            .collect();
        let mut values = vec![(
            "success".to_string(),
            Pinned::N(u64::from(solution.success)),
        )];
        record_states(&mut values, &solution.t, &states);
        cases.push((
            format!("{prefix}/{gravity_name}/fast_integrate_with_gravity"),
            values,
        ));
    }
    cases
}

fn probe_kernel_direct() -> Vec<Case> {
    kernels(
        "direct(1.1,335)",
        (1.1, 335.0, 0.0, 0.0),
        &kernel_inputs(0.0, 15.0, 1013.25),
    )
}

fn probe_kernel_standard_rh0() -> Vec<Case> {
    let mut cases = Vec::new();
    for (base_alt, base_t, base_p) in BASES {
        let (density, _) = calculate_atmosphere(base_alt, Some(base_t), Some(base_p), 0.0);
        cases.extend(kernels(
            &format!("standard/{base_alt}m/{base_t}C/{base_p}hPa/rh0"),
            (base_alt, base_t, base_p, density / 1.225),
            &kernel_inputs(base_alt, base_t, base_p),
        ));
    }
    cases
}

fn probe_local_air_dry() -> Vec<Case> {
    local_air(&[0.0])
}

fn probe_local_air_humid() -> Vec<Case> {
    local_air(&[37.5, 50.0, 100.0])
}

#[test]
#[ignore = "generator: prints the pinned constants from the current engine"]
fn print_current_pin_values() {
    print("LOCAL_AIR_DRY_GRAVITY", &probe_local_air_dry());
    print("LOCAL_AIR_HUMID_GRAVITY", &probe_local_air_humid());
    print("VAPOR_MOLE_FRACTION", &probe_vapor_mole_fraction());
    print("KERNEL_DIRECT_GRAVITY", &probe_kernel_direct());
    print("KERNEL_STANDARD_RH0_GRAVITY", &probe_kernel_standard_rh0());
}

#[test]
fn local_air_dry_at_non_earth_gravity_is_byte_identical() {
    check(
        "LOCAL_AIR_DRY_GRAVITY",
        &probe_local_air_dry(),
        LOCAL_AIR_DRY_GRAVITY,
    );
}

#[test]
fn local_air_humid_at_non_earth_gravity_is_byte_identical() {
    check(
        "LOCAL_AIR_HUMID_GRAVITY",
        &probe_local_air_humid(),
        LOCAL_AIR_HUMID_GRAVITY,
    );
}

#[test]
fn vapor_mole_fraction_is_byte_identical() {
    check(
        "VAPOR_MOLE_FRACTION",
        &probe_vapor_mole_fraction(),
        VAPOR_MOLE_FRACTION,
    );
}

#[test]
fn kernel_direct_mode_at_non_earth_gravity_is_byte_identical() {
    check(
        "KERNEL_DIRECT_GRAVITY",
        &probe_kernel_direct(),
        KERNEL_DIRECT_GRAVITY,
    );
}

#[test]
fn kernel_standard_mode_at_rh0_and_non_earth_gravity_is_byte_identical() {
    check(
        "KERNEL_STANDARD_RH0_GRAVITY",
        &probe_kernel_standard_rh0(),
        KERNEL_STANDARD_RH0_GRAVITY,
    );
}

// ----- PINNED CONSTANTS (generated; see the module docs) -----------------------------------------

/// `get_local_atmosphere_humid_with_gravity` at 1.62 and 0 m/s², RH 0.
///
/// IMMUTABLE: no later step may change it, S8 included.
#[rustfmt::skip]
const LOCAL_AIR_DRY_GRAVITY: Pins = &[
    (
        "g1.62/0m/15C/1013.25hPa/rh0",
        &[
            0x3ff398c24c5f27a9, 0x407550f97c7d0dde, 0x3ff39be2780d2f46, 0x4075441783de95ea,
            0x3ff3a33e22100ec4, 0x407525ea95679bf5, 0x3ff3ca45deb5a0b1, 0x40748956c61b0f3a,
            0x3ff3d8eee23ffb31, 0x4072711d212dfca7,
        ],
    ),
    (
        "g1.62/1500m/25C/830hPa/rh0",
        &[
            0x3fef06312a50b35b, 0x4075ae6ed93deaf1, 0x3fef0af8f3959620, 0x4075a1c6117d757d,
            0x3fef1638f5b3e046, 0x407584201ecddd6a, 0x3fef51d0a6f1bac3, 0x4074ea609b63c1f3,
            0x3fee2c125fdf2cb5, 0x40734396676250bd,
        ],
    ),
    (
        "g0/0m/15C/1013.25hPa/rh0",
        &[
            0x3ff38530904bd755, 0x407550f97c7d0dde, 0x3ff39cde36735ce2, 0x4075441783de95ea,
            0x3ff3d4fe5a157d7b, 0x407525ea95679bf5, 0x3ff507e8c978c5df, 0x40748956c61b0f3a,
            0x3ffa1460243f8b92, 0x4072711d212dfca7,
        ],
    ),
    (
        "g0/1500m/25C/830hPa/rh0",
        &[
            0x3feee842a053af88, 0x4075ae6ed93deaf1, 0x3fef0c79ef6e2ca4, 0x4075a1c6117d757d,
            0x3fef62462bf48442, 0x407584201ecddd6a, 0x3ff09b26305c61af, 0x4074ea609b63c1f3,
            0x3ff3934e005d5088, 0x40734396676250bd,
        ],
    ),
];

/// As `LOCAL_AIR_DRY_GRAVITY`, at RH 37.5, 50 and 100.
///
/// MAY CHANGE ONLY IN S8 (the unified aloft vapour model, plan C7): byte-identical through every
/// other step.
///
/// Re-pinned in S8 (0.50.0, plan C7, owner-approved): the unified aloft vapour model holds the
/// governing air's (the station's or the zone's) water-vapour mole fraction for the lapse, the
/// density and the sound speed, with a virtual-temperature lapse. 90 of 120 values moved, all in
/// cases with a humidity above 0 and all off the station's own height; every RH-0 and
/// station-height value is unchanged. The largest float shift of each moved case, pinned before S8
/// -> pinned now:
/// - `g1.62/0m/15C/1013.25hPa/rh37.5` +12000m.speed_of_sound: 295.07029682688534 -> 295.3647829860783 (+9.98e-04)
/// - `g1.62/0m/15C/1013.25hPa/rh50` +12000m.speed_of_sound: 295.07052519493107 -> 295.4630462416408 (+1.33e-03)
/// - `g1.62/0m/15C/1013.25hPa/rh100` +12000m.speed_of_sound: 295.0714386636874 -> 295.85546512189603 (+2.66e-03)
/// - `g1.62/1500m/25C/830hPa/rh37.5` +12000m.speed_of_sound: 308.23181071906737 -> 308.9226251975872 (+2.24e-03)
/// - `g1.62/1500m/25C/830hPa/rh50` +12000m.speed_of_sound: 308.2343409820569 -> 309.1547458201768 (+2.99e-03)
/// - `g1.62/1500m/25C/830hPa/rh100` +12000m.speed_of_sound: 308.2444616313066 -> 310.0798426319318 (+5.95e-03)
/// - `g0/0m/15C/1013.25hPa/rh37.5` +12000m.speed_of_sound: 295.0701335029799 -> 295.3647829860783 (+9.99e-04)
/// - `g0/0m/15C/1013.25hPa/rh50` +12000m.speed_of_sound: 295.07030743001155 -> 295.4630462416408 (+1.33e-03)
/// - `g0/0m/15C/1013.25hPa/rh100` +12000m.speed_of_sound: 295.0710031361508 -> 295.85546512189603 (+2.66e-03)
/// - `g0/1500m/25C/830hPa/rh37.5` +12000m.speed_of_sound: 308.2300733449599 -> 308.9226251975872 (+2.25e-03)
/// - `g0/1500m/25C/830hPa/rh50` +12000m.speed_of_sound: 308.2320245158974 -> 309.1547458201768 (+2.99e-03)
/// - `g0/1500m/25C/830hPa/rh100` +12000m.speed_of_sound: 308.23982896017765 -> 310.0798426319318 (+5.97e-03)
#[rustfmt::skip]
const LOCAL_AIR_HUMID_GRAVITY: Pins = &[
    (
        "g1.62/0m/15C/1013.25hPa/rh37.5",
        &[
            0x3ff38cc3cc22bfbf, 0x4075566ef3fa83bd, 0x3ff38feea94a838b, 0x40754989aec9bf48,
            0x3ff39763879643f5, 0x40752b55060f3fc5, 0x3ff3bef1b23fbc67, 0x40748e991dc52527,
            0x3ff3d02681454e9a, 0x407275d626af357f,
        ],
    ),
    (
        "g1.62/0m/15C/1013.25hPa/rh50",
        &[
            0x3ff388c5a001bc86, 0x407558402c591d60, 0x3ff38bf40a709a62, 0x40754b59cdff6554,
            0x3ff3937147fa67ce, 0x40752d2292b33d9c, 0x3ff3bb2c1b3305ae, 0x40749059512b6a46,
            0x3ff3cd3999d541bc, 0x40727768a32d0621,
        ],
    ),
    (
        "g1.62/0m/15C/1013.25hPa/rh100",
        &[
            0x3ff378d38868b1d0, 0x40755f820d3cdb32, 0x3ff37c1017eb7a73, 0x407552974c0fd5aa,
            0x3ff383aeaceb0569, 0x40753455cabd120d, 0x3ff3ac1b60e634a3, 0x407497573a4bf840,
            0x3ff3c188ba4f8910, 0x40727daffc321698,
        ],
    ),
    (
        "g1.62/1500m/25C/830hPa/rh37.5",
        &[
            0x3feedb1acdca4501, 0x4075bb027ddc6474, 0x3feee0078b186721, 0x4075ae525e3a4c80,
            0x3feeeb9e9f7f0900, 0x4075909b38d608ec, 0x3fef29060a8f463c, 0x4074f682861ff1bc,
            0x3fee0d34f7bba4f5, 0x40734ec312a3a1a4,
        ],
    ),
    (
        "g1.62/1500m/25C/830hPa/rh50",
        &[
            0x3feeccc549ead76a, 0x4075bf3090cb53ad, 0x3feed1be3ebd9fe1, 0x4075b27e005e177b,
            0x3feedd721ccf4696, 0x407594c123b39ea0, 0x3fef1b72e35edc72, 0x4074fa8accd3ce53,
            0x3fee02ec4cf5da55, 0x40735279d6c0cda2,
        ],
    ),
    (
        "g1.62/1500m/25C/830hPa/rh100",
        &[
            0x3fee93930f829b77, 0x4075cfd940e49d0e, 0x3fee98bc61b65928, 0x4075c31cf667c74b,
            0x3feea4e2395588ea, 0x4075a54951fc4f79, 0x3feee5425bb1beb8, 0x40750a9cd923a6e1,
            0x3fedd9cede71f4f5, 0x4073614709114f94,
        ],
    ),
    (
        "g0/0m/15C/1013.25hPa/rh37.5",
        &[
            0x3ff37949fd914cde, 0x4075566ef3fa83bd, 0x3ff390e9343182ac, 0x40754989aec9bf48,
            0x3ff3c8e7202f926c, 0x40752b55060f3fc5, 0x3ff4fb167288e05f, 0x40748e991dc52527,
            0x3ffa0479d72d08ad, 0x407275d626af357f,
        ],
    ),
    (
        "g0/0m/15C/1013.25hPa/rh50",
        &[
            0x3ff37553c6b349b6, 0x407558402c591d60, 0x3ff38cee2f128acb, 0x40754b59cdff6554,
            0x3ff3c4e0b6f38904, 0x40752d2292b33d9c, 0x3ff4f6d1bee29982, 0x40749059512b6a46,
            0x3ff9ff2ed2e39b1f, 0x40727768a32d0621,
        ],
    ),
    (
        "g0/0m/15C/1013.25hPa/rh100",
        &[
            0x3ff365816d5cfda6, 0x40755f820d3cdb32, 0x3ff37d08a49d7561, 0x407552974c0fd5aa,
            0x3ff3b4cdaec0b68f, 0x40753455cabd120d, 0x3ff4e5c5f35b4836, 0x407497573a4bf840,
            0x3ff9ea0b73bcbfa7, 0x40727daffc321698,
        ],
    ),
    (
        "g0/1500m/25C/830hPa/rh37.5",
        &[
            0x3feebd7f2a768141, 0x4075bb027ddc6474, 0x3feee1845db76168, 0x4075ae525e3a4c80,
            0x3fef36d9e45d739f, 0x4075909b38d608ec, 0x3ff0842c44aa6505, 0x4074f682861ff1bc,
            0x3ff3783852f6c8ba, 0x40734ec312a3a1a4,
        ],
    ),
    (
        "g0/1500m/25C/830hPa/rh50",
        &[
            0x3feeaf45283e1484, 0x4075bf3090cb53ad, 0x3feed339afeafd11, 0x4075b27e005e177b,
            0x3fef2867b87114f5, 0x407594c123b39ea0, 0x3ff07c8777800b33, 0x4074fa8accd3ce53,
            0x3ff36f359f9efcce, 0x40735279d6c0cda2,
        ],
    ),
    (
        "g0/1500m/25C/830hPa/rh100",
        &[
            0x3fee7680529c3a2a, 0x4075cfd940e49d0e, 0x3fee9a325538553c, 0x4075c31cf667c74b,
            0x3feeeec2c6f5847f, 0x4075a54951fc4f79, 0x3ff05e072c759dff, 0x40750a9cd923a6e1,
            0x3ff34b411d9beb3f, 0x4073614709114f94,
        ],
    ),
];

/// The uncapped CIPM-2007 vapour mole fraction (the density-altitude refusal's predicate).
///
/// IMMUTABLE: no later step may change it, S8 included.
#[rustfmt::skip]
const VAPOR_MOLE_FRACTION: Pins = &[
    (
        "water_vapor_mole_fraction_uncapped",
        &[
            0x0000000000000000, 0x3f5170ffe77f3ca4, 0x3f57415534a9a631, 0x3f67415534a9a631,
            0x0000000000000000, 0x3f79f54988d905be, 0x3f814e3105e603d5, 0x3f914e3105e603d5,
            0x0000000000000000, 0x3f9987fb5cd42c90, 0x3fa105523de2c861, 0x3fb105523de2c861,
            0x0000000000000000, 0x3fe47408783481ea, 0x3feb4560a046028e, 0x3ffb4560a046028e,
        ],
    ),
];

/// Kernel direct (uniform) mode at 1.62 and 0 m/s²: `compute_derivatives_with_gravity`,
/// `integrate_trajectory_with_gravity` (RK45) and `fast_integrate_with_gravity`.
///
/// IMMUTABLE: no later step may change it, S8 included.
#[rustfmt::skip]
const KERNEL_DIRECT_GRAVITY: Pins = &[
    (
        "direct(1.1,335)/g1.62/compute_derivatives_with_gravity",
        &[
            0x4085e00000000000, 0x4000000000000000, 0x3fe0000000000000, 0xc076ac6e3f7438fc,
            0xc0054085f12bf796, 0xc002a8379c1390e7, 0x4085e00000000000, 0x4000000000000000,
            0x3fe0000000000000, 0xc076ac6e3f7438fc, 0xc0054085f12bf796, 0xc002a8379c1390e7,
        ],
    ),
    (
        "direct(1.1,335)/g1.62/integrate_trajectory_with_gravity/RK45",
        &[
            0x000000000000002a, 0x3f9db22d0e56041a, 0x4037ac2ea77075ef, 0x3fa7e4d5028bfec6,
            0xbf4ef476ae978e42, 0x40894d312b41f181, 0x3ff929c9409a9896, 0xbfb09c1355ef35e4,
            0x3fd27ef9db22d0e7, 0x406b7c4bc55d3baf, 0x3fd80882cb2bd362, 0xbfb6596675a24980,
            0x4085fc2dbe7e961c, 0x3fef24e03155a578, 0xbfe294fc0cc8f31d, 0x3fe3ced916872b05,
            0x407b1cbcf8107343, 0x3fe2cfa471a24569, 0xbfd78742a6210a2a, 0x4082abdd8bfa9c01,
            0x3fd54b02351468a6, 0xbff1897cd72ebec5, 0x3fee0c49ba5e3545, 0x40831c3ceb211fd1,
            0x3fe372d6cfdd1ce5, 0xbfe914db86fb5f9f, 0x408021eebf189619, 0xbfc8fd3de38a2c93,
            0xbff7dad21f9b0645, 0x3ff56170aaee56a1, 0x4089000000000000, 0x3fda9875ad2796d3,
            0xbff74f9ae937431f, 0x407b3c9091bce0fb, 0xbfe83ceda5dd69ea, 0xbffe1c851df775b6,
        ],
    ),
    (
        "direct(1.1,335)/g1.62/fast_integrate_with_gravity",
        &[
            0x0000000000000001, 0x000000000000053a, 0x3f50624dd2f1a9fc, 0x3fea54181aca7a60,
            0x3f5af2780841326d, 0xbeb2fed8d436f095, 0x4089b44295d3a864, 0x3ffa4b8fca7ad4bf,
            0xbf628c33edc16f93, 0x3fd5604189374bc7, 0x406f65aa487e006f, 0x3fdab31d11c631a6,
            0xbfbd7e6759d7bf37, 0x40857af10494d003, 0x3fec1f7811286601, 0xbfe518301c2ba87c,
            0x3fe56872b020c49c, 0x407cf523cba676cf, 0x3fe345dda03603e7, 0xbfdb24f3ec5dcd26,
            0x40823cd28471327e, 0x3fcf596e23a66931, 0xbff29dd0539fd8ba, 0x3ff00c49ba5e353f,
            0x408420dd2d947411, 0x3fe2f2fc69f67cdb, 0xbfec34e4aaee5a1f, 0x407f5ecb160dc9d0,
            0xbfd2b0fceb96ed8e, 0xbff8f7d6b1b747da, 0x3ff5616c72ccde9a, 0x4089000000000000,
            0x3fda98b829824c28, 0xbff74f89e967cdbd, 0x407b3c8f3030b45f, 0xbfe83cec6b73f526,
            0xbffe1c86ad725447,
        ],
    ),
    (
        "direct(1.1,335)/g0/compute_derivatives_with_gravity",
        &[
            0x4085e00000000000, 0x4000000000000000, 0x3fe0000000000000, 0xc076ac6e3f7438fc,
            0xbff09586c39f9d3f, 0xc002a8379c1390e7, 0x4085e00000000000, 0x4000000000000000,
            0x3fe0000000000000, 0xc076ac6e3f7438fc, 0xbff09586c39f9d3f, 0xc002a8379c1390e7,
        ],
    ),
    (
        "direct(1.1,335)/g0/integrate_trajectory_with_gravity/RK45",
        &[
            0x000000000000002a, 0x3f9db22d0e56041a, 0x4037ac2ea7604870, 0x3fa83da2881d2a56,
            0xbf4ef476b8a80b8a, 0x40894d312b0d9863, 0x3ff9e8a74cd73468, 0xbfb09c135e1378aa,
            0x3fd27ef9db22d0e7, 0x406b7c4bbe8f012c, 0x3fdc252d34a255e1, 0xbfb65966b9602a78,
            0x4085fc2dade07b11, 0x3ff683432cb93e84, 0xbfe294fc5f7d4f83, 0x3fe3ced916872b05,
            0x407b1cbcdd0e9c57, 0x3febc3532ea46315, 0xbfd787432c8ab1b9, 0x4082abdd555aa692,
            0x3ff31e967dfa73d2, 0xbff1897d5f1d4bb5, 0x3fee0c49ba5e3545, 0x40831c3cc7dcb4fc,
            0x3ff391a865bcce46, 0xbfe914dc36811f03, 0x408021ee7074c07c, 0x3ff0850e2e6c7164,
            0xbff7dad2e34c5490, 0x3ff56170fb5d0ac3, 0x4089000000000000, 0x3ff9999bd6436285,
            0xbff74f9c2af21557, 0x407b3c8fb1f77b5d, 0x3febe3e98c887c90, 0xbffe1c863463ba4d,
        ],
    ),
    (
        "direct(1.1,335)/g0/fast_integrate_with_gravity",
        &[
            0x0000000000000001, 0x000000000000053a, 0x3f50624dd2f1a9fc, 0x3fea54181aca7507,
            0x3f5af5dd99581e23, 0xbeb2fed8d46c38c5, 0x4089b44295d398b4, 0x3ffa523200d6cafc,
            0xbf628c33ee0f81ff, 0x3fd5604189374bc7, 0x406f65aa3e43fe34, 0x3fe0134a3bb41f8a,
            0xbfbd7e67bfa301f0, 0x40857af0ef1a5c78, 0x3ff5feec59b60e5a, 0xbfe51830871107e6,
            0x3fe56872b020c49c, 0x407cf523aaee7627, 0x3feda71096392fc6, 0xbfdb24f48f3486d1,
            0x40823cd248f25d71, 0x3ff2ace127d96501, 0xbff29dd0e7ad3109, 0x3ff00c49ba5e353f,
            0x408420dd052f4823, 0x3ff49c8a0388430d, 0xbfec34e573f99309, 0x407f5eca7326544f,
            0x3ff00fc537c0b314, 0xbff8f7d77c680391, 0x3ff5616cc350b3b9, 0x4089000000000000,
            0x3ff9999bd643628b, 0xbff74f8b2b772247, 0x407b3c8e5070f58d, 0x3febe3e82285c79b,
            0xbffe1c87c3d78f43,
        ],
    ),
];

/// Kernel standard mode at RH 0 and 1.62 and 0 m/s², the same families, at both station bases.
///
/// IMMUTABLE: no later step may change it, S8 included.
///
/// Compared to a relative 1e-15 until 0.50.0's R5, which made it bits like every other group: no
/// last-ulp movement ever happened.
#[rustfmt::skip]
const KERNEL_STANDARD_RH0_GRAVITY: Pins = &[
    (
        "standard/0m/15C/1013.25hPa/rh0/g1.62/compute_derivatives_with_gravity",
        &[
            0x4085e00000000000, 0x4000000000000000, 0x3fe0000000000000, 0xc0796c4b1919ea8e,
            0xc00641ef9426bc2e, 0xc004eb654ac7cb3f, 0x4085e00000000000, 0x4000000000000000,
            0x3fe0000000000000, 0xc079693185c7f941, 0xc00640cd5b593983, 0xc004e8d84af9653d,
        ],
    ),
    (
        "standard/0m/15C/1013.25hPa/rh0/g1.62/integrate_trajectory_with_gravity/RK45",
        &[
            0x000000000000002b, 0x3f9db22d0e56041a, 0x4037a624f7eef5dd, 0x3fa7deb5469e5084,
            0xbf515b077f2c77d3, 0x4089404c80002c6d, 0x3ff91cc5b604872c, 0xbfb29d6cc4502be6,
            0x3fd27ef9db22d0e7, 0x406b3c2be5d7d254, 0x3fd7cd0b248860ec, 0xbfb8d7af9170b29e,
            0x40859697b8ad0720, 0x3fee731a01b3a9e1, 0xbfe48e920cfd9132, 0x3fe3ced916872b05,
            0x407aa255a3cde66b, 0x3fe26b7e28d5509e, 0xbfd9e8743402b335, 0x408205b7a6686ceb,
            0x3fd37d50eb1d4676, 0xbff326f0d3ff0554, 0x3ff024dd2f1a9fc1, 0x4083b1c1360fbc48,
            0x3fe22b76dc51e265, 0xbfef227bee6a9df9, 0x407da8fbe92bcbad, 0xbfd4c020c3b4b6f3,
            0xbffb18933638cd81, 0x3ff6684e20ff8958, 0x4089000000000000, 0x3fd56d6806d35052,
            0xbffb6b10c6ff4cf7, 0x4078c580025b63d9, 0xbfeb42aff270f22a, 0xc00096dae971068c,
        ],
    ),
    (
        "standard/0m/15C/1013.25hPa/rh0/g1.62/fast_integrate_with_gravity",
        &[
            0x0000000000000001, 0x000000000000057a, 0x3f50624dd2f1a9fc, 0x3fea53dc3975da61,
            0x3f5af23abc261c95, 0xbeb552e2a54b17cb, 0x4089b3cdab8ccf8e, 0x3ffa4b1820eb3e5d,
            0xbf64d214a1c6faeb, 0x3fd6666666666667, 0x407034b17f39ea17, 0x3fdb42d1a8f5af6a,
            0xbfc1e3bd591a3d26, 0x4084d9964fefb999, 0x3fea4cdb0c5d86e8, 0xbfe83b3c6e8220c3,
            0x3fe66e978d4fdf3c, 0x407d871e4f32fb4c, 0x3fe307c5a2942d5d, 0xbfe03fca51e42737,
            0x4081463821ce10e9, 0x3fc4b63b19df9256, 0xbff50379e1ad95e5, 0x3ff0d0e560418937,
            0x40844f9dcbecb058, 0x3fe1b16540795d42, 0xbff0b88979aef9cb, 0x407d122b6bd235a2,
            0xbfd8a64a8a8e2a19, 0xbffbd438f245cd44, 0x3ff6684b4fc64fbd, 0x4089000000000000,
            0x3fd56d919dc39f1c, 0xbffb6b055d4d928b, 0x4078c57ee39c17eb, 0xbfeb42aeb650a6b6,
            0xc00096db8590ced8,
        ],
    ),
    (
        "standard/0m/15C/1013.25hPa/rh0/g0/compute_derivatives_with_gravity",
        &[
            0x4085e00000000000, 0x4000000000000000, 0x3fe0000000000000, 0xc0796c4b1919ea8e,
            0xbff2985a09952671, 0xc004eb654ac7cb3f, 0x4085e00000000000, 0x4000000000000000,
            0x3fe0000000000000, 0xc0798fa831ca3abe, 0xbff2b237bb1792d5, 0xc005087eb27a852f,
        ],
    ),
    (
        "standard/0m/15C/1013.25hPa/rh0/g0/integrate_trajectory_with_gravity/RK45",
        &[
            0x000000000000002b, 0x3f9db22d0e56041a, 0x4037a624f6bd8653, 0x3fa83773be37f880,
            0xbf515b07de2e8d83, 0x4089404c7c2d231e, 0x3ff9db736569c357, 0xbfb29d6d5c94e351,
            0x3fd27ef9db22d0e7, 0x406b3c2b6ff4eddb, 0x3fdbe382e280ff65, 0xbfb8d7b426dccfa3,
            0x408596969f4ef59a, 0x3ff61b3be9b51732, 0xbfe48e978555bfd7, 0x3fe3ced916872b05,
            0x407aa253e892876c, 0x3feb45fa18f3e651, 0xbfd9e87cd1f1311f, 0x408205b42337d2c6,
            0x3ff274705871738c, 0xbff326f9919892c4, 0x3ff024dd2f1a9fc1, 0x4083b1be6a9e25f4,
            0x3ff42ac0a716c31e, 0xbfef2289d7201dc8, 0x407da8ef50052561, 0x3fee5f2cf53ffe89,
            0xbffb18a2e30e05ad, 0x3ff6685582617e79, 0x4089000000000000, 0x3ff9999bd6436287,
            0xbffb6b2e4c874923, 0x4078c56a1fac221a, 0x3fe95d9dfdbdda1f, 0xc00096e886f55fc9,
        ],
    ),
    (
        "standard/0m/15C/1013.25hPa/rh0/g0/fast_integrate_with_gravity",
        &[
            0x0000000000000001, 0x000000000000057a, 0x3f50624dd2f1a9fc, 0x3fea53dc39757137,
            0x3f5af5a048165417, 0xbeb552e2a961f559, 0x4089b3cdab8b9b86, 0x3ffa51ba4830b1dd,
            0xbf64d214a7c405e9, 0x3fd6666666666667, 0x407034b11bae42ab, 0x3fe098441eec9eb5,
            0xbfc1e3c137f6a3f6, 0x4084d994ce3f2020, 0x3ff559b0c546a4e6, 0xbfe83b43ee0ea5f4,
            0x3fe66e978d4fdf3c, 0x407d871bf367df54, 0x3fee3c89c27d7e93, 0xbfe03fd0306875b7,
            0x40814633ff5274ad, 0x3ff1b0578f1830a8, 0xbff503842bac18ee, 0x3ff0d0e560418937,
            0x40844f9abb52ba25, 0x3ff4cc66eb2f339f, 0xbff0b8911a22592b, 0x407d121e4b4f06a2,
            0x3fedc4bd453c11a3, 0xbffbd4494786b171, 0x3ff66852ae3dd436, 0x4089000000000000,
            0x3ff9999bd6436290, 0xbffb6b22d72ba419, 0x4078c568ffc66a16, 0x3fe95d9cd6ef30e1,
            0xc00096e923cc3fea,
        ],
    ),
    (
        "standard/1500m/25C/830hPa/rh0/g1.62/compute_derivatives_with_gravity",
        &[
            0x4085e00000000000, 0x4000000000000000, 0x3fe0000000000000, 0xc07442c83e1e18ac,
            0xc0045ea3edc1dd24, 0xc000abfb1464d566, 0x4085e00000000000, 0x4000000000000000,
            0x3fe0000000000000, 0xc074408e2e02815b, 0xc0045dd372d67e4e, 0xc000aa25ffd34005,
        ],
    ),
    (
        "standard/1500m/25C/830hPa/rh0/g1.62/integrate_trajectory_with_gravity/RK45",
        &[
            0x000000000000002a, 0x3f9db22d0e56041a, 0x4037b156bee4edf2, 0x3fa7ea0feb429aaa,
            0xbf4bbf36d12bbba3, 0x408958393fa605e0, 0x3ff934ebb8d7ce10, 0xbfadc9b2fda9bf18,
            0x3fd27ef9db22d0e7, 0x406bb450d84ab8eb, 0x3fd83c6ff0abd722, 0xbfb42bc9964793d1,
            0x4086561afa91ac66, 0x3fefc2311a10d127, 0xbfe0d56d1cdbb76f, 0x3fe3ced916872b05,
            0x407b8a47c1402897, 0x3fe3291dba3b4f3d, 0xbfd56613d5216808, 0x4083431dc777e3fd,
            0x3fd6f0908e63c260, 0xbff0111b0316066d, 0x3fed6872b020c4a1, 0x408332a9e81e0f1a,
            0x3fe42782a15782b8, 0xbfe615e2e4290b27, 0x4080fe29a327bc4f, 0xbfc20748f7c90456,
            0xbff5b6c964c5c763, 0x3ff49369f2405fe8, 0x4089000000000000, 0x3fde9292b001f169,
            0xbff41780022d4fde, 0x407d6959a0621c1e, 0xbfe5a0804d7279bf, 0xbffb67c016c209c8,
        ],
    ),
    (
        "standard/1500m/25C/830hPa/rh0/g1.62/fast_integrate_with_gravity",
        &[
            0x0000000000000001, 0x0000000000000507, 0x3f50624dd2f1a9fc, 0x3fea544b088bde43,
            0x3f5af2ac2a8dd124, 0xbeb103e9005ac72b, 0x4089b4a6081d10c1, 0x3ffa4bf592f84b22,
            0xbf609d43fa7c22cd, 0x3fd48b4395810625, 0x406e8ac90037b8cc, 0x3fda335fceef6ffc,
            0xbfb8b13adbd55db1, 0x40860150d7f7b661, 0x3feda41d1548e23c, 0xbfe27b6acd7a4c72,
            0x3fe49374bc6a7efa, 0x407c75c1d702b1cc, 0x3fe36b6a6ad19ec3, 0xbfd6f6fd915a94bf,
            0x4083108d69eb67e7, 0x3fd43a9b19b29875, 0xbff08eee7cbc4532, 0x3feee147ae147ae1,
            0x4083f8eb0e904508, 0x3fe3e4f8cf737a1d, 0xbfe81e8430337438, 0x4080adf3bfe78a14,
            0xbfcb258145a8ca4e, 0xbff67e631367f643, 0x3ff49365badd684d, 0x4089000000000000,
            0x3fde92da2927bca7, 0xbff4176f09a9f484, 0x407d695863b4be9e, 0xbfe5a07f655416b6,
            0xbffb67c17bfeab2b,
        ],
    ),
    (
        "standard/1500m/25C/830hPa/rh0/g0/compute_derivatives_with_gravity",
        &[
            0x4085e00000000000, 0x4000000000000000, 0x3fe0000000000000, 0xc07442c83e1e18ac,
            0xbfeda3857996d0b6, 0xc000abfb1464d566, 0x4085e00000000000, 0x4000000000000000,
            0x3fe0000000000000, 0xc0745e296bad5eea, 0xbfedcb92e6a8c55d, 0xc000c282a1beef04,
        ],
    ),
    (
        "standard/1500m/25C/830hPa/rh0/g0/integrate_trajectory_with_gravity/RK45",
        &[
            0x000000000000002a, 0x3f9db22d0e56041a, 0x4037b156bdf7ab86, 0x3fa842ea4dc9f413,
            0xbf4bbf3764c60cd2, 0x408958393cac9289, 0x3ff9f3f32687133a, 0xbfadc9b3ea84a6f7,
            0x3fd27ef9db22d0e7, 0x406bb45078d53e9b, 0x3fdc5e8a2159fe6f, 0xbfb42bcd4c768562,
            0x4086561a132adde3, 0x3ff6df58175f12aa, 0xbfe0d5719c86c2a4, 0x3fe3ced916872b05,
            0x407b8a4648fa1d90, 0x3fec337da2db7be4, 0xbfd5661b25d0cf6a, 0x4083431ab3469432,
            0x3ff3b97520ef4cc0, 0xbff01122ac7934c6, 0x3fed6872b020c4a1, 0x408332a7ec8e5229,
            0x3ff3a89d49bbe4a8, 0xbfe615ecc2415907, 0x4080fe2468de71ec, 0x3ff1668d3436ec14,
            0xbff5b6d6670c15fe, 0x3ff4936eac099566, 0x4089000000000000, 0x3ff9999bd6436288,
            0xbff41792e9523f41, 0x407d694723d2db81, 0x3fee1dfda7ac426b, 0xbffb67d71727c831,
        ],
    ),
    (
        "standard/1500m/25C/830hPa/rh0/g0/fast_integrate_with_gravity",
        &[
            0x0000000000000001, 0x0000000000000507, 0x3f50624dd2f1a9fc, 0x3fea544b088b8ce8,
            0x3f5af611c005baf1, 0xbeb103e903849bcb, 0x4089b4a6081c2275, 0x3ffa5297d627ec18,
            0xbf609d43ff1e1faf, 0x3fd48b4395810625, 0x406e8ac8803aacd5, 0x3fdf4671d67902c9,
            0xbfb8b13fd5cfa3e1, 0x4086014fc2a2f28a, 0x3ff68884cba822e7, 0xbfe27b7031bbc793,
            0x3fe49374bc6a7efa, 0x407c75c037d72b3c, 0x3fed249e6afcd7ed, 0xbfd6f705a39f93d5,
            0x4083108a29688b10, 0x3ff385ade74c54f3, 0xbff08ef6946861c9, 0x3feee147ae147ae1,
            0x4083f8e8d391ab59, 0x3ff473a054519d36, 0xbfe81e8f49ff48e7, 0x4080adee3462e2ee,
            0x3ff1146a2666296c, 0xbff67e70dfd2aaa2, 0x3ff4936a7551f7c1, 0x4089000000000000,
            0x3ff9999bd6436291, 0xbff41781f37c3253, 0x407d6945e75de3da, 0x3fee1dfc639eded7,
            0xbffb67d87c1e1ec4,
        ],
    ),
];

#[test]
fn ulps_apart_counts_representable_values() {
    assert_eq!(ulps_apart(1.0, 1.0), 0);
    assert_eq!(ulps_apart(0.0, -0.0), 0);
    assert_eq!(ulps_apart(1.0, f64::from_bits(1.0f64.to_bits() + 2)), 2);
    assert_eq!(ulps_apart(-1.0, f64::from_bits((-1.0f64).to_bits() + 3)), 3);
    assert_eq!(ulps_apart(5e-324, -5e-324), 2);
}
