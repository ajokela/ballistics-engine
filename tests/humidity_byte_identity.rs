//! Byte-identity pins for the 0.50.0 humidity work, captured on 0.49.0 before any of it.
//!
//! # Why this file exists
//!
//! 0.50.0 routes every solver family through one air lookup and gives the kernel API a type
//! that cannot describe an atmosphere without its relative humidity. Its promise is that every
//! `TrajectorySolver`-backed surface stays BYTE-IDENTICAL at every valid RH, and so do direct
//! ("uniform") mode on every family, station density and sound speed, the muzzle Miller Sg,
//! Litz spin drift, aerodynamic jump and density altitude. These pins hold 0.50.0 to that
//! promise with `f64::to_bits` equality. A refactor that moves one of them has changed physics,
//! whatever its description says.
//!
//! # Rules
//!
//! - A pin is never edited. It may only be MIGRATED: when an API it calls is renamed or
//!   retyped (for example the kernel atmosphere tuple becoming `KernelAtmosphere`), the call
//!   changes and the pinned constant does not.
//! - Each group's header says what in 0.50.0 may change it. Two words are used:
//!   - `IMMUTABLE`: must pass unchanged through all of 0.50.0, the aloft vapour model included.
//!   - `MAY CHANGE ONLY WITH THE ALOFT VAPOUR MODEL`: the unified aloft vapour model may re-pin
//!     it, and no other change may. Every such group involves a relative humidity above zero,
//!     at the station or in a zone; at RH 0 the aloft model is bit-identical, so the RH-0 cases
//!     live in separate `IMMUTABLE` groups.
//!   - The aloft vapour model re-pinned every such group that moved (`ZERO_HUMID` did not: its
//!     trial heights kept their signs), and each re-pinned group carries a note with its count
//!     of moved values and, for every moved case, its largest shift, old -> new. No RH-0 value
//!     and no value at the station's own height moved, and every `IMMUTABLE` group passed it
//!     unchanged.
//! - Kernel outputs at RH > 0 (`compute_derivatives`, `try_integrate_trajectory`,
//!   `try_solve_trajectory_rust`, `fast_integrate*`, the lean Monte Carlo) are deliberately NOT
//!   pinned: they flew the dry speed of sound at any RH until 0.50.0 fixed that.
//! - Kernel STANDARD mode at RH 0 is pinned bit for bit like everything else. It was first
//!   allowed a relative 1e-15 (0.50.0 moved the kernels onto `air_at`, and fast_integrate lost
//!   a `ratio * 1.225 / 1.225` round trip), but the bits held through every change on every
//!   platform measured, so 0.50.0 dropped the tolerance.
//!
//! In-crate pins for the `pub(crate)` gravity-aware entry points live in
//! `src/humidity_byte_identity_tests.rs` under the same rules.
//!
//! # How each group is checked, and on which platforms
//!
//! A libm may round a `pow`, `exp` or `log` a last bit apart from another, and one such bit can
//! reach a pinned value. Bits captured on one machine are therefore not bits on every machine:
//! the BSDs' libms moved some of these pins even where 0.50.0 computes exactly what 0.49.0 does
//! (FreeBSD's even on the unchanged 0.49.0 code), and musl's moved two of the aloft model's
//! groups. So each group
//! is held in one of two ways (`Check`):
//!
//! - **`Check::Released049`: against the released 0.49.0, on every platform.** A group whose
//!   promise is "0.50.0 equals 0.49.0", and whose computation 0.49.0 exposes publicly, runs
//!   twice: on this crate, and on the released 0.49.0 crate from crates.io (the
//!   `ballistics_engine_049` dev-dependency, `ffi` off so the two C ABIs cannot collide). Both
//!   runs use the same probe source, `humidity_byte_identity/probes.rs`, `include!`d once per
//!   engine. They must agree bit for bit in every value, digests included, on whatever platform
//!   runs the test, so these groups need no table of libm alternates. On an exact-bit platform
//!   (below) they must also match their literal pins bit for bit.
//! - **`Check::Pinned`: against the literal pins only.** These are the groups the aloft vapour
//!   model moved (no engine computes their new values but this one), plus `ZERO_HUMID`, which
//!   it was allowed to move, and the IMMUTABLE groups 0.49.0 has no public entry for:
//!   `FFI_TRAJECTORY_DRY` (0.49.0's C ABI cannot be linked beside this crate's). On an exact-bit
//!   platform every value must match its pin bit for bit, or, on glibc, the glibc bits recorded in
//!   `LIBM_ALTERNATES`. On any other platform each float must lie within `UNVERIFIED_LIBM_ULPS` of
//!   its pin, every count, flag and text hash must match exactly, and digests are not compared (a
//!   digest cannot be compared within a tolerance; its trajectory's explicit indices still are).
//!
//! The exact-bit platforms are the ones the literal pins were captured or checked on (`Libm`):
//! macOS aarch64 (Apple's libm, where they were captured) and Linux glibc on aarch64 and x86_64
//! (the CI platform). Every other platform is unverified, macOS x86_64 included: Apple's x86_64
//! libm rounds some values apart from its aarch64 one. The bound there, `UNVERIFIED_LIBM_ULPS`,
//! is 2 ulps: the largest distance any `Check::Pinned` float showed when 0.50.0 ran these
//! pins on FreeBSD 14.3, NetBSD 10.1 and OpenBSD 7.8 (amd64), Alpine Linux musl (aarch64 and
//! x86_64), macOS x86_64 (Rosetta 2) and Windows 11 (x86_64 MSVC, the UCRT's own x64 code). No
//! count, flag or text hash moved on any of them. Its doc has the details.
//!
//! # How the pinned values were generated
//!
//! By running this file's own probes against the engine at 261c1eb (ballistics-engine
//! 0.49.0, the commit the humidity work started from) and pasting the printed constants verbatim:
//!
//! ```text
//! cargo test --test humidity_byte_identity -- --ignored --nocapture print_current_pin_values
//! ```
//!
//! They were captured on macOS (aarch64), and the identical file then passed in Linux glibc
//! containers on aarch64 and x86_64 (the CI platform). The aloft vapour model re-pinned its
//! groups the same way, and glibc's bits for the seven values where it differed are in
//! `LIBM_ALTERNATES`.
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

// =================================================================================================
// Pin machinery
// =================================================================================================

/// One pinned quantity: the bits of an `f64`, an exact integer (a count, a flag, a text hash), or
/// a digest over a whole trajectory.
#[derive(Clone, Copy)]
enum Pin {
    F(f64),
    N(u64),
    /// A digest: compared bit for bit wherever bits are compared, and not at all where floats
    /// are compared within ulps (one ulp anywhere changes every bit of it).
    D(u64),
}

impl Pin {
    fn bits(self) -> u64 {
        match self {
            Pin::F(value) => value.to_bits(),
            Pin::N(value) | Pin::D(value) => value,
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
        self.entries
            .push((format!("{label} digest"), Pin::D(digest_f64s(values))));
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
        Pin::N(value) | Pin::D(value) => format!("{value}"),
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

/// A libm the literal pins are enforced against bit for bit: one they were captured or checked
/// on, on the architectures they were checked on. Each is named by the platform measured, not by
/// a `target_env` alone: a glibc on another architecture, or a `target_env = "gnu"` that is not
/// Linux (MinGW), is unverified and falls under [`UNVERIFIED_LIBM_ULPS`].
#[derive(Clone, Copy, PartialEq, Debug)]
enum Libm {
    /// Apple's libm on macOS aarch64, where the pins were captured. Not macOS x86_64: Apple's
    /// x86_64 libm (measured under Rosetta 2) rounds some values apart from its aarch64 one,
    /// up to 1137 ulps through `density_matched_altitude_m`'s bisection, so it is unverified.
    Apple,
    /// glibc, Linux aarch64 and x86_64 (the CI platform): every value matched except the seven
    /// in [`LIBM_ALTERNATES`].
    Glibc,
}

impl Libm {
    /// The exact-bit libm this build links, if it is one.
    fn of_this_build() -> Option<Libm> {
        if cfg!(all(target_os = "macos", target_arch = "aarch64")) {
            Some(Libm::Apple)
        } else if cfg!(all(
            target_os = "linux",
            target_env = "gnu",
            any(target_arch = "aarch64", target_arch = "x86_64")
        )) {
            Some(Libm::Glibc)
        } else {
            None
        }
    }
}

/// The widest distance, in units in the last place, a float of a `Check::Pinned` group may lie
/// from its pin on a platform that is not an exact-bit one ([`Libm`]).
///
/// Measured, not chosen. In 0.50.0 this file's pins (then all compared bit for bit) were run
/// with the 0.50.0 engine on FreeBSD 14.3, NetBSD 10.1 and OpenBSD 7.8 (amd64, rustc 1.93 and
/// 1.90), Alpine Linux musl on aarch64 and x86_64, macOS x86_64 (under Rosetta 2) and Windows 11
/// (x86_64 MSVC with the UCRT linked statically, so its own x64 math code ran, under x86_64
/// emulation on an ARM64 machine). On every one of them the farthest any float of a
/// `Check::Pinned` group landed from its pin was 2 ulps: a Litz spin drift
/// (`summary.spin_drift_m`, a `powf`), lateral and vertical positions late in the RK4 and RK45
/// lob, an impact energy, and densities at +12 km and +25 km. Every count, flag and text hash
/// matched; the other differences were digests. The rewritten file then passed on all of them,
/// and on Windows with the UCRT linked dynamically too. The `Check::Released049` groups' literal
/// pins moved farther on some (11 ulps on the BSDs and 1137 ulps on macOS x86_64, both through
/// `density_matched_altitude_m`'s bisection, which amplifies a last-bit difference), which is
/// why those groups are compared to 0.49.0 there instead; 0.49.0 matched this engine bit for bit
/// on every platform. Following be60e2b, the bound sits far below a physics change: every float
/// the aloft vapour model re-pinned moved by 327 ulps or more, and its only smaller changes were
/// two trajectory point counts, which every platform compares exactly.
const UNVERIFIED_LIBM_ULPS: u64 = 2;

/// The bit patterns another libm gives a few pinned values: `(group, case, label, libm, bits)`.
///
/// The unified aloft vapour model re-pinned its groups from macOS (aarch64). In the Linux
/// glibc containers, aarch64 and x86_64 alike, every re-pinned value then matched except these
/// seven: glibc's `pow` rounds the virtual-temperature lapse a last bit apart from Apple's at some
/// state of the RK4 30-degree lob at 50 %, and that ulp reaches the lob's final point (its speed,
/// energy, Cd and height, one ulp each, and so the position digest) and one Monte Carlo trial's
/// range (the ranges digest). Each value stays pinned bit for bit on every exact-bit platform: an
/// alternate is accepted only on the libm it was measured on, and every other exact-bit platform
/// must match the primary pin. (musl happens to give these seven bits too, among others; it is
/// not an exact-bit platform, so it is held to [`UNVERIFIED_LIBM_ULPS`] instead.)
/// [`every_libm_alternate_names_a_pinned_value_it_differs_from`] keeps the table honest.
const LIBM_ALTERNATES: &[(&str, &str, &str, Libm, u64)] = &[
    (
        "MONTE_CARLO_DEFAULT_HUMIDITY",
        "seeded/icao-sea-level/default-rh",
        "ranges digest",
        Libm::Glibc,
        0xaaf6b70a54d99f0d,
    ),
    (
        "TS_LOB_HUMID",
        "lob30/rk4/icao-sea-level/rh50/no-zones",
        "impact_velocity",
        Libm::Glibc,
        0x4061977687c0d8b1,
    ),
    (
        "TS_LOB_HUMID",
        "lob30/rk4/icao-sea-level/rh50/no-zones",
        "impact_energy",
        Libm::Glibc,
        0x405c132490c3a61f,
    ),
    (
        "TS_LOB_HUMID",
        "lob30/rk4/icao-sea-level/rh50/no-zones",
        "points[29902].z",
        Libm::Glibc,
        0xc0574cddecb698cf,
    ),
    (
        "TS_LOB_HUMID",
        "lob30/rk4/icao-sea-level/rh50/no-zones",
        "points[29902].speed",
        Libm::Glibc,
        0x4061977687c0d8b1,
    ),
    (
        "TS_LOB_HUMID",
        "lob30/rk4/icao-sea-level/rh50/no-zones",
        "points[29902].drag_coefficient",
        Libm::Glibc,
        0x3fc092c72722a004,
    ),
    (
        "TS_LOB_HUMID",
        "lob30/rk4/icao-sea-level/rh50/no-zones",
        "points (time, x, y, z) digest",
        Libm::Glibc,
        0xa9a2631de8fe39ed,
    ),
];

/// The alternate bit pattern this build's libm may give `group`/`case`/`label`, if any.
fn libm_alternate(group: &str, case: &str, label: &str) -> Option<u64> {
    let libm = Libm::of_this_build()?;
    LIBM_ALTERNATES
        .iter()
        .find(|&&(g, c, l, alternate_libm, _)| {
            g == group && c == case && l == label && alternate_libm == libm
        })
        .map(|&(_, _, _, _, bits)| bits)
}

/// Every quantity of every case must equal its pin bit for bit (or, on the one libm that rounds
/// it differently, that libm's own recorded bits: [`LIBM_ALTERNATES`]).
fn assert_bits(group: &str, actual: &Probe, expected: Pinned) {
    case_shape_or_panic(group, actual, expected);
    let mut moved = Vec::new();
    let mut total = 0;
    for ((name, pins), (_, bits)) in actual.iter().zip(expected) {
        for ((label, pin), &want) in pins.entries.iter().zip(bits.iter()) {
            total += 1;
            if pin.bits() != want && libm_alternate(group, name, label) != Some(pin.bits()) {
                let pinned = match pin {
                    Pin::F(_) => format!("{:?}", f64::from_bits(want)),
                    Pin::N(_) | Pin::D(_) => format!("{want}"),
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

/// The distance between two finite `f64`s in units in the last place: how many representable
/// values lie between them (0 for equal values, 0.0 and -0.0 included).
fn ulps_apart(a: f64, b: f64) -> u64 {
    // Map each bit pattern onto a line where adjacent floats are adjacent integers.
    let line = |x: f64| -> i128 {
        let bits = x.to_bits() as i64;
        i128::from(if bits < 0 { i64::MIN - bits } else { bits })
    };
    (line(a) - line(b)).unsigned_abs() as u64
}

/// On a platform whose libm the pins were not captured or checked on: each float within `ulps`
/// of its pin (a NaN pin must stay that NaN), each count, flag and text hash exactly. Digests are
/// skipped: one ulp anywhere in a trajectory changes every bit of its digest, so a digest can
/// only be compared bit for bit.
fn assert_within_ulps(group: &str, actual: &Probe, expected: Pinned, ulps: u64) {
    case_shape_or_panic(group, actual, expected);
    let mut moved = Vec::new();
    let mut total = 0;
    for ((name, pins), (_, bits)) in actual.iter().zip(expected) {
        for ((label, pin), &want) in pins.entries.iter().zip(bits.iter()) {
            total += 1;
            let ok = match *pin {
                Pin::D(_) => true,
                Pin::N(value) => value == want,
                Pin::F(value) => {
                    let pinned = f64::from_bits(want);
                    if pinned.is_nan() || value.is_nan() {
                        value.to_bits() == want
                    } else {
                        ulps_apart(value, pinned) <= ulps
                    }
                }
            };
            if !ok {
                let detail = match *pin {
                    Pin::F(value) => format!(
                        "pinned {:?}, now {value:?} ({} ulps apart, allowed {ulps})",
                        f64::from_bits(want),
                        ulps_apart(value, f64::from_bits(want))
                    ),
                    _ => format!("pinned {want}, now {}", describe(*pin)),
                };
                moved.push(format!("  {name} :: {label}: {detail}"));
            }
        }
    }
    assert!(
        moved.is_empty(),
        "{group}: {} of {total} pins moved beyond {ulps} ulps (counts exactly) on a libm the \
         pins were not captured on:\n{}",
        moved.len(),
        moved.join("\n")
    );
}

/// This crate and the released 0.49.0 crate, asked the same questions on this platform, must
/// give the same answers bit for bit: every case, every label, every value, digests included.
#[cfg(not(target_arch = "wasm32"))]
fn assert_same_as_released(group: &str, current: &Probe, released: &Probe) {
    let names =
        |probe: &Probe| -> Vec<String> { probe.iter().map(|(name, _)| name.clone()).collect() };
    assert_eq!(
        names(current),
        names(released),
        "{group}: this engine's and 0.49.0's probes name different cases"
    );
    let mut moved = Vec::new();
    let mut total = 0;
    for ((name, now), (_, then)) in current.iter().zip(released) {
        let labels = |pins: &Pins| -> Vec<String> {
            pins.entries
                .iter()
                .map(|(label, _)| label.clone())
                .collect()
        };
        assert_eq!(
            labels(now),
            labels(then),
            "{group}/{name}: this engine's and 0.49.0's probes record different quantities"
        );
        for ((label, pin), (_, was)) in now.entries.iter().zip(&then.entries) {
            total += 1;
            if pin.bits() != was.bits() {
                moved.push(format!(
                    "  {name} :: {label}: 0.49.0 {:#018x} ({}), this engine {:#018x} ({})",
                    was.bits(),
                    describe(*was),
                    pin.bits(),
                    describe(*pin)
                ));
            }
        }
    }
    assert!(
        moved.is_empty(),
        "{group}: {} of {total} values differ from the released 0.49.0 engine run beside this \
         one on this platform:\n{}",
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

/// How a group is held to its pins (see the module docs).
#[derive(Clone, Copy)]
enum Check {
    /// To the released 0.49.0 crate run beside this one, bit for bit on every platform (`None`
    /// on wasm32, where the old crate is not built and these tests do not run), and to the
    /// literal pins bit for bit on an exact-bit platform.
    Released049(Option<fn() -> Probe>),
    /// To the literal pins: bit for bit on an exact-bit platform, within
    /// [`UNVERIFIED_LIBM_ULPS`] (digests skipped) on any other.
    Pinned,
}

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

// =================================================================================================
// The two engines: this crate, and the released 0.49.0 beside it
// =================================================================================================

/// The probes on this crate, plus the FFI probe, which only this crate can link.
mod current {
    #![allow(dead_code)]
    use super::*;
    use ballistics_engine as engine;

    /// The calls whose signatures 0.50.0 changed, in 0.50.0's form.
    mod shim {
        use super::*;
        use engine::fast_trajectory::FastIntegrationParams;
        use engine::trajectory_integration::TrajectoryParams;
        use engine::{KernelAtmosphere, RelativeHumidity};
        use std::collections::HashMap;

        /// The pinned kernel atmospheres were written as legacy tuples; 0.50.0 moved the kernels
        /// onto `KernelAtmosphere`, and `from_legacy_params` at 0 % RH is the one-line migration
        /// (a dry station keeps its tuple's ratio verbatim; the direct sentinel becomes uniform
        /// air).
        fn kernel_atmosphere(atmos_params: (f64, f64, f64, f64)) -> KernelAtmosphere {
            KernelAtmosphere::from_legacy_params(atmos_params, RelativeHumidity::DRY)
                .expect("pinned atmospheres are valid")
        }

        pub(super) fn trajectory_params(
            atmos_params: (f64, f64, f64, f64),
            atmo_sock: Option<AtmoSock>,
        ) -> TrajectoryParams {
            TrajectoryParams {
                mass_kg: MASS_KG,
                bc: BC_G7,
                drag_model: DragModel::G7,
                wind_segments: crosswind_segments(),
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

        pub(super) fn fast_params(
            atmos_params: (f64, f64, f64, f64),
            atmo_sock: Option<AtmoSock>,
        ) -> FastIntegrationParams {
            FastIntegrationParams {
                horiz: 800.0,
                vert: 0.0,
                initial_state: launch_state(),
                t_span: (0.0, 30.0),
                atmosphere: kernel_atmosphere(atmos_params),
                atmo_sock,
            }
        }

        pub(super) fn compute_derivatives(
            position: Vector3<f64>,
            velocity: Vector3<f64>,
            inputs: &BallisticInputs,
            wind: Vector3<f64>,
            atmos_params: (f64, f64, f64, f64),
            bc: f64,
            atmo_sock: Option<&AtmoSock>,
        ) -> [f64; 6] {
            engine::derivatives::compute_derivatives(
                position,
                velocity,
                inputs,
                wind,
                &kernel_atmosphere(atmos_params),
                bc,
                None,
                0.0,
                atmo_sock,
            )
        }

        /// The legacy binding entry, RK45 to 3 s (it takes no zones).
        pub(super) fn try_solve_trajectory_rust(
            initial_state: [f64; 6],
            wind_segments: Vec<WindSegment>,
            atmos_params: (f64, f64, f64, f64),
        ) -> Result<Vec<HashMap<String, f64>>, engine::wind::WindSegmentError> {
            engine::trajectory_integration::try_solve_trajectory_rust(
                initial_state,
                (0.0, 3.0),
                MASS_KG,
                BC_G7,
                DragModel::G7,
                wind_segments,
                kernel_atmosphere(atmos_params),
                None,
                false,
                false,
                false,
                "RK45".to_string(),
                1e-6,
                0.01,
                800.0,
            )
        }

        pub(super) fn aerodynamic_jump_launch_offset_rad(
            inputs: &BallisticInputs,
            atmos_params: (f64, f64, f64, f64),
        ) -> f64 {
            engine::fast_trajectory::aerodynamic_jump_launch_offset_rad(
                inputs,
                &kernel_atmosphere(atmos_params),
            )
        }

        /// The inputs at 0 % humidity (0.50.0's field is a percent).
        pub(super) fn dry(inputs: BallisticInputs) -> BallisticInputs {
            BallisticInputs {
                humidity_percent: 0.0,
                ..inputs
            }
        }
    }

    include!("humidity_byte_identity/probes.rs");

    #[cfg(feature = "ffi")]
    fn front_end_ffi(humid: bool) -> Probe {
        use engine::ffi::{
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
                    let samples = std::slice::from_raw_parts(
                        r.sampled_points,
                        r.sampled_point_count as usize,
                    );
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

    #[cfg(feature = "ffi")]
    pub(super) fn probe_ffi_trajectory_dry() -> Probe {
        front_end_ffi(false)
    }

    #[cfg(feature = "ffi")]
    pub(super) fn probe_ffi_trajectory_humid() -> Probe {
        front_end_ffi(true)
    }
}

/// The same probes on the released 0.49.0 crate. Its kernels still take the legacy atmosphere
/// tuple and its inputs a humidity fraction, so its shim passes them through unchanged: exactly
/// the calls the byte-identity probes made when these pins were generated on 0.49.0.
#[cfg(not(target_arch = "wasm32"))]
mod released_0_49_0 {
    #![allow(dead_code)]
    use super::*;
    use ballistics_engine_049 as engine;

    /// The calls whose signatures 0.50.0 changed, in 0.49.0's form.
    mod shim {
        use super::*;
        use engine::fast_trajectory::FastIntegrationParams;
        use engine::trajectory_integration::TrajectoryParams;
        use std::collections::HashMap;

        pub(super) fn trajectory_params(
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

        pub(super) fn fast_params(
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

        pub(super) fn compute_derivatives(
            position: Vector3<f64>,
            velocity: Vector3<f64>,
            inputs: &BallisticInputs,
            wind: Vector3<f64>,
            atmos_params: (f64, f64, f64, f64),
            bc: f64,
            atmo_sock: Option<&AtmoSock>,
        ) -> [f64; 6] {
            engine::derivatives::compute_derivatives(
                position,
                velocity,
                inputs,
                wind,
                atmos_params,
                bc,
                None,
                0.0,
                atmo_sock,
            )
        }

        /// The legacy binding entry, RK45 to 3 s (it takes no zones).
        pub(super) fn try_solve_trajectory_rust(
            initial_state: [f64; 6],
            wind_segments: Vec<WindSegment>,
            atmos_params: (f64, f64, f64, f64),
        ) -> Result<Vec<HashMap<String, f64>>, engine::wind::WindSegmentError> {
            engine::trajectory_integration::try_solve_trajectory_rust(
                initial_state,
                (0.0, 3.0),
                MASS_KG,
                BC_G7,
                DragModel::G7,
                wind_segments,
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
        }

        pub(super) fn aerodynamic_jump_launch_offset_rad(
            inputs: &BallisticInputs,
            atmos_params: (f64, f64, f64, f64),
        ) -> f64 {
            engine::fast_trajectory::aerodynamic_jump_launch_offset_rad(inputs, atmos_params)
        }

        /// The inputs at 0 % humidity (0.49.0's field is a fraction).
        pub(super) fn dry(inputs: BallisticInputs) -> BallisticInputs {
            BallisticInputs {
                humidity: 0.0,
                ..inputs
            }
        }
    }

    include!("humidity_byte_identity/probes.rs");
}

/// `Check::Released049` with this group's 0.49.0 probe (none on wasm32, where the 0.49.0 crate is
/// not a dependency and these tests do not run).
#[cfg(not(target_arch = "wasm32"))]
macro_rules! released {
    ($probe:ident) => {
        Check::Released049(Some(released_0_49_0::$probe))
    };
}
#[cfg(target_arch = "wasm32")]
macro_rules! released {
    ($probe:ident) => {
        Check::Released049(None)
    };
}

// =================================================================================================
// The groups, their pins, and the tests
// =================================================================================================

/// One pinned group, as the generator prints it and the tests check it.
struct Group {
    /// The constant the group is pinned in.
    name: &'static str,
    /// The doc comment printed above the constant: what is pinned, and which later step may
    /// change it.
    header: &'static str,
    probe: fn() -> Probe,
    pinned: Pinned,
    check: Check,
    /// The cargo feature the group needs, if any.
    feature: Option<&'static str>,
}

const IMMUTABLE: &str = "IMMUTABLE: nothing in 0.50.0 may change it, the aloft vapour model \
                         included.";
const ALOFT_MODEL_ONLY: &str = "MAY CHANGE ONLY WITH THE ALOFT VAPOUR MODEL: byte-identical \
                                through every other change in 0.50.0.";

/// Every pinned group, in file order.
fn groups() -> Vec<Group> {
    let group = |name, header, probe, pinned, check| Group {
        name,
        header,
        probe,
        pinned,
        check,
        feature: None,
    };
    #[allow(unused_mut)]
    let mut groups = vec![
        group(
            "TS_FLAT_DRY",
            "`TrajectorySolver::solve`, flat fire to 1000 m with sampling: Euler, RK4 and RK45 \
             at ICAO sea level (presence-blind) and at 1500 m, 25 C, 830 hPa (presence-aware), \
             without zones and with two dry zones. RH 0 throughout.",
            current::probe_ts_flat_dry,
            TS_FLAT_DRY,
            released!(probe_ts_flat_dry),
        ),
        group(
            "TS_FLAT_HUMID",
            "`TrajectorySolver::solve`, flat fire, the same methods and stations: station RH 0 \
             under humid zones (RH 20 then 80), and station RH 50 and 100 with and without \
             those zones.",
            current::probe_ts_flat_humid,
            TS_FLAT_HUMID,
            Check::Pinned,
        ),
        group(
            "TS_LOB_DRY",
            "`TrajectorySolver::solve`, the 30-degree lob to ground impact, Euler, RK4 and RK45 \
             at both stations, RH 0, no zones.",
            current::probe_ts_lob_dry,
            TS_LOB_DRY,
            released!(probe_ts_lob_dry),
        ),
        group(
            "TS_LOB_HUMID",
            "`TrajectorySolver::solve`, the 30-degree lob: RH 50 and 100 without zones, and RH 50 \
             under humid zones.",
            current::probe_ts_lob_humid,
            TS_LOB_HUMID,
            Check::Pinned,
        ),
        group(
            "TS_EFFECTS_DRY",
            "`TrajectorySolver::solve` (RK4, flat) with Litz spin drift plus aerodynamic jump, \
             and with Magnus, at both stations, RH 0.",
            current::probe_ts_effects_dry,
            TS_EFFECTS_DRY,
            released!(probe_ts_effects_dry),
        ),
        group(
            "TS_EFFECTS_HUMID",
            "As `TS_EFFECTS_DRY`, at RH 50 and 100.",
            current::probe_ts_effects_humid,
            TS_EFFECTS_HUMID,
            Check::Pinned,
        ),
        group(
            "ZERO_DRY",
            "`calculate_zero_angle_with_station_presence` to 300 m (every zero trial flies \
             `calculate_acceleration`), Euler, RK4 and RK45 at both stations, RH 0. The \
             bisection stops at 1e-7 rad, so a zero pin moves only when a trial height changes \
             sign; the flight groups are the sensitive ones.",
            current::probe_zero_dry,
            ZERO_DRY,
            released!(probe_zero_dry),
        ),
        group(
            "ZERO_HUMID",
            "As `ZERO_DRY`, at RH 50 and 100.",
            current::probe_zero_humid,
            ZERO_HUMID,
            Check::Pinned,
        ),
        group(
            "MONTE_CARLO_DEFAULT_HUMIDITY",
            "Seeded `TrajectorySolver` Monte Carlo at the `BallisticInputs` default humidity \
             (50 %): the CLI's seeded entry at ICAO sea level, and the presence-aware entry the \
             FFI Monte Carlo runs on, at 1500 m.",
            current::probe_monte_carlo_default_humidity,
            MONTE_CARLO_DEFAULT_HUMIDITY,
            Check::Pinned,
        ),
        group(
            "STATION_AIR",
            "Station density and sound speed: `calculate_atmosphere` over altitude, T/P override \
             and RH, plus `calculate_air_density_cimp` and `moist_speed_of_sound`. Station air \
             is not part of the aloft model, so this holds through it at every RH.",
            current::probe_station_air,
            STATION_AIR,
            released!(probe_station_air),
        ),
        group(
            "LOCAL_AIR_DRY",
            "`get_local_atmosphere_humid` from 200 m below each station base to 25 km above it, \
             RH 0.",
            current::probe_local_air_dry,
            LOCAL_AIR_DRY,
            released!(probe_local_air_dry),
        ),
        group(
            "LOCAL_AIR_HUMID",
            "As `LOCAL_AIR_DRY`, at RH 37.5, 50 and 100. The legacy wrapper shares `air_at`'s \
             composition, so the aloft model reaches it.",
            current::probe_local_air_humid,
            LOCAL_AIR_HUMID,
            Check::Pinned,
        ),
        group(
            "DENSITY_ALTITUDE",
            "NWS pressure altitude, the FAA-rule and density-matched density altitudes, the \
             density-altitude entry inverse, and the `atmosphere.density_altitude` service.",
            current::probe_density_altitude,
            DENSITY_ALTITUDE,
            released!(probe_density_altitude),
        ),
        group(
            "STABILITY_SPIN_JUMP",
            "Muzzle Miller Sg (the published dry T/P correction, kept deliberately; a humid \
             correction was rejected), Litz spin drift, the fast path's aerodynamic-jump offset, \
             and `TrajectorySolver`'s aerodynamic-jump components, which depend on station T/P \
             only and so are pinned at every RH.",
            current::probe_stability_spin_jump,
            STABILITY_SPIN_JUMP,
            released!(probe_stability_spin_jump),
        ),
        group(
            "KERNEL_DIRECT",
            "Kernel direct (uniform) mode, density 1.1 kg/m3 and sound speed 335 m/s: \
             `compute_derivatives`, `try_integrate_trajectory` RK4 and RK45, \
             `try_solve_trajectory_rust`, `fast_integrate` and `fast_integrate_with_segments`. \
             0.50.0 migrates the tuple to the uniform `KernelAtmosphere`; the values do not move.",
            current::probe_kernel_direct,
            KERNEL_DIRECT,
            released!(probe_kernel_direct),
        ),
        group(
            "KERNEL_STANDARD_RH0",
            "Kernel standard mode at RH 0, station and zones dry, the same families plus the lean \
             Monte Carlo, at both stations.",
            current::probe_kernel_standard_rh0,
            KERNEL_STANDARD_RH0,
            released!(probe_kernel_standard_rh0),
        ),
        group(
            "FRONT_ENDS_DRY",
            "One in-process call through each JSON front end, RH 0: solve-json decoded and run \
             by `solve_v1` (Euler at sea level, RK45 at 1500 m), and `card_service::come_ups_v1` \
             at both stations.",
            current::probe_front_ends_dry,
            FRONT_ENDS_DRY,
            released!(probe_front_ends_dry),
        ),
        group(
            "FRONT_ENDS_HUMID",
            "As `FRONT_ENDS_DRY` with humidity: solve-json at RH 0.5 (RK4, sea level) and 0.9 \
             (RK45, 1500 m), and the come-ups card at 50 and 100 %.",
            current::probe_front_ends_humid,
            FRONT_ENDS_HUMID,
            Check::Pinned,
        ),
    ];
    #[cfg(feature = "bridge")]
    groups.extend([
        Group {
            feature: Some("bridge"),
            ..group(
                "BRIDGE_SOLVE_DRY",
                "The bridge `solve` command (`bridge_call`) on the `FRONT_ENDS_DRY` requests.",
                current::probe_bridge_solve_dry,
                BRIDGE_SOLVE_DRY,
                released!(probe_bridge_solve_dry),
            )
        },
        Group {
            feature: Some("bridge"),
            ..group(
                "BRIDGE_SOLVE_HUMID",
                "The bridge `solve` command on the `FRONT_ENDS_HUMID` requests.",
                current::probe_bridge_solve_humid,
                BRIDGE_SOLVE_HUMID,
                Check::Pinned,
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
                current::probe_ffi_trajectory_dry,
                FFI_TRAJECTORY_DRY,
                // 0.49.0's C ABI exports the same #[no_mangle] symbols as this crate's, so the
                // two cannot be linked into one test binary: no runtime oracle.
                Check::Pinned,
            )
        },
        Group {
            feature: Some("ffi"),
            ..group(
                "FFI_TRAJECTORY_HUMID",
                "As `FFI_TRAJECTORY_DRY`, at RH 50 and 100.",
                current::probe_ffi_trajectory_humid,
                FFI_TRAJECTORY_HUMID,
                Check::Pinned,
            )
        },
    ]);
    groups
}

/// Whether a group's constant may change with the aloft vapour model: every group whose cases
/// involve an RH above
/// zero, at the station or in a zone, except the station-air and stability groups.
fn aloft_model_may_change(name: &str) -> bool {
    name.ends_with("_HUMID") || name == "MONTE_CARLO_DEFAULT_HUMIDITY"
}

fn group_named(name: &str) -> Group {
    groups()
        .into_iter()
        .find(|group| group.name == name)
        .unwrap_or_else(|| panic!("no group {name}"))
}

/// Hold one group to its pins the way its [`Check`] says, on this platform.
fn check_group(name: &str) {
    let group = group_named(name);
    let probe = (group.probe)();
    let exact_bit_platform = Libm::of_this_build().is_some();
    match group.check {
        Check::Released049(released) => {
            #[cfg(not(target_arch = "wasm32"))]
            assert_same_as_released(
                name,
                &probe,
                &released.expect("the 0.49.0 probe exists off wasm32")(),
            );
            #[cfg(target_arch = "wasm32")]
            let _ = released;
            if exact_bit_platform {
                assert_bits(name, &probe, group.pinned);
            }
        }
        Check::Pinned => {
            if exact_bit_platform {
                assert_bits(name, &probe, group.pinned);
            } else {
                assert_within_ulps(name, &probe, group.pinned, UNVERIFIED_LIBM_ULPS);
            }
        }
    }
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
        print_doc(if aloft_model_may_change(group.name) {
            ALOFT_MODEL_ONLY
        } else {
            IMMUTABLE
        });
        if let Some(feature) = group.feature {
            println!("#[cfg(feature = {feature:?})]");
        }
        // rustfmt would put every 18-character literal on its own line.
        println!("#[rustfmt::skip]");
        print_probe(group.name, &probe);
    }
}

macro_rules! pin_test {
    ($test:ident, $group:ident) => {
        #[test]
        fn $test() {
            check_group(stringify!($group));
        }
    };
}

pin_test!(trajectory_solver_flat_dry_is_byte_identical, TS_FLAT_DRY);
pin_test!(
    trajectory_solver_flat_humid_is_byte_identical,
    TS_FLAT_HUMID
);
pin_test!(trajectory_solver_lob_dry_is_byte_identical, TS_LOB_DRY);
pin_test!(trajectory_solver_lob_humid_is_byte_identical, TS_LOB_HUMID);
pin_test!(
    trajectory_solver_effects_dry_are_byte_identical,
    TS_EFFECTS_DRY
);
pin_test!(
    trajectory_solver_effects_humid_are_byte_identical,
    TS_EFFECTS_HUMID
);
pin_test!(zero_angles_dry_are_byte_identical, ZERO_DRY);
pin_test!(zero_angles_humid_are_byte_identical, ZERO_HUMID);
pin_test!(
    seeded_monte_carlo_at_default_humidity_is_byte_identical,
    MONTE_CARLO_DEFAULT_HUMIDITY
);
pin_test!(station_air_is_byte_identical, STATION_AIR);
pin_test!(local_air_dry_is_byte_identical, LOCAL_AIR_DRY);
pin_test!(local_air_humid_is_byte_identical, LOCAL_AIR_HUMID);
pin_test!(density_altitude_is_byte_identical, DENSITY_ALTITUDE);
pin_test!(
    stability_spin_drift_and_jump_are_byte_identical,
    STABILITY_SPIN_JUMP
);
pin_test!(kernel_direct_mode_is_byte_identical, KERNEL_DIRECT);
pin_test!(
    kernel_standard_mode_at_rh0_is_byte_identical,
    KERNEL_STANDARD_RH0
);
pin_test!(front_ends_dry_are_byte_identical, FRONT_ENDS_DRY);
pin_test!(front_ends_humid_are_byte_identical, FRONT_ENDS_HUMID);
#[cfg(feature = "bridge")]
pin_test!(bridge_solve_dry_is_byte_identical, BRIDGE_SOLVE_DRY);
#[cfg(feature = "bridge")]
pin_test!(bridge_solve_humid_is_byte_identical, BRIDGE_SOLVE_HUMID);
#[cfg(feature = "ffi")]
pin_test!(ffi_trajectory_dry_is_byte_identical, FFI_TRAJECTORY_DRY);
#[cfg(feature = "ffi")]
pin_test!(ffi_trajectory_humid_is_byte_identical, FFI_TRAJECTORY_HUMID);

/// Every group is checked by exactly one test: a group added to `groups()` without a
/// `pin_test!` would be generated and never compared.
#[test]
fn every_group_has_a_test() {
    let source = include_str!("humidity_byte_identity.rs");
    for group in groups() {
        let needle = format!(", {});", group.name);
        let tests = source.matches(&needle).count()
            + source.matches(&format!("    {}\n);", group.name)).count();
        assert_eq!(tests, 1, "{}: {tests} pin_test! entries", group.name);
    }
}

/// The groups held to 0.49.0 are exactly the ones whose promise is "0.50.0 equals 0.49.0": never
/// a group the aloft vapour model may change (0.49.0 is not its oracle), and every IMMUTABLE
/// group except the one
/// 0.49.0 cannot be linked for (its C ABI).
#[test]
fn every_immutable_group_with_a_0_49_0_entry_is_held_to_0_49_0() {
    for group in groups() {
        let held_to_released = matches!(group.check, Check::Released049(_));
        if aloft_model_may_change(group.name) {
            assert!(
                !held_to_released,
                "{}: the aloft model may change it",
                group.name
            );
        } else if group.name != "FFI_TRAJECTORY_DRY" {
            assert!(
                held_to_released,
                "{}: IMMUTABLE, and 0.49.0 has it",
                group.name
            );
        }
    }
}

/// The released crate really is another engine: where the aloft vapour model moved the physics,
/// it disagrees with
/// this one. (Were the dev-dependency ever resolved to this crate itself, every
/// `Check::Released049` comparison would pass vacuously.)
#[cfg(not(target_arch = "wasm32"))]
#[test]
fn the_released_engine_is_0_49_0_and_not_this_one() {
    let now = current::probe_local_air_humid();
    let then = released_0_49_0::probe_local_air_humid();
    let differing = now
        .iter()
        .zip(&then)
        .flat_map(|((_, a), (_, b))| a.entries.iter().zip(&b.entries))
        .filter(|((_, a), (_, b))| a.bits() != b.bits())
        .count();
    assert!(
        differing > 0,
        "0.49.0's aloft air equals this engine's: the oracle is not 0.49.0"
    );
    // And it agrees exactly with 0.49.0's own pins where those were captured.
    if Libm::of_this_build().is_some() {
        assert_bits(
            "LOCAL_AIR_DRY",
            &released_0_49_0::probe_local_air_dry(),
            LOCAL_AIR_DRY,
        );
    }
}

/// The released 0.49.0 engine is in Cargo.lock as a second `ballistics-engine` package, and the
/// release bump moves only this crate. At the 0.49.0 -> 0.50.0 bump both copies are 0.49.0, so a
/// root `cargo update -p ballistics-engine@0.49.0` (bump.sh's old line, and its fallback without
/// the version) is ambiguous and would stop `scripts/release/bump.sh`; `--workspace` is not.
#[test]
fn the_release_bump_moves_only_this_crate_and_never_the_0_49_0_oracle() {
    // Compiled in, so the test needs no checkout where it runs (a cross-built test binary).
    let lock = include_str!("../Cargo.lock");
    let engines: Vec<&str> = lock
        .split("[[package]]")
        .filter(|package| package.contains("\nname = \"ballistics-engine\"\n"))
        .collect();
    assert_eq!(
        engines.len(),
        2,
        "Cargo.lock: this crate and the 0.49.0 oracle"
    );
    assert_eq!(
        engines
            .iter()
            .filter(|package| package.contains("\nversion = \"0.49.0\"\n")
                && package.contains("\nsource = \"registry+"))
            .count(),
        1,
        "Cargo.lock: the oracle is the released 0.49.0 from crates.io"
    );
    let bump = include_str!("../scripts/release/bump.sh");
    let root_updates: Vec<&str> = bump
        .lines()
        .map(str::trim)
        .filter(|line| line.starts_with("cargo update") && !line.contains("--manifest-path"))
        .collect();
    assert_eq!(
        root_updates,
        ["cargo update --workspace"],
        "scripts/release/bump.sh must move Cargo.lock's own crate with --workspace alone"
    );
}

/// Each `LIBM_ALTERNATES` entry names a value the probe records, in a group the aloft vapour
/// model may re-pin, and
/// differs from that value's primary pin, so no alternate can widen a pin that does not need it
/// or outlive the pin it was recorded for.
#[test]
fn every_libm_alternate_names_a_pinned_value_it_differs_from() {
    for &(group_name, case, label, _, bits) in LIBM_ALTERNATES {
        assert!(
            aloft_model_may_change(group_name),
            "{group_name} is IMMUTABLE"
        );
        let group = group_named(group_name);
        let probe = (group.probe)();
        let (case_index, (_, pins)) = probe
            .iter()
            .enumerate()
            .find(|(_, (name, _))| name == case)
            .unwrap_or_else(|| panic!("{group_name}: no case {case}"));
        let label_index = pins
            .entries
            .iter()
            .position(|(entry, _)| entry == label)
            .unwrap_or_else(|| panic!("{group_name}/{case}: no {label}"));
        let (pinned_case, pinned_bits) = group.pinned[case_index];
        assert_eq!(pinned_case, case);
        assert_ne!(
            pinned_bits[label_index], bits,
            "{group_name}/{case} :: {label}: the alternate is the primary pin"
        );
    }
}

/// The ulp measure the unverified-libm comparison uses.
#[test]
fn ulps_apart_counts_representable_values() {
    assert_eq!(ulps_apart(1.0, 1.0), 0);
    assert_eq!(ulps_apart(0.0, -0.0), 0);
    assert_eq!(ulps_apart(1.0, f64::from_bits(1.0f64.to_bits() + 2)), 2);
    assert_eq!(ulps_apart(-1.0, f64::from_bits((-1.0f64).to_bits() + 3)), 3);
    assert_eq!(
        ulps_apart(f64::MIN_POSITIVE, -f64::MIN_POSITIVE),
        2 * f64::MIN_POSITIVE.to_bits()
    );
    assert_eq!(ulps_apart(5e-324, -5e-324), 2);
}

// ----- PINNED CONSTANTS (generated; see the module docs) -----------------------------------------

/// `TrajectorySolver::solve`, flat fire to 1000 m with sampling: Euler, RK4 and RK45 at ICAO sea
/// level (presence-blind) and at 1500 m, 25 C, 830 hPa (presence-aware), without zones and with two
/// dry zones. RH 0 throughout.
///
/// IMMUTABLE: nothing in 0.50.0 may change it, the aloft vapour model included.
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
/// MAY CHANGE ONLY WITH THE ALOFT VAPOUR MODEL: byte-identical through every other change in
/// 0.50.0.
///
/// Re-pinned by the aloft vapour model (0.50.0): the unified aloft vapour model holds the
/// governing air's (the station's or the zone's) water-vapour mole fraction for the lapse, the
/// density and the sound speed, with a virtual-temperature lapse. 925 of 1620 values moved, all in
/// cases with a humidity above 0 and all off the station's own height; every RH-0 and
/// station-height value is unchanged. The largest float shift of each moved case,
/// pinned before the aloft vapour model -> pinned now:
/// - `flat/euler/icao-sea-level/rh0/zones-rh20-rh80` points[1964].drag_coefficient: 0.21625498159624057 -> 0.21625065852669897 (-2.00e-05)
/// - `flat/euler/icao-sea-level/rh50/no-zones` points[1964].drag_coefficient: 0.22335507451629613 -> 0.22335131439114955 (-1.68e-05)
/// - `flat/euler/icao-sea-level/rh50/zones-rh20-rh80` points[1964].drag_coefficient: 0.21316524983232532 -> 0.21316121741649008 (-1.89e-05)
/// - `flat/euler/icao-sea-level/rh100/no-zones` points[1961].drag_coefficient: 0.2257665106154504 -> 0.2257587849895006 (-3.42e-05)
/// - `flat/euler/icao-sea-level/rh100/zones-rh20-rh80` points[1964].drag_coefficient: 0.2102937767718652 -> 0.2102900132956857 (-1.79e-05)
/// - `flat/euler/1500m-25C-830hPa/rh0/zones-rh20-rh80` impact_energy: 898.7868653722531 -> 898.793490741712 (+7.37e-06)
/// - `flat/euler/1500m-25C-830hPa/rh50/no-zones` impact_energy: 902.8418441050804 -> 902.8467451689961 (+5.43e-06)
/// - `flat/euler/1500m-25C-830hPa/rh50/zones-rh20-rh80` impact_energy: 898.7868653722531 -> 898.793490741712 (+7.37e-06)
/// - `flat/euler/1500m-25C-830hPa/rh100/no-zones` impact_energy: 912.0087179015902 -> 912.0184825693569 (+1.07e-05)
/// - `flat/euler/1500m-25C-830hPa/rh100/zones-rh20-rh80` impact_energy: 898.7868653722531 -> 898.793490741712 (+7.37e-06)
/// - `flat/rk4/icao-sea-level/rh0/zones-rh20-rh80` points[1963].drag_coefficient: 0.217314962807575 -> 0.21731058558099348 (-2.01e-05)
/// - `flat/rk4/icao-sea-level/rh50/no-zones` points[1963].drag_coefficient: 0.22457758794752264 -> 0.22457375325337803 (-1.71e-05)
/// - `flat/rk4/icao-sea-level/rh50/zones-rh20-rh80` points[1963].drag_coefficient: 0.21415394908923926 -> 0.21414986626609883 (-1.91e-05)
/// - `flat/rk4/icao-sea-level/rh100/no-zones` points[1960].drag_coefficient: 0.22706225120487242 -> 0.22705440142370198 (-3.46e-05)
/// - `flat/rk4/icao-sea-level/rh100/zones-rh20-rh80` points[1963].drag_coefficient: 0.211216467020616 -> 0.21121265703273817 (-1.80e-05)
/// - `flat/rk4/1500m-25C-830hPa/rh0/zones-rh20-rh80` impact_energy: 899.6211801264712 -> 899.6277959956163 (+7.35e-06)
/// - `flat/rk4/1500m-25C-830hPa/rh50/no-zones` impact_energy: 903.684569602166 -> 903.6894628571841 (+5.41e-06)
/// - `flat/rk4/1500m-25C-830hPa/rh50/zones-rh20-rh80` impact_energy: 899.6211801264712 -> 899.6277959956163 (+7.35e-06)
/// - `flat/rk4/1500m-25C-830hPa/rh100/no-zones` impact_energy: 912.8490094509566 -> 912.8587585221065 (+1.07e-05)
/// - `flat/rk4/1500m-25C-830hPa/rh100/zones-rh20-rh80` impact_energy: 899.6211801264712 -> 899.6277959956163 (+7.35e-06)
/// - `flat/rk45/icao-sea-level/rh0/zones-rh20-rh80` points[200].drag_coefficient: 0.21731082853108347 -> 0.2173064482450904 (-2.02e-05)
/// - `flat/rk45/icao-sea-level/rh50/no-zones` points[200].drag_coefficient: 0.22457914193212117 -> 0.22457530522455865 (-1.71e-05)
/// - `flat/rk45/icao-sea-level/rh50/zones-rh20-rh80` points[200].drag_coefficient: 0.21415009287585954 -> 0.21414600719892168 (-1.91e-05)
/// - `flat/rk45/icao-sea-level/rh100/no-zones` points[199].drag_coefficient: 0.22707757431859077 -> 0.2270697331623272 (-3.45e-05)
/// - `flat/rk45/icao-sea-level/rh100/zones-rh20-rh80` points[200].drag_coefficient: 0.2112128684993169 -> 0.21120905584658373 (-1.81e-05)
/// - `flat/rk45/1500m-25C-830hPa/rh0/zones-rh20-rh80` impact_energy: 899.6170526085074 -> 899.6236682738398 (+7.35e-06)
/// - `flat/rk45/1500m-25C-830hPa/rh50/no-zones` impact_energy: 903.6857804389596 -> 903.6906733452946 (+5.41e-06)
/// - `flat/rk45/1500m-25C-830hPa/rh50/zones-rh20-rh80` impact_energy: 899.6170526085074 -> 899.6236682738398 (+7.35e-06)
/// - `flat/rk45/1500m-25C-830hPa/rh100/no-zones` impact_energy: 912.8500651048803 -> 912.8598149896424 (+1.07e-05)
/// - `flat/rk45/1500m-25C-830hPa/rh100/zones-rh20-rh80` impact_energy: 899.6170526085074 -> 899.6236682738398 (+7.35e-06)
#[rustfmt::skip]
const TS_FLAT_HUMID: Pinned = &[
    (
        "flat/euler/icao-sea-level/rh0/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc06e2668e4c2c1, 0x3fff69982140ab5b, 0x40742705aa527cf0,
            0x40826bef80ee7ffd, 0x407544b4bbbd89be, 0x0000000000000001, 0x408830faa1bc8c87,
            0x0000000000000001, 0x408d39ec719e29a9, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ad,
            0x3f50624dd2f1a9fc, 0x3fea51c1c06bf31a, 0x3f5aca57130f1490, 0xbec5216177e29fc7,
            0x4089b3da8a188477, 0x0000000000000001, 0x3fd30dc2e8c659d5, 0x3fdf6c8b4395810d,
            0x4075e18ed24d1683, 0xbfd84ad18a6e2c07, 0xbfd0cba65bbd9fcc, 0x4083560a1f175047,
            0x0000000000000001, 0x3fd58cf605f7dc58, 0x3fef6c8b4395810d, 0x40835108c1df67d3,
            0xc0060ec47226dea4, 0xbfed8f1e3dc8317c, 0x407e15d478cd5914, 0x0000000000000001,
            0x3fd8b936baacb482, 0x3ff7916872b01fdd, 0x4089e8fe8fdf9734, 0xc01af68c8ed4dbfd,
            0xbffdcbadb767af74, 0x4077ff00c5a2823b, 0x0000000000000001, 0x3fdbab56b02436e7,
            0x3fff69982140ab5b, 0x408f400000000000, 0xc028394aa818fb7d, 0xc007f15755620e81,
            0x40742705aa527cf0, 0x0000000000000001, 0x3fcbae1a010e294f, 0xfe49040c2f14f96e,
            0x0000000000000015, 0x6b27d4b4a197a35a,
        ],
    ),
    (
        "flat/euler/icao-sea-level/rh50/no-zones",
        &[
            0x408f400000000000, 0x3fc06c4a2b3ff6ad, 0x3fff69ff2f0497e2, 0x40743c219f9dc3e4,
            0x4082929b62dcdb2c, 0x40754bf73b88949e, 0x0000000000000001, 0x40882c3eb1ca7a25,
            0x0000000000000001, 0x408d42b5c5188542, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ad,
            0x3f50624dd2f1a9fc, 0x3fea51ba3e4e3ed6, 0x3f5aca4f62cf3933, 0xbec546bfb16f49c3,
            0x4089b3d335068e3e, 0x0000000000000001, 0x3fd310d7d819611a, 0x3fdf6c8b4395810d,
            0x4075dc5bc4af9a7e, 0xbfd84b1ea6090229, 0xbfd0e586e892219d, 0x40834ce90c6079fd,
            0x0000000000000001, 0x3fd59516b20f8447, 0x3fef6c8b4395810d, 0x40834ba7dfab929e,
            0xc00610a3dd7bcc69, 0xbfeda9e2e2cdbfc6, 0x407e14639199e769, 0x0000000000000001,
            0x3fd8beafa587f911, 0x3ff7916872b01fdd, 0x4089e4eb8a69989c, 0xc01afe4a19811d74,
            0xbffdd5d13c36bd13, 0x407809d542e75e4c, 0x0000000000000001, 0x3fdbaa4c2766eaa0,
            0x3fff69ff2f0497e2, 0x408f400000000000, 0xc028450c97b06431, 0xc007f22570e9e789,
            0x40743c219f9dc3e4, 0x0000000000000001, 0x3fcc96c69f6a0f9d, 0x06609c64211bb41a,
            0x0000000000000015, 0x36d8182c7d344e2f,
        ],
    ),
    (
        "flat/euler/icao-sea-level/rh50/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc06e2668e4c2c1, 0x3fff69982140ab5b, 0x40742705aa527cf0,
            0x40826bef80ee7ffd, 0x40754bf73b88949e, 0x0000000000000001, 0x4088272dd26cd351,
            0x0000000000000001, 0x408d2f0967953393, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ad,
            0x3f50624dd2f1a9fc, 0x3fea51c1c06bf31a, 0x3f5aca57130f1490, 0xbec5216177e29fc7,
            0x4089b3da8a188477, 0x0000000000000001, 0x3fd310d546ce290d, 0x3fdf6c8b4395810d,
            0x4075e18ed24d1683, 0xbfd84ad18a6e2c07, 0xbfd0cba65bbd9fcc, 0x4083560a1f175047,
            0x0000000000000001, 0x3fd5905c171f68c8, 0x3fef6c8b4395810d, 0x40835108c1df67d3,
            0xc0060ec47226dea4, 0xbfed8f1e3dc8317c, 0x407e15d478cd5914, 0x0000000000000001,
            0x3fd8be0323fd235d, 0x3ff7916872b01fdd, 0x4089e8fe8fdf9734, 0xc01af68c8ed4dbfd,
            0xbffdcbadb767af74, 0x4077ff00c5a2823b, 0x0000000000000001, 0x3fdbae8c37d14fc0,
            0x3fff69982140ab5b, 0x408f400000000000, 0xc028394aa818fb7d, 0xc007f15755620e81,
            0x40742705aa527cf0, 0x0000000000000001, 0x3fcb48dde4ca28d0, 0xfe49040c2f14f96e,
            0x0000000000000015, 0x6b27d4b4a197a35a,
        ],
    ),
    (
        "flat/euler/icao-sea-level/rh100/no-zones",
        &[
            0x408f400000000000, 0x3fc06cfb6bbd0f41, 0x3fff5ddd350c5d0f, 0x407447624a01d315,
            0x4082a749318eeb75, 0x40755334ef1fc968, 0x0000000000000001, 0x408830eb0b2f7656,
            0x0000000000000001, 0x408d4c41ee4bf5d0, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007aa,
            0x3f50624dd2f1a9fc, 0x3fea51bd09223bce, 0x3f5aca523ecb5d7b, 0xbec538da0c5e0036,
            0x4089b3d5ef19ee17, 0x0000000000000001, 0x3fd313e577ef400d, 0x3fdf5c28f5c28f63,
            0x4075d468fa986e4a, 0xbfd81a610f6d750c, 0xbfd0cb8caf224857, 0x408352faab8389cb,
            0x0000000000000001, 0x3fd5955848e08846, 0x3fef645a1cac0838, 0x40834acac88e518c,
            0xc00604908a038129, 0xbfed8d6aa0a0c74b, 0x407e214eb56d5f3f, 0x0000000000000001,
            0x3fd8bd6dedd16f93, 0x3ff789374bc6a709, 0x4089e41d261446b3, 0xc01aed3e6664cb9e,
            0xbffdb70e39c3db4a, 0x407818fd06909e1a, 0x0000000000000001, 0x3fdba787b5d4d37a,
            0x3fff5ddd350c5d0f, 0x408f400000000000, 0xc02834aee34646dc, 0xc007d9e17cf971c8,
            0x407447624a01d315, 0x0000000000000001, 0x3fcce5a9f32844d7, 0xfc2c4a1a0e858782,
            0x0000000000000015, 0x318ae756e1616dc7,
        ],
    ),
    (
        "flat/euler/icao-sea-level/rh100/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc06e2668e4c2c1, 0x3fff69982140ab5b, 0x40742705aa527cf0,
            0x40826bef80ee7ffd, 0x40755334ef1fc968, 0x0000000000000001, 0x40881d5d9eaeeb84,
            0x0000000000000001, 0x408d26dc5fabd125, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ad,
            0x3f50624dd2f1a9fc, 0x3fea51c1c06bf31a, 0x3f5aca57130f1490, 0xbec5216177e29fc7,
            0x4089b3da8a188477, 0x0000000000000001, 0x3fd313e3db5dfc03, 0x3fdf6c8b4395810d,
            0x4075e18ed24d1683, 0xbfd84ad18a6e2c07, 0xbfd0cba65bbd9fcc, 0x4083560a1f175047,
            0x0000000000000001, 0x3fd593c1e49c26e3, 0x3fef6c8b4395810d, 0x40835108c1df67d3,
            0xc0060ec47226dea4, 0xbfed8f1e3dc8317c, 0x407e15d478cd5914, 0x0000000000000001,
            0x3fd8c2cac128e4e2, 0x3ff7916872b01fdd, 0x4089e8fe8fdf9734, 0xc01af68c8ed4dbfd,
            0xbffdcbadb767af74, 0x4077ff00c5a2823b, 0x0000000000000001, 0x3fdbb1b4dd842df7,
            0x3fff69982140ab5b, 0x408f400000000000, 0xc028394aa818fb7d, 0xc007f15755620e81,
            0x40742705aa527cf0, 0x0000000000000001, 0x3fcaeac87ce3e356, 0xfe49040c2f14f96e,
            0x0000000000000015, 0x6b27d4b4a197a35a,
        ],
    ),
    (
        "flat/euler/1500m-25C-830hPa/rh0/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc0a5dc4acc5ccf, 0x3ffbebc5b9c50b37, 0x4078e24744fd68a8,
            0x408c165911ac8aa8, 0x4075a2608de6e31f, 0x0000000000000001, 0x408dc0324e6c6ff9,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d3,
            0x3f50624dd2f1a9fc, 0x3fea529f2159466f, 0x3f5acb39c436ada5, 0xbec0d398a2a9eea0,
            0x4089b4b2badc6d3b, 0x0000000000000001, 0x3fd3348d527fcc08, 0x3fdbe76c8b439587,
            0x407437d3dc6db77b, 0xbfccf8af7783b24d, 0xbfc5fbfe86ac25e8, 0x4084ec66fc813ab5,
            0x0000000000000001, 0x3fd4fa2f98be2246, 0x3febef9db22d0e5c, 0x40826f9b23fdf37f,
            0xc000efee95ea4c55, 0xbfe3fd5853dd6622, 0x40815a5f4acf57ee, 0x0000000000000001,
            0x3fd6fe8de3a6ade4, 0x3ff4f5c28f5c285f, 0x408963d758ab2599, 0xc015a853b5a7b54b,
            0xbff4a86ee65d1197, 0x407d315ac9067ec2, 0x0000000000000001, 0x3fd969137b1f2f42,
            0x3ffbebc5b9c50b37, 0x408f400000000000, 0xc023d2198f04e217, 0xc000f5b2866acd63,
            0x4078e24744fd68a8, 0x0000000000000001, 0x3fdb78d4e52ad476, 0xc6b3687204f95b8a,
            0x0000000000000015, 0xef8105e6b8e9a2df,
        ],
    ),
    (
        "flat/euler/1500m-25C-830hPa/rh50/no-zones",
        &[
            0x408f400000000000, 0x3fc0a49e7aafae4e, 0x3ffbebdc49dc8e2e, 0x4078f0a04a9ca907,
            0x408c36c62254c713, 0x4075b318f42cf8e0, 0x0000000000000001, 0x408db55540556944,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d3,
            0x3f50624dd2f1a9fc, 0x3fea529a4a656351, 0x3f5acb34cf863b98, 0xbec0ebaecef0e832,
            0x4089b4ae00f1a3cd, 0x0000000000000001, 0x3fd33b66fe995b33, 0x3fdbe76c8b439587,
            0x407434f96d5480bd, 0xbfccfa49051362ad, 0xbfc618652786ab50, 0x4084e685f419f85b,
            0x0000000000000001, 0x3fd5031267a231a7, 0x3febef9db22d0e5c, 0x40826c3e568dd676,
            0xc000f0acec967ee3, 0xbfe40e1479d5041b, 0x4081596904ea1b47, 0x0000000000000001,
            0x3fd7092ecf45cab5, 0x3ff4f5c28f5c285f, 0x4089614288ed78e7, 0xc015ac591633a4de,
            0xbff4aedb4dd4ebdb, 0x407d398cf9dfb329, 0x0000000000000001, 0x3fd9700aff6c0c22,
            0x3ffbebdc49dc8e2e, 0x408f400000000000, 0xc023d80a7ca670d0, 0xc000f5dfa699d32d,
            0x4078f0a04a9ca907, 0x0000000000000001, 0x3fdb7af9e60c95d4, 0x8cc5c6dafe2c1675,
            0x0000000000000015, 0xddc6a919748a9d07,
        ],
    ),
    (
        "flat/euler/1500m-25C-830hPa/rh50/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc0a5dc4acc5ccf, 0x3ffbebc5b9c50b37, 0x4078e24744fd68a8,
            0x408c165911ac8aa8, 0x4075b318f42cf8e0, 0x0000000000000001, 0x408da244674db21c,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d3,
            0x3f50624dd2f1a9fc, 0x3fea529f2159466f, 0x3f5acb39c436ada5, 0xbec0d398a2a9eea0,
            0x4089b4b2badc6d3b, 0x0000000000000001, 0x3fd33b655d058486, 0x3fdbe76c8b439587,
            0x407437d3dc6db77b, 0xbfccf8af7783b24d, 0xbfc5fbfe86ac25e8, 0x4084ec66fc813ab5,
            0x0000000000000001, 0x3fd500b127908a44, 0x3febef9db22d0e5c, 0x40826f9b23fdf37f,
            0xc000efee95ea4c55, 0xbfe3fd5853dd6622, 0x40815a5f4acf57ee, 0x0000000000000001,
            0x3fd7087775d585a8, 0x3ff4f5c28f5c285f, 0x408963d758ab2599, 0xc015a853b5a7b54b,
            0xbff4a86ee65d1197, 0x407d315ac9067ec2, 0x0000000000000001, 0x3fd974012d108fad,
            0x3ffbebc5b9c50b37, 0x408f400000000000, 0xc023d2198f04e217, 0xc000f5b2866acd63,
            0x4078e24744fd68a8, 0x0000000000000001, 0x3fdb81388955b4b5, 0xc6b3687204f95b8a,
            0x0000000000000015, 0xef8105e6b8e9a2df,
        ],
    ),
    (
        "flat/euler/1500m-25C-830hPa/rh100/no-zones",
        &[
            0x408f400000000000, 0x3fc0a5ed39d0f6e5, 0x3ffbd9b96d8b7ac8, 0x407910f951a9a35c,
            0x408c8025da307779, 0x4075c3b860e9c586, 0x0000000000000001, 0x408dc5812cf98898,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006ce,
            0x3f50624dd2f1a9fc, 0x3fea529f6372da24, 0x3f5acb3a07e66570, 0xbec0d24fa9443d64,
            0x4089b4b2fb696fe3, 0x0000000000000001, 0x3fd34230642e05f3, 0x3fdbd70a3d70a3dd,
            0x40742d475c528ef6, 0xbfcca59fb49013fb, 0xbfc5e1ec7687a00e, 0x4084ef25cc58363c,
            0x0000000000000001, 0x3fd50610da3ba399, 0x3febdf3b645a1cb2, 0x408268532f9efe50,
            0xc000d9da60c7b42a, 0xbfe3e00bdf0b0966, 0x4081664d395f00da, 0x0000000000000001,
            0x3fd7097aa5590e35, 0x3ff4e560418936b7, 0x40895ba36dd38d1c, 0xc0158b33aa9973d5,
            0xbff47b4f3f80ddc1, 0x407d5a4770dbb57e, 0x0000000000000001, 0x3fd96b0bccd1d91d,
            0x3ffbd9b96d8b7ac8, 0x408f400000000000, 0xc023c20a6126d0f1, 0xc000d199edf7ac75,
            0x407910f951a9a35c, 0x0000000000000001, 0x3fdb752af2aa5e9a, 0x41dce89e6a3f6b2f,
            0x0000000000000015, 0xc31fe96a37d15305,
        ],
    ),
    (
        "flat/euler/1500m-25C-830hPa/rh100/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc0a5dc4acc5ccf, 0x3ffbebc5b9c50b37, 0x4078e24744fd68a8,
            0x408c165911ac8aa8, 0x4075c3b860e9c586, 0x0000000000000001, 0x408d8794dcee1aec,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d3,
            0x3f50624dd2f1a9fc, 0x3fea529f2159466f, 0x3f5acb39c436ada5, 0xbec0d398a2a9eea0,
            0x4089b4b2badc6d3b, 0x0000000000000001, 0x3fd342307a83c9c7, 0x3fdbe76c8b439587,
            0x407437d3dc6db77b, 0xbfccf8af7783b24d, 0xbfc5fbfe86ac25e8, 0x4084ec66fc813ab5,
            0x0000000000000001, 0x3fd5072ebd5cd110, 0x3febef9db22d0e5c, 0x40826f9b23fdf37f,
            0xc000efee95ea4c55, 0xbfe3fd5853dd6622, 0x40815a5f4acf57ee, 0x0000000000000001,
            0x3fd71265f2b15f6e, 0x3ff4f5c28f5c285f, 0x408963d758ab2599, 0xc015a853b5a7b54b,
            0xbff4a86ee65d1197, 0x407d315ac9067ec2, 0x0000000000000001, 0x3fd97ebd612242fe,
            0x3ffbebc5b9c50b37, 0x408f400000000000, 0xc023d2198f04e217, 0xc000f5b2866acd63,
            0x4078e24744fd68a8, 0x0000000000000001, 0x3fdb895ca2686553, 0xc6b3687204f95b8a,
            0x0000000000000015, 0xef8105e6b8e9a2df,
        ],
    ),
    (
        "flat/rk4/icao-sea-level/rh0/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc08a6b00758fe8, 0x3fff64805716bcfa, 0x4074294661cf1bb2,
            0x4082700e1a93c7aa, 0x407544b4bbbd89be, 0x0000000000000001, 0x40883366a41454d4,
            0x0000000000000001, 0x408d3fa8404a6374, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ac,
            0x3f50624dd2f1a9fc, 0x3fea53e15f4bb1bb, 0x3f5ae115c0abe4b4, 0xbeb51fa51d33e34e,
            0x4089b3db0cdf2bba, 0x0000000000000001, 0x3fd30dc2baedb49a, 0x3fdf6c8b4395810d,
            0x4075e389b2ce609b, 0xbfd81fb74c6379a9, 0xbfd0c1cbab9c2531, 0x408356a1eb40d93f,
            0x0000000000000001, 0x3fd58ca7f588d8fc, 0x3fef6c8b4395810d, 0x408352e6e46d24bb,
            0xc0060439faf6497a, 0xbfed85d2999d1a0d, 0x407e17490872d6d2, 0x0000000000000001,
            0x3fd8b8888b8a5d0d, 0x3ff7916872b01fdd, 0x4089eb9e17f94a79, 0xc01aeed5186acb14,
            0xbffdc526a2b3f863, 0x4078007dc7106e9d, 0x0000000000000001, 0x3fdbaac04aa82659,
            0x3fff64805716bcfa, 0x408f400000000000, 0xc0282c2f9aa338f2, 0xc007e727c10e31b2,
            0x4074294661cf1bb2, 0x0000000000000001, 0x3fcbd0d551128e94, 0x392f07b4fa65a45e,
            0x0000000000000015, 0x2a32fcc0f491b817,
        ],
    ),
    (
        "flat/rk4/icao-sea-level/rh50/no-zones",
        &[
            0x408f400000000000, 0x3fc0888f190e8a73, 0x3fff64e2a47bf614, 0x40743e60da479f63,
            0x408296bb918d6f63, 0x40754bf73b88949e, 0x0000000000000001, 0x40882eb006c73652,
            0x0000000000000001, 0x408d487a4f6fc2d6, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ac,
            0x3f50624dd2f1a9fc, 0x3fea53dd9ed22ed7, 0x3f5ae111eb18ed2d, 0xbeb544fd883c9b30,
            0x4089b3d3b983115d, 0x0000000000000001, 0x3fd310d7a9b4d2de, 0x3fdf6c8b4395810d,
            0x4075de59fc7957b6, 0xbfd81ffeb7744a2b, 0xbfd0db9b97c304bb, 0x40834d826cc26d46,
            0x0000000000000001, 0x3fd594c6fc2008d3, 0x3fef6c8b4395810d, 0x40834d8a0217fe95,
            0xc0060619942e599c, 0xbfeda08356e9c989, 0x407e15e4fe449134, 0x0000000000000001,
            0x3fd8bdfb6a67a75d, 0x3ff7916872b01fdd, 0x4089e790d3d5de2d, 0xc01af692146c0d5e,
            0xbffdcf3bd53472aa, 0x40780b5a7f0edaab, 0x0000000000000001, 0x3fdba9b233582fa9,
            0x3fff64e2a47bf614, 0x408f400000000000, 0xc02837e090ad0272, 0xc007e7ec5bd8a3d2,
            0x40743e60da479f63, 0x0000000000000001, 0x3fccbed52ee1b198, 0xd2619817cee01999,
            0x0000000000000015, 0xb189af175b696b30,
        ],
    ),
    (
        "flat/rk4/icao-sea-level/rh50/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc08a6b00758fe8, 0x3fff64805716bcfa, 0x4074294661cf1bb2,
            0x4082700e1a93c7aa, 0x40754bf73b88949e, 0x0000000000000001, 0x40882998b0d07957,
            0x0000000000000001, 0x408d34c502b2e97a, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ac,
            0x3f50624dd2f1a9fc, 0x3fea53e15f4bb1bb, 0x3f5ae115c0abe4b4, 0xbeb51fa51d33e34e,
            0x4089b3db0cdf2bba, 0x0000000000000001, 0x3fd310d51902f22b, 0x3fdf6c8b4395810d,
            0x4075e389b2ce609b, 0xbfd81fb74c6379a9, 0xbfd0c1cbab9c2531, 0x408356a1eb40d93f,
            0x0000000000000001, 0x3fd5900dbfdc4a6f, 0x3fef6c8b4395810d, 0x408352e6e46d24bb,
            0xc0060439faf6497a, 0xbfed85d2999d1a0d, 0x407e17490872d6d2, 0x0000000000000001,
            0x3fd8bd54f4f14ffc, 0x3ff7916872b01fdd, 0x4089eb9e17f94a79, 0xc01aeed5186acb14,
            0xbffdc526a2b3f863, 0x4078007dc7106e9d, 0x0000000000000001, 0x3fdbadf795748f00,
            0x3fff64805716bcfa, 0x408f400000000000, 0xc0282c2f9aa338f2, 0xc007e727c10e31b2,
            0x4074294661cf1bb2, 0x0000000000000001, 0x3fcb694348072022, 0x392f07b4fa65a45e,
            0x0000000000000015, 0x2a32fcc0f491b817,
        ],
    ),
    (
        "flat/rk4/icao-sea-level/rh100/no-zones",
        &[
            0x408f400000000000, 0x3fc0894039dc54ee, 0x3fff58c694d8f3ab, 0x407449a5434d892c,
            0x4082ab728fd22b45, 0x40755334ef1fc968, 0x0000000000000001, 0x408836a145084d18,
            0x0000000000000001, 0x408d52064cd8235f, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007a9,
            0x3f50624dd2f1a9fc, 0x3fea53df04066c75, 0x3f5ae1135825fc8b, 0xbeb53719fa38dcdc,
            0x4089b3d672f8c9a7, 0x0000000000000001, 0x3fd313e549c9f2f9, 0x3fdf5c28f5c28f63,
            0x4075d664f797206c, 0xbfd7ef59489be294, 0xbfd0c1ac771e67c8, 0x408353934da87267,
            0x0000000000000001, 0x3fd5950908e5a8f6, 0x3fef5c28f5c28f63, 0x408348cfc84e1f26,
            0xc005ec9191b55190, 0xbfed7681531b5e41, 0x407e26842564b55a, 0x0000000000000001,
            0x3fd8baff95c1c29e, 0x3ff7851eb851ea9f, 0x4089e3aa9b2260f1, 0xc01adbbb82f96142,
            0xbffda7c8f4f877c8, 0x40781d27e54c0dd4, 0x0000000000000001, 0x3fdba5ddf05ff9e2,
            0x3fff58c694d8f3ab, 0x408f400000000000, 0xc028279222ed7eb8, 0xc007cfb43c929f16,
            0x407449a5434d892c, 0x0000000000000001, 0x3fcd101e5e438a3a, 0x03e0ad065b41d6ca,
            0x0000000000000015, 0x6186f69bfb179d86,
        ],
    ),
    (
        "flat/rk4/icao-sea-level/rh100/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc08a6b00758fe8, 0x3fff64805716bcfa, 0x4074294661cf1bb2,
            0x4082700e1a93c7aa, 0x40755334ef1fc968, 0x0000000000000001, 0x4088230d81736a92,
            0x0000000000000001, 0x408d29dd93b83171, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ac,
            0x3f50624dd2f1a9fc, 0x3fea53e15f4bb1bb, 0x3f5ae115c0abe4b4, 0xbeb51fa51d33e34e,
            0x4089b3db0cdf2bba, 0x0000000000000001, 0x3fd313e3ad9ac1c2, 0x3fdf6c8b4395810d,
            0x4075e389b2ce609b, 0xbfd81fb74c6379a9, 0xbfd0c1cbab9c2531, 0x408356a1eb40d93f,
            0x0000000000000001, 0x3fd59373439ced00, 0x3fef6c8b4395810d, 0x408352e6e46d24bb,
            0xc0060439faf6497a, 0xbfed85d2999d1a0d, 0x407e17490872d6d2, 0x0000000000000001,
            0x3fd8c21c91eb9cba, 0x3ff7916872b01fdd, 0x4089eb9e17f94a79, 0xc01aeed5186acb14,
            0xbffdc526a2b3f863, 0x4078007dc7106e9d, 0x0000000000000001, 0x3fdbb121ff69dca1,
            0x3fff64805716bcfa, 0x408f400000000000, 0xc0282c2f9aa338f2, 0xc007e727c10e31b2,
            0x4074294661cf1bb2, 0x0000000000000001, 0x3fcb09042f3a7ae0, 0x392f07b4fa65a45e,
            0x0000000000000015, 0x2a32fcc0f491b817,
        ],
    ),
    (
        "flat/rk4/1500m-25C-830hPa/rh0/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc0c1e16df8fc77, 0x3ffbe88d5da11641, 0x4078e53bfc2fb538,
            0x408c1d05b9e82ddb, 0x4075a2608de6e31f, 0x0000000000000001, 0x408dc28facc5bc3f,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d2,
            0x3f50624dd2f1a9fc, 0x3fea544fff9bec29, 0x3f5ae186cf1e4b54, 0xbeb0d27d0982f173,
            0x4089b4b30e53a878, 0x0000000000000001, 0x3fd3348d35b3db19, 0x3fdbe76c8b439587,
            0x4074393ee0296704, 0xbfccad387f77f5c2, 0xbfc5ede122629260, 0x4084ecca28328707,
            0x0000000000000001, 0x3fd4fa07bcb26ee8, 0x3febef9db22d0e5c, 0x408270f98d83588b,
            0xc000e69ce4bb66f7, 0xbfe3f6885a551954, 0x40815ae666be2340, 0x0000000000000001,
            0x3fd6fe2a6631fae2, 0x3ff4f1a9fbe76bf5, 0x408962127d3f7be4, 0xc015983d7b8887b4,
            0xbff49c738374b635, 0x407d35627235069b, 0x0000000000000001, 0x3fd9671c32b24bbb,
            0x3ffbe88d5da11641, 0x408f400000000000, 0xc023c8ebda976b49, 0xc000ef41ce22e373,
            0x4078e53bfc2fb538, 0x0000000000000001, 0x3fdb7786eaa42959, 0x6ebd632fe11449e8,
            0x0000000000000015, 0x19fb630b6f496002,
        ],
    ),
    (
        "flat/rk4/1500m-25C-830hPa/rh50/no-zones",
        &[
            0x408f400000000000, 0x3fc0c0a6a6cb2acd, 0x3ffbe89f6f022d9b, 0x4078f39aeb232dc1,
            0x408c3d84051a3b4d, 0x4075b318f42cf8e0, 0x0000000000000001, 0x408db7b65f523169,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d2,
            0x3f50624dd2f1a9fc, 0x3fea544d947445be, 0x3f5ae184565cac97, 0xbeb0ea9002a6cd7e,
            0x4089b4ae555a0c98, 0x0000000000000001, 0x3fd33b66e177a72a, 0x3fdbe76c8b439587,
            0x407436666b31716f, 0xbfccaecaf2aeedf9, 0xbfc60a34154c8251, 0x4084e6ea05929d70,
            0x0000000000000001, 0x3fd502e9d3a603fd, 0x3febef9db22d0e5c, 0x40826d9f4624de0f,
            0xc000e75b62ec8818, 0xbfe40737f0dd37a0, 0x408159f59e789004, 0x0000000000000001,
            0x3fd708c61c220105, 0x3ff4f1a9fbe76bf5, 0x40895f808429e2b0, 0xc0159c3f94877869,
            0xbff4a2d8daaabfcc, 0x407d3d97ccad90ef, 0x0000000000000001, 0x3fd96e15efd0c538,
            0x3ffbe89f6f022d9b, 0x408f400000000000, 0xc023ced39cceff8d, 0xc000ef65f0e51231,
            0x4078f39aeb232dc1, 0x0000000000000001, 0x3fdb79ab4d3d818c, 0x0eccd0920bdcc7af,
            0x0000000000000015, 0x9317f593ab2a2458,
        ],
    ),
    (
        "flat/rk4/1500m-25C-830hPa/rh50/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc0c1e16df8fc77, 0x3ffbe88d5da11641, 0x4078e53bfc2fb538,
            0x408c1d05b9e82ddb, 0x4075b318f42cf8e0, 0x0000000000000001, 0x408da7f3de6276bd,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d2,
            0x3f50624dd2f1a9fc, 0x3fea544fff9bec29, 0x3f5ae186cf1e4b54, 0xbeb0d27d0982f173,
            0x4089b4b30e53a878, 0x0000000000000001, 0x3fd33b6540371536, 0x3fdbe76c8b439587,
            0x4074393ee0296704, 0xbfccad387f77f5c2, 0xbfc5ede122629260, 0x4084ecca28328707,
            0x0000000000000001, 0x3fd5008915f13a30, 0x3febef9db22d0e5c, 0x408270f98d83588b,
            0xc000e69ce4bb66f7, 0xbfe3f6885a551954, 0x40815ae666be2340, 0x0000000000000001,
            0x3fd70812f456fdd9, 0x3ff4f1a9fbe76bf5, 0x408962127d3f7be4, 0xc015983d7b8887b4,
            0xbff49c738374b635, 0x407d35627235069b, 0x0000000000000001, 0x3fd9720ecd5389b9,
            0x3ffbe88d5da11641, 0x408f400000000000, 0xc023c8ebda976b49, 0xc000ef41ce22e373,
            0x4078e53bfc2fb538, 0x0000000000000001, 0x3fdb7ff0f1515707, 0x6ebd632fe11449e8,
            0x0000000000000015, 0x19fb630b6f496002,
        ],
    ),
    (
        "flat/rk4/1500m-25C-830hPa/rh100/no-zones",
        &[
            0x408f400000000000, 0x3fc0c1f233818054, 0x3ffbd68644ba56d5, 0x407913ede7c90340,
            0x408c86debcc9bcdf, 0x4075c3b860e9c586, 0x0000000000000001, 0x408dc7dfaf742037,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006cd,
            0x3f50624dd2f1a9fc, 0x3fea545020a2eeef, 0x3f5ae186f0df0b4e, 0xbeb0d134499fc422,
            0x4089b4b34ecfa914, 0x0000000000000001, 0x3fd342304752b474, 0x3fdbd70a3d70a3dd,
            0x40742eb175326d03, 0xbfcc5a54eb741991, 0xbfc5d3d833fd5c50, 0x4084ef88c91a87b4,
            0x0000000000000001, 0x3fd505e8a1d1029c, 0x3febd70a3d70a3dd, 0x4082653ce0544dd1,
            0xc000c44380e69407, 0xbfe3cea49311eb8f, 0x408168a7688265ce, 0x0000000000000001,
            0x3fd707bbcdf7f25e, 0x3ff4e147ae147a4d, 0x408959d952a179ef, 0xc0157b23f694353f,
            0xbff46f60ec8756b5, 0x407d5e509a606a41, 0x0000000000000001, 0x3fd9691776f39575,
            0x3ffbd68644ba56d5, 0x408f400000000000, 0xc023b8e5899dc759, 0xc000cb339c556496,
            0x407913ede7c90340, 0x0000000000000001, 0x3fdb73de208aead7, 0xced996a228229a14,
            0x0000000000000015, 0x85544df478312dc8,
        ],
    ),
    (
        "flat/rk4/1500m-25C-830hPa/rh100/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc0c1e16df8fc77, 0x3ffbe88d5da11641, 0x4078e53bfc2fb538,
            0x408c1d05b9e82ddb, 0x4075c3b860e9c586, 0x0000000000000001, 0x408d89ed0118b7fe,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d2,
            0x3f50624dd2f1a9fc, 0x3fea544fff9bec29, 0x3f5ae186cf1e4b54, 0xbeb0d27d0982f173,
            0x4089b4b30e53a878, 0x0000000000000001, 0x3fd342305da2951e, 0x3fdbe76c8b439587,
            0x4074393ee0296704, 0xbfccad387f77f5c2, 0xbfc5ede122629260, 0x4084ecca28328707,
            0x0000000000000001, 0x3fd507065d63d594, 0x3febef9db22d0e5c, 0x408270f98d83588b,
            0xc000e69ce4bb66f7, 0xbfe3f6885a551954, 0x40815ae666be2340, 0x0000000000000001,
            0x3fd7120058b44e2b, 0x3ff4f1a9fbe76bf5, 0x408962127d3f7be4, 0xc015983d7b8887b4,
            0xbff49c738374b635, 0x407d35627235069b, 0x0000000000000001, 0x3fd97ccf308d16a1,
            0x3ffbe88d5da11641, 0x408f400000000000, 0xc023c8ebda976b49, 0xc000ef41ce22e373,
            0x4078e53bfc2fb538, 0x0000000000000001, 0x3fdb881bff35d219, 0x6ebd632fe11449e8,
            0x0000000000000015, 0x19fb630b6f496002,
        ],
    ),
    (
        "flat/rk45/icao-sea-level/rh0/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc089ab613b8bb0, 0x3fff6485bd0a65cb, 0x4074294427e2af95,
            0x4082700a082c38c3, 0x407544b4bbbd89be, 0x0000000000000001, 0x40884d0e66f23c65,
            0x0000000000000001, 0x408d4a2b40fa3946, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000c9,
            0x3f50624dd2f1a9fc, 0x3fea53e15f4bb11a, 0x3f5ae115c0abe47b, 0xbeb51fa51d3a2461,
            0x4089b3db0cdf2bb9, 0x0000000000000001, 0x3fd30dc2baedb49b, 0x3fde33269df97ab0,
            0x4075252c697cfd04, 0xbfd4963b3290ce0e, 0xbfcf1f4ddc26a2aa, 0x408389c145af20b9,
            0x0000000000000001, 0x3fd572cfdb3b372d, 0x3fef19934efcbd5c, 0x40832bcb84ba52e5,
            0xc0057c83a440701d, 0xbfecfc949f28c9a3, 0x407e3d4715481288, 0x0000000000000001,
            0x3fd8a6d0cc6f26d2, 0x3ff78cc9a77e5eb0, 0x4089e82474e389e3, 0xc01ae3c8a85b4630,
            0xbffdbb5149db2ad2, 0x4078037b4595251b, 0x0000000000000001, 0x3fdba991452e6ef5,
            0x3fff6485bd0a65cb, 0x408f400000000000, 0xc0282c37f6bff4ad, 0xc007e7328cf57f98,
            0x4074294427e2af95, 0x0000000000000001, 0x3fcbd0b29c361719, 0x6aab593f2039841d,
            0x0000000000000015, 0x7e334f66ed409f2e,
        ],
    ),
    (
        "flat/rk45/icao-sea-level/rh50/no-zones",
        &[
            0x408f400000000000, 0x3fc087c778eabfb6, 0x3fff64e2e1c2dcff, 0x40743e6192986670,
            0x408296bce40c8dbb, 0x40754bf73b88949e, 0x0000000000000001, 0x408848624ca5f016,
            0x0000000000000001, 0x408d48210cecf44b, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000c9,
            0x3f50624dd2f1a9fc, 0x3fea53dd9ed22e2b, 0x3f5ae111eb18ecdb, 0xbeb544fd88434b0a,
            0x4089b3d3b983115d, 0x0000000000000001, 0x3fd310d7a9b4d2de, 0x3fde33269df97ab0,
            0x40752055055616e5, 0xbfd496dd27b0cd5e, 0xbfcf4f7e922bf8e0, 0x408380d7da7444af,
            0x0000000000000001, 0x3fd57a977e9e2116, 0x3fef19934efcbd5c, 0x40832671fc3ad60b,
            0xc0057e3f48d3c3bd, 0xbfed1734af14db96, 0x407e3b91582045db, 0x0000000000000001,
            0x3fd8ac68f0410f16, 0x3ff78cc9a77e5eb0, 0x4089e4182a93eaaa, 0xc01aeb8192b2b348,
            0xbffdc5640eac135a, 0x40780e56c6143b99, 0x0000000000000001, 0x3fdba8830e2dee1a,
            0x3fff64e2e1c2dcff, 0x408f400000000000, 0xc02837e21c7ddf58, 0xc007e7ecd6666e10,
            0x40743e6192986670, 0x0000000000000001, 0x3fccbee233b6df8e, 0x8c695b77715f8bcb,
            0x0000000000000015, 0xc8abf4335c4c3468,
        ],
    ),
    (
        "flat/rk45/icao-sea-level/rh50/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc089ab613b8bb0, 0x3fff6485bd0a65cb, 0x4074294427e2af95,
            0x4082700a082c38c3, 0x40754bf73b88949e, 0x0000000000000001, 0x40882c705b408fc6,
            0x0000000000000001, 0x408d4a2b40fa3946, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000c9,
            0x3f50624dd2f1a9fc, 0x3fea53e15f4bb11a, 0x3f5ae115c0abe47b, 0xbeb51fa51d3a2461,
            0x4089b3db0cdf2bb9, 0x0000000000000001, 0x3fd310d51902f22c, 0x3fde33269df97ab0,
            0x4075252c697cfd04, 0xbfd4963b3290ce0e, 0xbfcf1f4ddc26a2aa, 0x408389c145af20b9,
            0x0000000000000001, 0x3fd5762236264600, 0x3fef19934efcbd5c, 0x40832bcb84ba52e5,
            0xc0057c83a440701d, 0xbfecfc949f28c9a3, 0x407e3d4715481288, 0x0000000000000001,
            0x3fd8ab9d467b10b3, 0x3ff78cc9a77e5eb0, 0x4089e82474e389e3, 0xc01ae3c8a85b4630,
            0xbffdbb5149db2ad2, 0x4078037b4595251b, 0x0000000000000001, 0x3fdbaccc1867aad4,
            0x3fff6485bd0a65cb, 0x408f400000000000, 0xc0282c37f6bff4ad, 0xc007e7328cf57f98,
            0x4074294427e2af95, 0x0000000000000001, 0x3fcb6922e8be823a, 0x6aab593f2039841d,
            0x0000000000000015, 0x67ae91e9d0fb7ec2,
        ],
    ),
    (
        "flat/rk45/icao-sea-level/rh100/no-zones",
        &[
            0x408f400000000000, 0x3fc0887b94722144, 0x3fff58c889bd1602, 0x407449ac06aa6ffa,
            0x4082ab7f026636ff, 0x40755334ef1fc968, 0x0000000000000001, 0x40884d1738c80c12,
            0x0000000000000001, 0x408d6a2fa538c985, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000c8,
            0x3f50624dd2f1a9fc, 0x3fea53df04066bcf, 0x3f5ae1135825fc45, 0xbeb53719fa3f5b50,
            0x4089b3d672f8c9a7, 0x0000000000000001, 0x3fd313e549c9f2f9, 0x3fde33269df97ab0,
            0x40752223775d70c8, 0xbfd496a296bce3e3, 0xbfcf3d8375b61eea, 0x4083842bd73f5e41,
            0x0000000000000001, 0x3fd57c3f769aea72, 0x3fef19934efcbd5c, 0x40832962a2484be1,
            0xc0057f9cc1488d91, 0xbfed0892cf1e4178, 0x407e44bfba0bf65e, 0x0000000000000001,
            0x3fd8aceb94e2f0f1, 0x3ff78cc9a77e5eb0, 0x4089e970b8315da3, 0xc01aee168270ee20,
            0xbffdb81677061637, 0x4078183242923ef1, 0x0000000000000001, 0x3fdba7d862886433,
            0x3fff58c889bd1602, 0x408f400000000000, 0xc028279eb4f172d8, 0xc007cfb8265ae002,
            0x407449ac06aa6ffa, 0x0000000000000001, 0x3fcd109efaebfc76, 0xc5e5fce28fea924d,
            0x0000000000000015, 0xe97af33d26e1013a,
        ],
    ),
    (
        "flat/rk45/icao-sea-level/rh100/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc089ab613b8bb0, 0x3fff6485bd0a65cb, 0x4074294427e2af95,
            0x4082700a082c38c3, 0x40755334ef1fc968, 0x0000000000000001, 0x40882c705b408fc6,
            0x0000000000000001, 0x408d2ef56bb5537d, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000c9,
            0x3f50624dd2f1a9fc, 0x3fea53e15f4bb11a, 0x3f5ae115c0abe47b, 0xbeb51fa51d3a2461,
            0x4089b3db0cdf2bb9, 0x0000000000000001, 0x3fd313e3ad9ac1c3, 0x3fde33269df97ab0,
            0x4075252c697cfd04, 0xbfd4963b3290ce0e, 0xbfcf1f4ddc26a2aa, 0x408389c145af20b9,
            0x0000000000000001, 0x3fd5797343b8b216, 0x3fef19934efcbd5c, 0x40832bcb84ba52e5,
            0xc0057c83a440701d, 0xbfecfc949f28c9a3, 0x407e3d4715481288, 0x0000000000000001,
            0x3fd8b064ec94d8ba, 0x3ff78cc9a77e5eb0, 0x4089e82474e389e3, 0xc01ae3c8a85b4630,
            0xbffdbb5149db2ad2, 0x4078037b4595251b, 0x0000000000000001, 0x3fdbaffa0d188736,
            0x3fff6485bd0a65cb, 0x408f400000000000, 0xc0282c37f6bff4ad, 0xc007e7328cf57f98,
            0x4074294427e2af95, 0x0000000000000001, 0x3fcb08e5f9bd7466, 0x6aab593f2039841d,
            0x0000000000000015, 0x67ae91e9d0fb7ec2,
        ],
    ),
    (
        "flat/rk45/1500m-25C-830hPa/rh0/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc0c1d3bc40b5eb, 0x3ffbe894029742ec, 0x4078e5383dfb8bae,
            0x408c1cfd45cabd8f, 0x4075a2608de6e31f, 0x0000000000000001, 0x408dd2b80b32eb3d,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000b3,
            0x3f50624dd2f1a9fc, 0x3fea544fff9bebef, 0x3f5ae186cf1e4b5b, 0xbeb0d27d09852cb2,
            0x4089b4b30e53a878, 0x0000000000000001, 0x3fd3348d35b3db19, 0x3fda5c1c6088d6d8,
            0x4073355a4b3c6c39, 0xbfc52e415cf1c88b, 0xbfc3ae4ef67262fa, 0x408526e0788c8e39,
            0x0000000000000001, 0x3fd4e2e0f272f3b6, 0x3feb947496aad1d6, 0x40823f6bb2029c93,
            0xc0005e7eff55bf54, 0xbfe3808477a42bb7, 0x40816f438d0b521b, 0x0000000000000001,
            0x3fd6ef4fe5207cd0, 0x3ff4fd6d7e889c20, 0x40896ccb6e8922bc, 0xc015b2b3d7c79867,
            0xbff4b0d2deae05ed, 0x407d2d2da1372e6e, 0x0000000000000001, 0x3fd96b1c5b1c1fb2,
            0x3ffbe894029742ec, 0x408f400000000000, 0xc023c8fb765006df, 0xc000ef4f180f39de,
            0x4078e5383dfb8bae, 0x0000000000000001, 0x3fdb778891c09c3b, 0xd3da8f0dcdbaeb76,
            0x0000000000000015, 0x2d8559414f02b2a0,
        ],
    ),
    (
        "flat/rk45/1500m-25C-830hPa/rh50/no-zones",
        &[
            0x408f400000000000, 0x3fc0c0967b5ab7d9, 0x3ffbe8a2c2c91292, 0x4078f39c03818162,
            0x408c3d867fbf3213, 0x4075b318f42cf8e0, 0x0000000000000001, 0x408dd1e663e3f17c,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000b3,
            0x3f50624dd2f1a9fc, 0x3fea544d94744582, 0x3f5ae184565cac97, 0xbeb0ea9002a9249e,
            0x4089b4ae555a0c98, 0x0000000000000001, 0x3fd33b66e177a72a, 0x3fda5c1c6088d6d8,
            0x407332c90129cdb3, 0xbfc5302ce53dccfc, 0xbfc3c7dd8540858f, 0x4085213b72ca3339,
            0x0000000000000001, 0x3fd4eb942f962533, 0x3feb947496aad1d6, 0x40823c152d5f7183,
            0xc0005f27cf5b2387, 0xbfe3912156bb1c65, 0x40816e2a82be0cdc, 0x0000000000000001,
            0x3fd6f9e1a24cc1c7, 0x3ff4fd6d7e889c20, 0x40896a3e5052c1da, 0xc015b6c092c58bdc,
            0xbff4b72c210d2ace, 0x407d3573feec8072, 0x0000000000000001, 0x3fd9720651c50aa7,
            0x3ffbe8a2c2c91292, 0x408f400000000000, 0xc023cedf46405d17, 0xc000ef6c9872d92c,
            0x4078f39c03818162, 0x0000000000000001, 0x3fdb79aad2102810, 0xd0b845658856206f,
            0x0000000000000015, 0x6773f99007fa4b86,
        ],
    ),
    (
        "flat/rk45/1500m-25C-830hPa/rh50/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc0c1d3bc40b5eb, 0x3ffbe894029742ec, 0x4078e5383dfb8bae,
            0x408c1cfd45cabd8f, 0x4075b318f42cf8e0, 0x0000000000000001, 0x408db18131aa27ae,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000b3,
            0x3f50624dd2f1a9fc, 0x3fea544fff9bebef, 0x3f5ae186cf1e4b5b, 0xbeb0d27d09852cb2,
            0x4089b4b30e53a878, 0x0000000000000001, 0x3fd33b6540371536, 0x3fda5c1c6088d6d8,
            0x4073355a4b3c6c39, 0xbfc52e415cf1c88b, 0xbfc3ae4ef67262fa, 0x408526e0788c8e39,
            0x0000000000000001, 0x3fd4e95683bce2c5, 0x3feb947496aad1d6, 0x40823f6bb2029c93,
            0xc0005e7eff55bf54, 0xbfe3808477a42bb7, 0x40816f438d0b521b, 0x0000000000000001,
            0x3fd6f91493552b4d, 0x3ff4fd6d7e889c20, 0x40896ccb6e8922bc, 0xc015b2b3d7c79867,
            0xbff4b0d2deae05ed, 0x407d2d2da1372e6e, 0x0000000000000001, 0x3fd9760518f193f3,
            0x3ffbe894029742ec, 0x408f400000000000, 0xc023c8fb765006df, 0xc000ef4f180f39de,
            0x4078e5383dfb8bae, 0x0000000000000001, 0x3fdb7ff290b78c84, 0xd3da8f0dcdbaeb76,
            0x0000000000000015, 0x2d8559414f02b2a0,
        ],
    ),
    (
        "flat/rk45/1500m-25C-830hPa/rh100/no-zones",
        &[
            0x408f400000000000, 0x3fc0c1e4a3d6d9d3, 0x3ffbd6892bef2cd2, 0x407913eedb3f7dc9,
            0x408c86e0e6ae68ff, 0x4075c3b860e9c586, 0x0000000000000001, 0x408ddecea66b906c,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000b2,
            0x3f50624dd2f1a9fc, 0x3fea545020a2eeb6, 0x3f5ae186f0df0b5a, 0xbeb0d13449a1f3b7,
            0x4089b4b34ecfa914, 0x0000000000000001, 0x3fd342304752b474, 0x3fda5c1c6088d6d8,
            0x40733580010ccb4b, 0xbfc52e2c882ca3f9, 0xbfc3acd799db325a, 0x40852737fa7ce649,
            0x0000000000000001, 0x3fd4ef9e6a019591, 0x3feb947496aad1d6, 0x408240f4da6ea061,
            0xc00060e5ec1cd6b7, 0xbfe378dfc0205c28, 0x408177634059ca55, 0x0000000000000001,
            0x3fd6fce4bc689043, 0x3ff4d477bbf93ff7, 0x40894e15e691cca1, 0xc0155e63c3fab967,
            0xbff45966f33673b2, 0x407d673378b26328, 0x0000000000000001, 0x3fd964c6f79483e3,
            0x3ffbd6892bef2cd2, 0x408f400000000000, 0xc023b8efbcbebee4, 0xc000cb396abf0da4,
            0x407913eedb3f7dc9, 0x0000000000000001, 0x3fdb73ddb564ae6c, 0x33e6aa6c8588f471,
            0x0000000000000015, 0xcdccca1fe7eb995b,
        ],
    ),
    (
        "flat/rk45/1500m-25C-830hPa/rh100/zones-rh20-rh80",
        &[
            0x408f400000000000, 0x3fc0c1d3bc40b5eb, 0x3ffbe894029742ec, 0x4078e5383dfb8bae,
            0x408c1cfd45cabd8f, 0x4075c3b860e9c586, 0x0000000000000001, 0x408d902b929f2ec3,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000000b3,
            0x3f50624dd2f1a9fc, 0x3fea544fff9bebef, 0x3f5ae186cf1e4b5b, 0xbeb0d27d09852cb2,
            0x4089b4b30e53a878, 0x0000000000000001, 0x3fd342305da2951e, 0x3fda5c1c6088d6d8,
            0x4073355a4b3c6c39, 0xbfc52e415cf1c88b, 0xbfc3ae4ef67262fa, 0x408526e0788c8e39,
            0x0000000000000001, 0x3fd4efc11978188e, 0x3feb947496aad1d6, 0x40823f6bb2029c93,
            0xc0005e7eff55bf54, 0xbfe3808477a42bb7, 0x40816f438d0b521b, 0x0000000000000001,
            0x3fd702d9f8577610, 0x3ff4fd6d7e889c20, 0x40896ccb6e8922bc, 0xc015b2b3d7c79867,
            0xbff4b0d2deae05ed, 0x407d2d2da1372e6e, 0x0000000000000001, 0x3fd980bd16c64856,
            0x3ffbe894029742ec, 0x408f400000000000, 0xc023c8fb765006df, 0xc000ef4f180f39de,
            0x4078e5383dfb8bae, 0x0000000000000001, 0x3fdb881d95d0ae86, 0xd3da8f0dcdbaeb76,
            0x0000000000000015, 0x2d8559414f02b2a0,
        ],
    ),
];

/// `TrajectorySolver::solve`, the 30-degree lob to ground impact, Euler, RK4 and RK45 at both
/// stations, RH 0, no zones.
///
/// IMMUTABLE: nothing in 0.50.0 may change it, the aloft vapour model included.
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
/// MAY CHANGE ONLY WITH THE ALOFT VAPOUR MODEL: byte-identical through every other change in
/// 0.50.0.
///
/// Re-pinned by the aloft vapour model (0.50.0): the unified aloft vapour model holds the
/// governing air's (the station's or the zone's) water-vapour mole fraction for the lapse, the
/// density and the sound speed, with a virtual-temperature lapse. 589 of 954 values moved, all in
/// cases with a humidity above 0 and all off the station's own height; every RH-0 and
/// station-height value is unchanged. The largest float shift of each moved case,
/// pinned before the aloft vapour model -> pinned now:
/// - `lob30/euler/icao-sea-level/rh50/no-zones` mach_0_9_distance_m: 965.3861958282084 -> 965.6329488435597 (+2.56e-04)
/// - `lob30/euler/icao-sea-level/rh100/no-zones` points[22446].speed: 127.5594013596002 -> 127.52526950945204 (-2.68e-04)
/// - `lob30/euler/icao-sea-level/rh50/zones-rh20-rh80` mach_0_9_distance_m: 957.0403448161148 -> 957.288978457813 (+2.60e-04)
/// - `lob30/euler/1500m-25C-830hPa/rh50/no-zones` mach_1_2_distance_m: 830.2790841155103 -> 829.8635264947768 (-5.01e-04)
/// - `lob30/euler/1500m-25C-830hPa/rh100/no-zones` impact_energy: 138.74260738670813 -> 138.6534135669379 (-6.43e-04)
/// - `lob30/euler/1500m-25C-830hPa/rh50/zones-rh20-rh80` impact_energy: 136.72679730916047 -> 136.6695828704991 (-4.18e-04)
/// - `lob30/rk4/icao-sea-level/rh50/no-zones` mach_1_0_distance_m: 820.1693390213044 -> 819.8409845907621 (-4.00e-04)
/// - `lob30/rk4/icao-sea-level/rh100/no-zones` mach_1_2_distance_m: 676.2747456535427 -> 675.8844782945004 (-5.77e-04)
/// - `lob30/rk4/icao-sea-level/rh50/zones-rh20-rh80` points[14911].speed: 131.037273180594 -> 131.01333338506552 (-1.83e-04)
/// - `lob30/rk4/1500m-25C-830hPa/rh50/no-zones` impact_energy: 137.85552498537697 -> 137.81154304391646 (-3.19e-04)
/// - `lob30/rk4/1500m-25C-830hPa/rh100/no-zones` mach_1_0_distance_m: 1017.9663353405539 -> 1017.1698100019853 (-7.82e-04)
/// - `lob30/rk4/1500m-25C-830hPa/rh50/zones-rh20-rh80` mach_1_2_distance_m: 829.7986799955099 -> 829.3842364257149 (-4.99e-04)
/// - `lob30/rk45/icao-sea-level/rh50/no-zones` points[1497].speed: 131.82263824903768 -> 131.80415102367417 (-1.40e-04)
/// - `lob30/rk45/icao-sea-level/rh100/no-zones` points[1499].speed: 132.03292171636494 -> 131.99578090414283 (-2.81e-04)
/// - `lob30/rk45/icao-sea-level/rh50/zones-rh20-rh80` points[1493].speed: 131.06719266451245 -> 131.04325583766737 (-1.83e-04)
/// - `lob30/rk45/1500m-25C-830hPa/rh50/no-zones` mach_0_9_distance_m: 1197.4854956343006 -> 1200.1874036118031 (+2.26e-03)
/// - `lob30/rk45/1500m-25C-830hPa/rh100/no-zones` mach_1_2_distance_m: 835.9277836191676 -> 832.1649155434807 (-4.50e-03)
/// - `lob30/rk45/1500m-25C-830hPa/rh50/zones-rh20-rh80` mach_1_0_distance_m: 1014.7578097506399 -> 1011.5831210537965 (-3.13e-03)
///
/// glibc rounds a few of these a last bit apart; its bits are in `LIBM_ALTERNATES`.
#[rustfmt::skip]
const TS_LOB_HUMID: Pinned = &[
    (
        "lob30/euler/icao-sea-level/rh50/no-zones",
        &[
            0x40b25fc50ae59127, 0x4091a93883ec2b36, 0x403de541b8e36c38, 0x40619720c38db839,
            0x405c1212d2ed542e, 0x40754bf73b88949e, 0x0000000000000001, 0x408518fa0c5abf29,
            0x0000000000000001, 0x40899c11d4fd4b9c, 0x0000000000000001, 0x408e2d10477bb90a,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x00000000000074c9,
            0x3f50624dd2f1a9fc, 0x3fe6cb0e0b537889, 0x3fda51948f9fd0d9, 0xbec546bfb16f49c3,
            0x4089b3c9348be756, 0x0000000000000001, 0x3fd310db58b42b79, 0x401de56041893af3,
            0x40a06f27c38640c1, 0x408fa900bf3a01da, 0xc03217246f60a8c1, 0x4067cb009a40b264,
            0x0000000000000001, 0x3fc09339c63c50d4, 0x402de5604189310a, 0x40a95bb94e1ab7c4,
            0x409117ccf5e89314, 0xc044c9b15e32c29f, 0x406078af195a92d4, 0x0000000000000001,
            0x3fc093e727942290, 0x40366c083126ee55, 0x40afe8504591feb0, 0x40852d131e5de2e9,
            0xc050b117f0532430, 0x405fd624d5c78d41, 0x0000000000000001, 0x3fc095030fe22004,
            0x403de541b8e36c38, 0x40b25fc50ae59127, 0xc059000000000000, 0xc0574bc45bebcaff,
            0x40619720c38db839, 0x0000000000000001, 0x3fc092c74984e08a, 0x40d1f49cf81486af,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/euler/icao-sea-level/rh100/no-zones",
        &[
            0x40b26b552d0863ed, 0x4091b316a8e7eea0, 0x403ded89b82f459f, 0x40619d58035df21a,
            0x405c25ecf2891825, 0x40755334ef1fc968, 0x0000000000000001, 0x40851cf3e5a725d1,
            0x0000000000000001, 0x4089a4446c533a99, 0x0000000000000001, 0x408e36e3fafc287f,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x00000000000074e9,
            0x3f50624dd2f1a9fc, 0x3fe6cb107662f65b, 0x3fda51975a742b82, 0xbec538da0c5e0036,
            0x4089b3cbee9f472f, 0x0000000000000001, 0x3fd313e8f7ee0281, 0x401ded916872b3cb,
            0x40a0780134efb429, 0x408fb991d22a42b5, 0xc03218fa317f1b41, 0x4067d4243d49eccc,
            0x0000000000000001, 0x3fc0933a615882fd, 0x402ded916872a9da, 0x40a96a90cccb6abb,
            0x409121a68c2a4a19, 0xc044cc8dd361b5eb, 0x40607f3dacb9df6f, 0x0000000000000001,
            0x3fc093e56ccaadb7, 0x4036722d0e5608fd, 0x40affbc2dd6dda6a, 0x40853abfc36416d6,
            0xc050b3beb3e1fd53, 0x405fe19e04012bae, 0x0000000000000001, 0x3fc095026305ddb8,
            0x403ded89b82f459f, 0x40b26b552d0863ed, 0xc059000000000000, 0xc0574fe5246103c5,
            0x40619d58035df21a, 0x0000000000000001, 0x3fc092c7316c9e37, 0xc7b165403ffd5e4d,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/euler/icao-sea-level/rh50/zones-rh20-rh80",
        &[
            0x40b23e54491cd4cf, 0x409192c8fa8117a8, 0x403dd070b15628c3, 0x4061828366a463cc,
            0x405bd06e300ec672, 0x40754bf73b88949e, 0x0000000000000001, 0x40851c2fcc099912,
            0x0000000000000001, 0x4089926a8c76483a, 0x0000000000000001, 0x408dea4fd3f00c71,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000007478,
            0x3f50624dd2f1a9fc, 0x3fe6cb148bee4258, 0x3fda519c11be810c, 0xbec5216177e29fc7,
            0x4089b3d0899ddd90, 0x0000000000000001, 0x3fd310d8c768a2b1, 0x401dd0e560418cd7,
            0x40a05b9b532bc589, 0x408f87aa54242115, 0xc03210b48599008b, 0x4067aad797fa1644,
            0x0000000000000001, 0x3fc0932a74188368, 0x402dd0e560418302, 0x40a93520b54c4d02,
            0x40910079cc1c9881, 0xc044c313610e0709, 0x406060591f2f7f4a, 0x0000000000000001,
            0x3fc09414760ce465, 0x40365cac08312bb1, 0x40afb2563e303dae, 0x40850a1f101460f1,
            0xc050ab6d57eb8ea4, 0x405fae6386cbddeb, 0x0000000000000001, 0x3fc0952c5b3ec9c6,
            0x403dd070b15628c3, 0x40b23e54491cd4cf, 0xc059000000000000, 0xc05742f629dbc3ee,
            0x4061828366a463cc, 0x0000000000000001, 0x3fc092d0368cefdc, 0xdbfba846b6525072,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/euler/1500m-25C-830hPa/rh50/no-zones",
        &[
            0x40b5a866e6381543, 0x4094540638e64e03, 0x403ff707eef5ae31, 0x40637c88ffbc6b40,
            0x4061396581ddc680, 0x4075b318f42cf8e0, 0x0000000000000001, 0x4089eee88094325e,
            0x0000000000000001, 0x408faf5ec243504f, 0x0000000000000001, 0x4092b44504723019,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000007cde,
            0x3f50624dd2f1a9fc, 0x3fe6cbd01341ba63, 0x3fda52749bd4531e, 0xbec0ebaecef0e832,
            0x4089b4a40076fcf5, 0x0000000000000001, 0x3fd33b6a724ce681, 0x401ff6c8b4395c7a,
            0x40a2cbce46fa24d5, 0x4092140e38e14f67, 0xc03276323f69dd51, 0x406a0d55b116fc80,
            0x0000000000000001, 0x3fc096b7e388392a, 0x402ff74bc6a7e819, 0x40ad61f966d9f159,
            0x4093c2109f9bde3c, 0xc0456991dc612f37, 0x40626e61f148c075, 0x0000000000000001,
            0x3fc092ad00e7ff9c, 0x4037f958106254c5, 0x40b2a94b9a0cf922, 0x4088e2301d673e22,
            0xc05145720e250932, 0x4061b82dd536c94f, 0x0000000000000001, 0x3fc092dfd5a910ca,
            0x403ff707eef5ae31, 0x40b5a866e6381543, 0xc059000000000000, 0xc0582f9ad6b01254,
            0x40637c88ffbc6b40, 0x0000000000000001, 0x3fc09315ec0edfb2, 0x859b002705d0a2ef,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/euler/1500m-25C-830hPa/rh100/no-zones",
        &[
            0x40b5c6c16800f8b9, 0x40946d099b7a7390, 0x40400561b7f0448b, 0x40638c12f81d2a9a,
            0x406154e8c391985a, 0x4075c3b860e9c586, 0x0000000000000001, 0x4089fcb2898e285a,
            0x0000000000000001, 0x408fc6c2911dcc58, 0x0000000000000001, 0x4092c4763e60f20f,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000007d2c,
            0x3f50624dd2f1a9fc, 0x3fe6cbd47d77bbe3, 0x3fda5279b4e274ff, 0xbec0d24fa9443d64,
            0x4089b4a8faeec90c, 0x0000000000000001, 0x3fd34233da20e3ed, 0x402005a1cac08545,
            0x40a2e27f35eac1ab, 0x40922905b3df879b, 0xc0327a603f6e9fea, 0x406a231793952f4c,
            0x0000000000000001, 0x3fc096bdd118cea0, 0x403005a1cac07f55, 0x40ad885d0c48d405,
            0x4093dafd63f0a959, 0xc0456fbfcbe15f99, 0x40627eb26c92cb39, 0x0000000000000001,
            0x3fc092acfd02e553, 0x40380872b020cba2, 0x40b2c2c84d2928cc, 0x408903f0c4f57453,
            0xc0514b6526e43f8c, 0x4061c6b26de9409b, 0x0000000000000001, 0x3fc092df51fd6abf,
            0x40400561b7f0448b, 0x40b5c6c16800f8b9, 0xc059000000000000, 0xc058386f5939369a,
            0x40638c12f81d2a9a, 0x0000000000000001, 0x3fc093166706bae6, 0x8f7d37f92b7232fa,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/euler/1500m-25C-830hPa/rh50/zones-rh20-rh80",
        &[
            0x40b5851036e1d6fa, 0x40943c54c8ea4071, 0x403fe2a8bf191178, 0x406368257931e4b4,
            0x4061156d390e582d, 0x4075b318f42cf8e0, 0x0000000000000001, 0x4089e90239afb859,
            0x0000000000000001, 0x408f9951e9442148, 0x0000000000000001, 0x40928ac99f564e7e,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000007c8f,
            0x3f50624dd2f1a9fc, 0x3fe6cbd444393069, 0x3fda527972c8d8a1, 0xbec0d398a2a9eea0,
            0x4089b4a8ba61c664, 0x0000000000000001, 0x3fd33b68d0b81e39, 0x401fe24dd2f1ae5e,
            0x40a2b7074e3b3c24, 0x4092023c9fa020e1, 0xc03270a447e8f515, 0x4069ef46c1afc841,
            0x0000000000000001, 0x3fc09657730bab1f, 0x402fe2d0e5603a11, 0x40ad39449e520b57,
            0x4093a9b31679163a, 0xc04563b5e51e43f2, 0x406256c18963eaa1, 0x0000000000000001,
            0x3fc092ae114cf994, 0x4037ea3d70a3dde8, 0x40b28cde44018f49, 0x4088bddeaf8a9a3b,
            0xc051408d49f57403, 0x4061a45127739354, 0x0000000000000001, 0x3fc092eb4bee64c9,
            0x403fe2a8bf191178, 0x40b5851036e1d6fa, 0xc059000000000000, 0xc05827ed0527e373,
            0x406368257931e4b4, 0x0000000000000001, 0x3fc093067dba96bd, 0x2cf50f2cfd2c3410,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk4/icao-sea-level/rh50/no-zones",
        &[
            0x40b2612e490bf54c, 0x4091ab336f7e8151, 0x403de6dd09c7ac69, 0x4061977687c0d8b2,
            0x405c132490c3a621, 0x40754bf73b88949e, 0x0000000000000001, 0x40851b1ab2ce9023,
            0x0000000000000001, 0x40899eba56210e1b, 0x0000000000000001, 0x408e324067883af6,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x00000000000074cf,
            0x3f50624dd2f1a9fc, 0x3fe6cce8183899fc, 0x3fda53cc84838c77, 0xbeb544e8dd30cfb1,
            0x4089b3c9bfecfcf4, 0x0000000000000001, 0x3fd310db27e5b0a1, 0x401de66666666a0e,
            0x40a070b03e137376, 0x408fac88dda3d045, 0xc0321710a3c86226, 0x4067cb792c27d1e0,
            0x0000000000000001, 0x3fc0933a0bb9d5c2, 0x402de6e978d4f7b1, 0x40a95e01ac44b961,
            0x409119b09ffd7100, 0xc044ca68b1273db1, 0x406078c3a43a827b, 0x0000000000000001,
            0x3fc093e702587c76, 0x40366d4fdf3b6938, 0x40afeb0694c37f80, 0x40852f2a8f444aa2,
            0xc050b1e2edb303d3, 0x405fd6ae7dd9f273, 0x0000000000000001, 0x3fc09502806b4e44,
            0x403de6dd09c7ac69, 0x40b2612e490bf54c, 0xc059000000000000, 0xc0574cddecb698ce,
            0x4061977687c0d8b2, 0x0000000000000001, 0x3fc092c72722a003, 0xe0db3ef21d03bad8,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk4/icao-sea-level/rh100/no-zones",
        &[
            0x40b26cbe6c065772, 0x4091b511753928bd, 0x403def248e75797a, 0x40619dadb45f58c6,
            0x405c26fed3c02650, 0x40755334ef1fc968, 0x0000000000000001, 0x40851f13695b2735,
            0x0000000000000001, 0x4089a6ebe85f3d22, 0x0000000000000001, 0x408e3e41f3e0c38d,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x00000000000074f0,
            0x3f50624dd2f1a9fc, 0x3fe6cce94d8f0784, 0x3fda53cde9b40d50, 0xbeb537056ef7b260,
            0x4089b3cc7958d2f3, 0x0000000000000001, 0x3fd313e8c762a381, 0x401def9db22d1201,
            0x40a079e77572fdd8, 0x408fbd84c5f8d90a, 0xc03219a92df2ceda, 0x4067d42c6a230eba,
            0x0000000000000001, 0x3fc0933a6617bb12, 0x402def9db22d080e, 0x40a96d1b8bfd1152,
            0x4091236fada5e65d, 0xc044cdb0613fdaf7, 0x40607f36ebe112ea, 0x0000000000000001,
            0x3fc093e578fd4e77, 0x403673b645a1cfa7, 0x40affeaafe9dd3f6, 0x40853c2c5a81a962,
            0xc050b4c247ff595a, 0x405fe23cbf4ed9f8, 0x0000000000000001, 0x3fc09501bdcd2a6b,
            0x403def248e75797a, 0x40b26cbe6c065772, 0xc059000000000000, 0xc05750fe3a4058a8,
            0x40619dadb45f58c6, 0x0000000000000001, 0x3fc092c70f2d7350, 0x1f133dad9b3ad421,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk4/icao-sea-level/rh50/zones-rh20-rh80",
        &[
            0x40b23fba1dc7a82b, 0x409194bfb578fdec, 0x403dd20a12014ea5, 0x406182d87cd77cad,
            0x405bd17c84438eb3, 0x40754bf73b88949e, 0x0000000000000001, 0x40851e4bc3dfd78f,
            0x0000000000000001, 0x40899772cae26c8c, 0x0000000000000001, 0x408def76acdbb7ac,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x000000000000747e,
            0x3f50624dd2f1a9fc, 0x3fe6cceb57ff745a, 0x3fda53d044f595bd, 0xbeb51f90ae03f360,
            0x4089b3d113360e07, 0x0000000000000001, 0x3fd310d8973a29ff, 0x401dd1eb851ebbf2,
            0x40a05d210f4febcc, 0x408f8b2bcdbeccd3, 0xc03210a2b295507d, 0x4067ab4e885ba9eb,
            0x0000000000000001, 0x3fc0932aa135a05d, 0x402dd26e978d49a9, 0x40a93764e0852075,
            0x40910258f3b6be9b, 0xc044c3cc361aeb8c, 0x4060606d3a2299a5, 0x0000000000000001,
            0x3fc094144fb40aea, 0x40365db22d0e5acd, 0x40afb4d63098bbde, 0x40850cd6d084b967,
            0xc050ac008f60fdc3, 0x405faed7f34a07a3, 0x0000000000000001, 0x3fc0952be2d5f49f,
            0x403dd20a12014ea5, 0x40b23fba1dc7a82b, 0xc059000000000000, 0xc057440f0429c752,
            0x406182d87cd77cad, 0x0000000000000001, 0x3fc092d00f10dd01, 0xf3439173982e2e2d,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk4/1500m-25C-830hPa/rh50/no-zones",
        &[
            0x40b5a9d162b9a869, 0x409455f3e8bad6a2, 0x403ff87fa4f7fc58, 0x40637cdbf42c6722,
            0x406139f8291e1d5c, 0x4075b318f42cf8e0, 0x0000000000000001, 0x4089f3e80a438605,
            0x0000000000000001, 0x408fb46ec3698834, 0x0000000000000001, 0x4092b6dd6164cb7e,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000007ce4,
            0x3f50624dd2f1a9fc, 0x3fe6cd490d8df337, 0x3fda543c797e4377, 0xbeb0ea802787ff89,
            0x4089b4a45a2f0a80, 0x0000000000000001, 0x3fd33b6a5355dfb8, 0x401ff8d4fdf3bab0,
            0x40a2cdb4df4ce691, 0x40921600f7be460e, 0xc03276e0fcc9502a, 0x406a0d532dfec014,
            0x0000000000000001, 0x3fc096b7db3fe5b2, 0x402ff8d4fdf3aec0, 0x40ad6448ab2a4199,
            0x4093c3e3c16bd24e, 0xc0456a46b4f90e4d, 0x40626e728566afcb, 0x0000000000000001,
            0x3fc092ad00bd32bd, 0x4037fa9fbe76cfa8, 0x40b2aab0e9fe3c9e, 0x4088e3f6a5dcfc7c,
            0xc051463965885e3a, 0x4061b87319b8c220, 0x0000000000000001, 0x3fc092dfafbea5fb,
            0x403ff87fa4f7fc58, 0x40b5a9d162b9a869, 0xc059000000000000, 0xc05830905a403b43,
            0x40637cdbf42c6722, 0x0000000000000001, 0x3fc093162d605051, 0xdc5f005783aa1e64,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk4/1500m-25C-830hPa/rh100/no-zones",
        &[
            0x40b5c82bda519754, 0x40946ef6f61bbe97, 0x4040061d114b1f9c, 0x40638c65bf4ab879,
            0x4061557b8f7f0898, 0x4075c3b860e9c586, 0x0000000000000001, 0x4089fec28cf49b17,
            0x0000000000000001, 0x408fc95bc558a87c, 0x0000000000000001, 0x4092c70e545c65e3,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000007d31,
            0x3f50624dd2f1a9fc, 0x3fe6cd4b4257d3c1, 0x3fda543f05a653d1, 0xbeb0d124a64ac850,
            0x4089b4a9539344f5, 0x0000000000000001, 0x3fd34233bb7509b0, 0x40200624dd2f1cd2,
            0x40a2e3ff15572440, 0x40922abd23883f9a, 0xc0327a52a27780d9, 0x406a237f2f978100,
            0x0000000000000001, 0x3fc096bf26e7b4e9, 0x40300624dd2f16e3, 0x40ad8a620cfe8d45,
            0x4093dcea9e8f4d52, 0xc045700c3e967de8, 0x40627ede60b84eb4, 0x0000000000000001,
            0x3fc092acfcdbf94e, 0x403809374bc6aef7, 0x40b2c3f50078ffde, 0x4089071fcf8e9994,
            0xc0514bbdc14767b2, 0x4061c6e40c5e93da, 0x0000000000000001, 0x3fc092df370ad25b,
            0x4040061d114b1f9c, 0x40b5c82bda519754, 0xc059000000000000, 0xc0583963dd25dc28,
            0x40638c65bf4ab879, 0x0000000000000001, 0x3fc09316a828b5b4, 0x2632a1281854f4a9,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk4/1500m-25C-830hPa/rh50/zones-rh20-rh80",
        &[
            0x40b586795bd7cb71, 0x40943e41340ba4fb, 0x403fe420875b9edf, 0x40636878415c5688,
            0x406115fef8eb50bc, 0x4075b318f42cf8e0, 0x0000000000000001, 0x4089eb12ea8c1305,
            0x0000000000000001, 0x408f9e5f4a6e948e, 0x0000000000000001, 0x40928d5ff308aefb,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000007c95,
            0x3f50624dd2f1a9fc, 0x3fe6cd4b25bfcc0b, 0x3fda543ee4a1f3ac, 0xbeb0d26d4c72cc2a,
            0x4089b4a9131ee34a, 0x0000000000000001, 0x3fd33b68b217bb5d, 0x401fe45a1cac0c94,
            0x40a2b8ec9bdfdade, 0x4092042e1c2b6576, 0xc0327153f2d42e6a, 0x4069ef440c3a97e3,
            0x0000000000000001, 0x3fc096576a984bf6, 0x402fe45a1cac00b8, 0x40ad3b9212509fdb,
            0x4093ab84e2cbd760, 0xc045646b647c1200, 0x406256d22c347fa2, 0x0000000000000001,
            0x3fc092ae0ff9240e, 0x4037eb4395810d04, 0x40b28e261ec75f06, 0x4088c058652b363b,
            0xc051411dacc06c1b, 0x4061a48cb59d613c, 0x0000000000000001, 0x3fc092eb27d69dc8,
            0x403fe420875b9edf, 0x40b586795bd7cb71, 0xc059000000000000, 0xc05828e3165d2137,
            0x40636878415c5688, 0x0000000000000001, 0x3fc09306b9ee45c4, 0xca0e948203a0ec77,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk45/icao-sea-level/rh50/no-zones",
        &[
            0x40b2612e4a3d8e84, 0x4091ab336332540a, 0x403de6dd0af37fa1, 0x4061977687f384c5,
            0x405c1324916562c8, 0x40754bf73b88949e, 0x0000000000000001, 0x408520775d6296e1,
            0x0000000000000001, 0x4089a0d254700781, 0x0000000000000001, 0x408e3aae053b0ce2,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000000bb2,
            0x3f50624dd2f1a9fc, 0x3fe6cce81838996a, 0x3fda53cc84838bcf, 0xbeb544e8dd375f42,
            0x4089b3c9bfecfcf4, 0x0000000000000001, 0x3fd310db27e5b0a1, 0x401dceb78897e915,
            0x40a0683760770127, 0x408fa2d8f9496bd5, 0xc0320577ca49b7a2, 0x4067d59fc5dc2ed0,
            0x0000000000000001, 0x3fc0934039f9cb88, 0x402de23d0bfa08ac, 0x40a95ba5c68e6fc2,
            0x40911aa0d8695bb4, 0xc044c6952c5d9d0c, 0x406079bb9aed7737, 0x0000000000000001,
            0x3fc093e54179d52a, 0x40366bffcdab19e5, 0x40afea08d1518b72, 0x40853290afee7a58,
            0xc050b0c06ed98053, 0x405fd641b6947aef, 0x0000000000000001, 0x3fc09502f1c99c87,
            0x403de6dd0af37fa1, 0x40b2612e4a3d8e84, 0xc059000000000000, 0xc0574cdded74b55f,
            0x4061977687f384c5, 0x0000000000000001, 0x3fc092c727228bba, 0xc8052fff95d015e1,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk45/icao-sea-level/rh100/no-zones",
        &[
            0x40b26cbe69e80d82, 0x4091b5117589f3ac, 0x403def248ac26a05, 0x40619dadb1fec992,
            0x405c26fecc270d12, 0x40755334ef1fc968, 0x0000000000000001, 0x408524720c089bab,
            0x0000000000000001, 0x4089be9e5f468d90, 0x0000000000000001, 0x408e4256977ebfd9,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000000bb6,
            0x3f50624dd2f1a9fc, 0x3fe6cce94d8f06f6, 0x3fda53cde9b40cac, 0xbeb537056efe0fc7,
            0x4089b3cc7958d2f3, 0x0000000000000001, 0x3fd313e8c762a381, 0x401dd8f4f93bc01f,
            0x40a071c978b4c3e5, 0x408fb43cf2c82cf2, 0xc03208d530de73f9, 0x4067dde22feec3f2,
            0x0000000000000001, 0x3fc093405372b59b, 0x402dec7a7c9ddfb6, 0x40a96b858358a1d1,
            0x409124111af59f1f, 0xc044cb1f024c955b, 0x40607fdd6fea28cb, 0x0000000000000001,
            0x3fc093e44c748e8e, 0x403673ade225fb2e, 0x40affea4a7ccbf9d, 0x40853c421cd91b7e,
            0xc050b4bb07f424c5, 0x405fe23a08b8a863, 0x0000000000000001, 0x3fc09501c0a02cbb,
            0x403def248ac26a05, 0x40b26cbe69e80d82, 0xc059000000000000, 0xc05750fe37501d2a,
            0x40619dadb1fec992, 0x0000000000000001, 0x3fc092c70f2e6637, 0x3f06fd9fc4090c84,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk45/icao-sea-level/rh50/zones-rh20-rh80",
        &[
            0x40b23fb5174dbcd2, 0x409194b872f31184, 0x403dd204c043943c, 0x406182d74e2d81f5,
            0x405bd178c2a2a218, 0x40754bf73b88949e, 0x0000000000000001, 0x4085267d9e5eb1eb,
            0x0000000000000001, 0x4089a57cb48bcc17, 0x0000000000000001, 0x408dfa035970b838,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000000baa,
            0x3f50624dd2f1a9fc, 0x3fe6cceb57ff73d1, 0x3fda53d044f5951f, 0xbeb51f90ae0a12e4,
            0x4089b3d113360e07, 0x0000000000000001, 0x3fd310d8973a29ff, 0x401dba3ca7503b01,
            0x40a054aeb33ce3b8, 0x408f817f7a015967, 0xc031ff052ca408fe, 0x4067b571afd0915c,
            0x0000000000000001, 0x3fc0932ecdd5de7e, 0x402dcdc22ab25a98, 0x40a93506ccff6119,
            0x409103420798a800, 0xc044bff979abb048, 0x406061625a11048f, 0x0000000000000001,
            0x3fc094127c9e1729, 0x40365ca3a4b55753, 0x40afb404c1a2b047, 0x40850f7f76c080de,
            0xc050ab17a3aef268, 0x405fae7f3ae760b0, 0x0000000000000001, 0x3fc0952c3e983c40,
            0x403dd204c043943c, 0x40b23fb5174dbcd2, 0xc059000000000000, 0xc057440b807ef90f,
            0x406182d74e2d81f5, 0x0000000000000001, 0x3fc092d00f9d29c8, 0xb63cf3e571120e65,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk45/1500m-25C-830hPa/rh50/no-zones",
        &[
            0x40b5a9d162038269, 0x409455f3ea204fb5, 0x403ff87fa360bb89, 0x40637cdbf31a085f,
            0x406139f827390c5b, 0x4075b318f42cf8e0, 0x0000000000000001, 0x4089f67528f74c12,
            0x0000000000000001, 0x408fb692af18dbf6, 0x0000000000000001, 0x4092c0bfe6bb7f63,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000000c81,
            0x3f50624dd2f1a9fc, 0x3fe6cd490d8df305, 0x3fda543c797e433b, 0xbeb0ea80278a454a,
            0x4089b4a45a2f0a80, 0x0000000000000001, 0x3fd33b6a5355dfb8, 0x401fe33269df971d,
            0x40a2c53d458d092a, 0x4092111f933718b1, 0xc0326753845f2b5b, 0x406a1623c5406547,
            0x0000000000000001, 0x3fc096d5260159e1, 0x402ff19934efcb2f, 0x40ad602eac522915,
            0x4093c559b9d7d615, 0xc04564840e89adb0, 0x40626ff9e2b43087, 0x0000000000000001,
            0x3fc092acfdb65c89, 0x4037f8cc9a77e6f6, 0x40b2a9e6cff004e4, 0x4088e901604feae6,
            0xc05144aed8b2d346, 0x4061b82d7bb7acae, 0x0000000000000001, 0x3fc092dfd5da1711,
            0x403ff87fa360bb89, 0x40b5a9d162038269, 0xc059000000000000, 0xc058309058ea7667,
            0x40637cdbf31a085f, 0x0000000000000001, 0x3fc093162d5f7825, 0xbbfcc2aef420d53d,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk45/1500m-25C-830hPa/rh100/no-zones",
        &[
            0x40b5c82bd5453e0e, 0x40946ef6f1db535e, 0x4040061d0de84d0e, 0x40638c65bc1b452c,
            0x4061557b89d8e2fc, 0x4075c3b860e9c586, 0x0000000000000001, 0x408a0151bf3d8ed0,
            0x0000000000000001, 0x408fdf2dfab08917, 0x0000000000000001, 0x4092cb624c7129aa,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000000c89,
            0x3f50624dd2f1a9fc, 0x3fe6cd4b4257d392, 0x3fda543f05a6539b, 0xbeb0d124a64ce5dd,
            0x4089b4a9539344f5, 0x0000000000000001, 0x3fd34233bb7509b0, 0x401ff7ad4b274531,
            0x40a2dbe754482bb1, 0x4092261295567a6d, 0xc0326b8674d8a75f, 0x406a2be4ca354092,
            0x0000000000000001, 0x3fc096db11f4b413, 0x4030030a0b1bbca2, 0x40ad86d97c30d5fe,
            0x4093de2bc5b5b747, 0xc0456b1b80cff650, 0x4062802eed1428e1, 0x0000000000000001,
            0x3fc092acfc6dbf18, 0x40380828c36da988, 0x40b2c37f819855f9, 0x40890a0ce990caed,
            0xc0514ad96c98bf02, 0x4061c6bbd8038d81, 0x0000000000000001, 0x3fc092df4cdff4b8,
            0x4040061d0de84d0e, 0x40b5c82bd5453e0e, 0xc059000000000000, 0xc0583963d830750f,
            0x40638c65bc1b452c, 0x0000000000000001, 0x3fc09316a82633b0, 0xcbf93d3fe086fce8,
            0xffffffffffffffff,
        ],
    ),
    (
        "lob30/rk45/1500m-25C-830hPa/rh50/zones-rh20-rh80",
        &[
            0x40b586761083b380, 0x40943e3c83bd5748, 0x403fe41d62478142, 0x40636877833d5b1d,
            0x406115fdaa2c4554, 0x4075b318f42cf8e0, 0x0000000000000001, 0x4089f7acee1b6ed8,
            0x0000000000000001, 0x408f9caa3b5efd52, 0x0000000000000001, 0x4092910d6eefc36e,
            0x0000000000000000, 0x0000000000000000, 0xe585c7f5c416fe29, 0x0000000000000c79,
            0x3f50624dd2f1a9fc, 0x3fe6cd4b25bfcbdd, 0x3fda543ee4a1f377, 0xbeb0d26d4c74f5de,
            0x4089b4a9131ee34a, 0x0000000000000001, 0x3fd33b68b217bb5d, 0x401fce261559e847,
            0x40a2b042e7b42656, 0x4091ff2e9088a554, 0xc0326159051a92fb, 0x4069f84eeb2d1628,
            0x0000000000000001, 0x3fc09673ddf54c6c, 0x402fdcd59a091cba, 0x40ad3750cc91c040,
            0x4093ad05502a8848, 0xc0455e6e1fed0cc3, 0x406258667aef7e76, 0x0000000000000001,
            0x3fc092adf0b12b6b, 0x4037e94c14b2a434, 0x40b28d4b7a8b9c46, 0x4088c5b9157d51f2,
            0xc0513f74b3ee006c, 0x4061a440eeef4d7c, 0x0000000000000001, 0x3fc092eb55c4cbd4,
            0x403fe41d62478142, 0x40b586761083b380, 0xc059000000000000, 0xc05828e12030b8d5,
            0x40636877833d5b1d, 0x0000000000000001, 0x3fc09306b963eb85, 0x9106aeae156784c9,
            0xffffffffffffffff,
        ],
    ),
];

/// `TrajectorySolver::solve` (RK4, flat) with Litz spin drift plus aerodynamic jump, and with
/// Magnus, at both stations, RH 0.
///
/// IMMUTABLE: nothing in 0.50.0 may change it, the aloft vapour model included.
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
/// MAY CHANGE ONLY WITH THE ALOFT VAPOUR MODEL: byte-identical through every other change in
/// 0.50.0.
///
/// Re-pinned by the aloft vapour model (0.50.0): the unified aloft vapour model holds the
/// governing air's (the station's or the zone's) water-vapour mole fraction for the lapse, the
/// density and the sound speed, with a virtual-temperature lapse. 260 of 432 values moved, all in
/// cases with a humidity above 0 and all off the station's own height; every RH-0 and
/// station-height value is unchanged. The largest float shift of each moved case,
/// pinned before the aloft vapour model -> pinned now:
/// - `spin-drift+aerodynamic-jump/rk4/icao-sea-level/rh50` points[1963].drag_coefficient: 0.22457961110519453 -> 0.22457579017393242 (-1.70e-05)
/// - `magnus/rk4/icao-sea-level/rh50` points[1963].drag_coefficient: 0.22457761861610434 -> 0.22457378373250123 (-1.71e-05)
/// - `spin-drift+aerodynamic-jump/rk4/icao-sea-level/rh100` points[1960].drag_coefficient: 0.2270644023671422 -> 0.2270565787292812 (-3.45e-05)
/// - `magnus/rk4/icao-sea-level/rh100` points[1960].drag_coefficient: 0.2270622832017817 -> 0.22705443303405887 (-3.46e-05)
/// - `spin-drift+aerodynamic-jump/rk4/1500m-25C-830hPa/rh50` impact_energy: 903.6885164129332 -> 903.6932692212704 (+5.26e-06)
/// - `magnus/rk4/1500m-25C-830hPa/rh50` impact_energy: 903.6845756762455 -> 903.6894691302815 (+5.42e-06)
/// - `spin-drift+aerodynamic-jump/rk4/1500m-25C-830hPa/rh100` impact_energy: 912.8529856128719 -> 912.862453914296 (+1.04e-05)
/// - `magnus/rk4/1500m-25C-830hPa/rh100` impact_energy: 912.8490154369877 -> 912.8587649036949 (+1.07e-05)
#[rustfmt::skip]
const TS_EFFECTS_HUMID: Pinned = &[
    (
        "spin-drift+aerodynamic-jump/rk4/icao-sea-level/rh50",
        &[
            0x408f400000000000, 0x3fc296032ae97b20, 0x3fff64dd321c9bd6, 0x40743e61cc300a18,
            0x408296bd4dd157b6, 0x40754bf73b88949e, 0x0000000000000001, 0x40882eb171145ab3,
            0x0000000000000001, 0x408d487d23cae347, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ac,
            0x3f50624dd2f1a9fc, 0x3fea53dd2d38b581, 0x3f5c8f673e64da92, 0xbeada2fd280448e2,
            0x4089b3d3b8df5ebe, 0x0000000000000001, 0x3fd310d7a9ee2551, 0x3fdf6c8b4395810d,
            0x4075de59eefc96f9, 0xbfd55519361cb35a, 0xbfcd93d5ece2dac8, 0x40834d823296023f,
            0x0000000000000001, 0x3fd594c71a59e7db, 0x3fef6c8b4395810d, 0x40834d8a956e9aa2,
            0xc00568590131866e, 0xbfe9f2e2653dc499, 0x407e15e57ff4a5db, 0x0000000000000001,
            0x3fd8bdfb2dc4428a, 0x3ff7916872b01fdd, 0x4089e792a78c4250, 0xc01a8cb989c0c9ca,
            0xbff9f274f793a850, 0x40780b5bd3a17387, 0x0000000000000001, 0x3fdba9b1ac827c6c,
            0x3fff64dd321c9bd6, 0x408f400000000000, 0xc027f80176d9ae25, 0xc004a46544936589,
            0x40743e61cc300a18, 0x0000000000000001, 0x3fccbee645229c66, 0x1b0183cc1489af44,
            0x0000000000000015, 0xc3ffda9c65778e2a,
        ],
    ),
    (
        "magnus/rk4/icao-sea-level/rh50",
        &[
            0x408f400000000000, 0x3fc0887808691bd0, 0x3fff64e2add69864, 0x40743e60dde64c96,
            0x408296bb98334813, 0x40754bf73b88949e, 0x0000000000000001, 0x40882eb005a022e8,
            0x0000000000000001, 0x408d487a4bc2ec3f, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007ac,
            0x3f50624dd2f1a9fc, 0x3fea53dd9ed23649, 0x3f5ae111cdf47712, 0xbeb544fd87f27ec9,
            0x4089b3d3b9831169, 0x0000000000000001, 0x3fd310d7a9b4d2da, 0x3fdf6c8b4395810d,
            0x4075de59fc716005, 0xbfd8206f55f5e529, 0xbfd0db9b97eaac3f, 0x40834d826cdcf7bd,
            0x0000000000000001, 0x3fd594c6fc123e87, 0x3fef6c8b4395810d, 0x40834d8a01bf1dd0,
            0xc00606560d8449da, 0xbfeda08358a420a9, 0x407e15e4ff193d85, 0x0000000000000001,
            0x3fd8bdfb6a0436c9, 0x3ff7916872b01fdd, 0x4089e790d2229853, 0xc01af6dc5c089428,
            0xbffdcf3bd96f9bc2, 0x40780b5a81877a7d, 0x0000000000000001, 0x3fdba9b2325dba28,
            0x3fff64e2add69864, 0x408f400000000000, 0xc028382e08c36af1, 0xc007e7ec6e8de883,
            0x40743e60dde64c96, 0x0000000000000001, 0x3fccbed57055c4d7, 0xb85beca88f1ee07b,
            0x0000000000000015, 0xab44e538aa017ab8,
        ],
    ),
    (
        "spin-drift+aerodynamic-jump/rk4/icao-sea-level/rh100",
        &[
            0x408f400000000000, 0x3fc296d48b074389, 0x3fff58c12e000d45, 0x407449a63935bae6,
            0x4082ab7454699b8b, 0x40755334ef1fc968, 0x0000000000000001, 0x408836a2ae27571d,
            0x0000000000000001, 0x408d52091f196239, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007a9,
            0x3f50624dd2f1a9fc, 0x3fea53de926cece6, 0x3f5c8f68c24087f2, 0xbead87360c00dfb0,
            0x4089b3d672551666, 0x0000000000000001, 0x3fd313e54a033ba9, 0x3fdf5c28f5c28f63,
            0x4075d664e96549ab, 0xbfd52577e56ee2ca, 0xbfcd63e9a0b35950, 0x408353931232fe65,
            0x0000000000000001, 0x3fd5950927c310fe, 0x3fef5c28f5c28f63, 0x408348d0599d8440,
            0xc0054ef79e59401e, 0xbfe9cc6215f4d39c, 0x407e2684a14f1397, 0x0000000000000001,
            0x3fd8baff5bf13d68, 0x3ff7851eb851ea9f, 0x4089e3ac6a36410b, 0xc01a71f2e0ecd427,
            0xbff9ceb0be914c3e, 0x40781d29333441d3, 0x0000000000000001, 0x3fdba5dd6ac71b9e,
            0x3fff58c12e000d45, 0x408f400000000000, 0xc027e7b31a66182a, 0xc0048e7a9d32b560,
            0x407449a63935bae6, 0x0000000000000001, 0x3fcd1030a1fdf0d4, 0x622e7db8161535eb,
            0x0000000000000015, 0x9b45f94b2f16ac37,
        ],
    ),
    (
        "magnus/rk4/icao-sea-level/rh100",
        &[
            0x408f400000000000, 0x3fc089292bf2c127, 0x3fff58c69e13a359, 0x407449a546df805b,
            0x4082ab729664522a, 0x40755334ef1fc968, 0x0000000000000001, 0x408836a143e152f6,
            0x0000000000000001, 0x408d5206492c57f7, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000007a9,
            0x3f50624dd2f1a9fc, 0x3fea53df040673e6, 0x3f5ae1133b04c6ef, 0xbeb53719f9eec8b9,
            0x4089b3d672f8c9b1, 0x0000000000000001, 0x3fd313e549c9f2f7, 0x3fdf5c28f5c28f63,
            0x4075d664f78f4a97, 0xbfd7efc95e051922, 0xbfd0c1ac7745667a, 0x408353934dc2f77d,
            0x0000000000000001, 0x3fd5950908d7e4c5, 0x3fef5c28f5c28f63, 0x408348cfc7f63fde,
            0xc005eccdbae07faf, 0xbfed768154d0b32d, 0x407e26842638ee3c, 0x0000000000000001,
            0x3fd8baff955ebeb3, 0x3ff7851eb851ea9f, 0x4089e3aa99740ad4, 0xc01adc055c829418,
            0xbffda7c8f92757c2, 0x40781d27e7c0ecc7, 0x0000000000000001, 0x3fdba5ddef645d28,
            0x3fff58c69e13a359, 0x408f400000000000, 0xc02827df243aaf8b, 0xc007cfb44f07fe63,
            0x407449a546df805b, 0x0000000000000001, 0x3fcd101ea2258476, 0x2c8ca9da48b1befc,
            0x0000000000000015, 0x774a224ff08633c5,
        ],
    ),
    (
        "spin-drift+aerodynamic-jump/rk4/1500m-25C-830hPa/rh50",
        &[
            0x408f400000000000, 0x3fc325d64f85c97b, 0x3ffbe89bc57194ce, 0x4078f39e5cc0738e,
            0x408c3d8bd0bbc570, 0x4075b318f42cf8e0, 0x0000000000000001, 0x408db7b8c02460aa,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d2,
            0x3f50624dd2f1a9fc, 0x3fea544d126f7afa, 0x3f5ccbfa418f132c, 0xbea29fa3ed806134,
            0x4089b4ae549f602d, 0x0000000000000001, 0x3fd33b66e1b8145d, 0x3fdbe76c8b439587,
            0x407436663eb27692, 0xbfc6ccb9a71805b4, 0xbfc21d7173accb21, 0x4084e6e9922f01d5,
            0x0000000000000001, 0x3fd502ea026ef1b8, 0x3febef9db22d0e5c, 0x40826d9f96b1e7a9,
            0xc0003bb82ef86c75, 0xbfe0884397e65e19, 0x408159f57a360820, 0x0000000000000001,
            0x3fd708c637204d18, 0x3ff4f1a9fbe76bf5, 0x40895f81d387ec2c, 0xc015261797f496dd,
            0xbff0f7bcd6650e5c, 0x407d3d987a895c05, 0x0000000000000001, 0x3fd96e159b954ac0,
            0x3ffbe89bc57194ce, 0x408f400000000000, 0xc023860c78f9fd30, 0xbffbaaa91b86832f,
            0x4078f39e5cc0738e, 0x0000000000000001, 0x3fdb79a9c9e8f4d7, 0x3290b02ebdde744b,
            0x0000000000000015, 0xe95887fa16b0ed66,
        ],
    ),
    (
        "magnus/rk4/1500m-25C-830hPa/rh50",
        &[
            0x408f400000000000, 0x3fc0c08f80dd5718, 0x3ffbe89f7283c58f, 0x4078f39aec972251,
            0x408c3d8408643149, 0x4075b318f42cf8e0, 0x0000000000000001, 0x408db7b65d3ecf23,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006d2,
            0x3f50624dd2f1a9fc, 0x3fea544d94744d24, 0x3f5ae1843965d216, 0xbeb0ea90025d2510,
            0x4089b4ae555a0ca3, 0x0000000000000001, 0x3fd33b66e177a727, 0x3fdbe76c8b439587,
            0x407436666b2ebdc0, 0xbfccaf788166b828, 0xbfc60a3415676703, 0x4084e6ea05a99a2c,
            0x0000000000000001, 0x3fd502e9d39cb202, 0x3febef9db22d0e5c, 0x40826d9f45f56560,
            0xc000e788df4941a5, 0xbfe40737f1c97aa9, 0x408159f59ecb2999, 0x0000000000000001,
            0x3fd708c61be48364, 0x3ff4f1a9fbe76bf5, 0x40895f808341308e, 0xc0159c758957766a,
            0xbff4a2d8dcedcdf3, 0x407d3d97ce3d89ef, 0x0000000000000001, 0x3fd96e15ef0efd17,
            0x3ffbe89f7283c58f, 0x408f400000000000, 0xc023cf06938239fb, 0xc000ef65f7e84223,
            0x4078f39aec972251, 0x0000000000000001, 0x3fdb79ab4c9a1767, 0xf8d105441fc6a817,
            0x0000000000000015, 0x14523b9061886c8b,
        ],
    ),
    (
        "spin-drift+aerodynamic-jump/rk4/1500m-25C-830hPa/rh100",
        &[
            0x408f400000000000, 0x3fc32769fa7eb1fc, 0x3ffbd682b0284951, 0x407913f13b632260,
            0x408c86e64e3ce1ab, 0x4075c3b860e9c586, 0x0000000000000001, 0x408dc7e20a360ede,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006cd,
            0x3f50624dd2f1a9fc, 0x3fea544f9e9e172d, 0x3f5ccbfd0b86404f, 0xbea26cec7b7a754c,
            0x4089b4b34e14fb64, 0x0000000000000001, 0x3fd3423047934bfd, 0x3fdbd70a3d70a3dd,
            0x40742eb147d612cc, 0xbfc67a81dec61e97, 0xbfc1eb4c2f47f3c4, 0x4084ef8853827d06,
            0x0000000000000001, 0x3fd505e8d1972a7a, 0x3febd70a3d70a3dd, 0x4082653d2dba72bc,
            0xc00018ee60ab763d, 0xbfe0554efff181ed, 0x408168a73e284bea, 0x0000000000000001,
            0x3fd707bbed5ac193, 0x3ff4e147ae147a4d, 0x408959da9a6155f7, 0xc0150516457774cb,
            0xbff0c983a63e61ff, 0x407d5e51367d99c3, 0x0000000000000001, 0x3fd969172b4cf3d6,
            0x3ffbd682b0284951, 0x408f400000000000, 0xc023701e7c25f380, 0xbffb699f5e068142,
            0x407913f13b632260, 0x0000000000000001, 0x3fdb73dca9bf498e, 0x25da1d1ae01e2c4c,
            0x0000000000000015, 0xf3bee55b885b6632,
        ],
    ),
    (
        "magnus/rk4/1500m-25C-830hPa/rh100",
        &[
            0x408f400000000000, 0x3fc0c1db1398e2e0, 0x3ffbd68648282306, 0x407913ede9417ee5,
            0x408c86dec0224294, 0x4075c3b860e9c586, 0x0000000000000001, 0x408dc7dfad634539,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0xd34b230aaf8a1446, 0x00000000000006cd,
            0x3f50624dd2f1a9fc, 0x3fea545020a2f653, 0x3f5ae186d3ef5e72, 0xbeb0d13449562df5,
            0x4089b4b34ecfa91d, 0x0000000000000001, 0x3fd342304752b46f, 0x3fdbd70a3d70a3dd,
            0x40742eb1752fd283, 0xbfcc5b0177f0bd52, 0xbfc5d3d834174663, 0x4084ef88c931883e,
            0x0000000000000001, 0x3fd505e8a1c7aa4f, 0x3febd70a3d70a3dd, 0x4082653ce025d652,
            0xc000c470988a653c, 0xbfe3cea493f92eb4, 0x408168a768d4fa4b, 0x0000000000000001,
            0x3fd707bbcdbabfaf, 0x3ff4e147ae147a4d, 0x408959d951bd1d71, 0xc0157b59766573bd,
            0xbff46f60eebf9b73, 0x407d5e509bef624f, 0x0000000000000001, 0x3fd9691776323f9d,
            0x3ffbd68648282306, 0x408f400000000000, 0xc023b9181713eb95, 0xc000cb33a330fcf8,
            0x407913ede9417ee5, 0x0000000000000001, 0x3fdb73de1fe539a4, 0x4169f6fdb3619709,
            0x0000000000000015, 0x4a98230d824aa06c,
        ],
    ),
];

/// `calculate_zero_angle_with_station_presence` to 300 m (every zero trial flies
/// `calculate_acceleration`), Euler, RK4 and RK45 at both stations, RH 0. The bisection stops at
/// 1e-7 rad, so a zero pin moves only when a trial height changes sign; the flight groups are the
/// sensitive ones.
///
/// IMMUTABLE: nothing in 0.50.0 may change it, the aloft vapour model included.
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
/// MAY CHANGE ONLY WITH THE ALOFT VAPOUR MODEL: byte-identical through every other change in
/// 0.50.0.
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
/// MAY CHANGE ONLY WITH THE ALOFT VAPOUR MODEL: byte-identical through every other change in
/// 0.50.0.
///
/// Re-pinned by the aloft vapour model (0.50.0): the unified aloft vapour model holds the
/// governing air's (the station's or the zone's) water-vapour mole fraction for the lapse, the
/// density and the sound speed, with a virtual-temperature lapse. 14 of 18 values moved, all in
/// cases with a humidity above 0 and all off the station's own height; every RH-0 and
/// station-height value is unchanged. The largest float shift of each moved case,
/// pinned before the aloft vapour model -> pinned now:
/// - `seeded/icao-sea-level/default-rh` impact_positions[0].z: 0.11563447090134216 -> 0.11563411086056596 (-3.11e-06)
/// - `seeded/1500m-25C-830hPa/default-rh` impact_velocities[0]: 563.6048728741576 -> 563.6045855884292 (-5.10e-07)
///
/// glibc rounds a few of these a last bit apart; its bits are in `LIBM_ALTERNATES`.
#[rustfmt::skip]
const MONTE_CARLO_DEFAULT_HUMIDITY: Pinned = &[
    (
        "seeded/icao-sea-level/default-rh",
        &[
            0x0000000000000028, 0x408283ebd09f0d48, 0x407de78b9d14e1e8, 0x3fe9a387cfcc2317,
            0x3fbd9a327472bb38, 0x3fb999999999999a, 0xe7a1e54bc5b05916, 0x3eb70a0652aa8412,
            0xdb57b2d42f00978a,
        ],
    ),
    (
        "seeded/1500m-25C-830hPa/default-rh",
        &[
            0x0000000000000028, 0x40811edd3e5d37c7, 0x40819cd630f80f7e, 0x3fc6881d9f864343,
            0xbfe4650065cacf50, 0x3fb999999999999a, 0xbeeb4722e7193213, 0x471e7f90deacdc40,
            0x287c34599195a896,
        ],
    ),
];

/// Station density and sound speed: `calculate_atmosphere` over altitude, T/P override and RH, plus
/// `calculate_air_density_cimp` and `moist_speed_of_sound`. Station air is not part of the aloft
/// model, so this holds through the aloft vapour model at every RH.
///
/// IMMUTABLE: nothing in 0.50.0 may change it, the aloft vapour model included.
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
/// IMMUTABLE: nothing in 0.50.0 may change it, the aloft vapour model included.
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
/// MAY CHANGE ONLY WITH THE ALOFT VAPOUR MODEL: byte-identical through every other change in
/// 0.50.0.
///
/// Re-pinned by the aloft vapour model (0.50.0): the unified aloft vapour model holds the
/// governing air's (the station's or the zone's) water-vapour mole fraction for the lapse, the
/// density and the sound speed, with a virtual-temperature lapse. 108 of 135 values moved, all in
/// cases with a humidity above 0 and all off the station's own height; every RH-0 and
/// station-height value is unchanged. The largest float shift of each moved case,
/// pinned before the aloft vapour model -> pinned now:
/// - `get_local_atmosphere_humid/0m/15C/1013.25hPa/rh37.5` +25000m.density: 0.04000548848073421 -> 0.04036000743432573 (+8.86e-03)
/// - `get_local_atmosphere_humid/0m/15C/1013.25hPa/rh50` +25000m.density: 0.0399736963957377 -> 0.04044670891224762 (+1.18e-02)
/// - `get_local_atmosphere_humid/0m/15C/1013.25hPa/rh100` +25000m.density: 0.03984673206914066 -> 0.04079533176434715 (+2.38e-02)
/// - `get_local_atmosphere_humid/0m/35C/1013.25hPa/rh37.5` +25000m.density: 0.049017317155805644 -> 0.050348589699717734 (+2.72e-02)
/// - `get_local_atmosphere_humid/0m/35C/1013.25hPa/rh50` +25000m.density: 0.048888234150328815 -> 0.05066655557263168 (+3.64e-02)
/// - `get_local_atmosphere_humid/0m/35C/1013.25hPa/rh100` +25000m.density: 0.048373670783280254 -> 0.051956884918299584 (+7.41e-02)
/// - `get_local_atmosphere_humid/1500m/25C/830hPa/rh37.5` +25000m.density: 0.038310529792086746 -> 0.03903155870695048 (+1.88e-02)
/// - `get_local_atmosphere_humid/1500m/25C/830hPa/rh50` +25000m.density: 0.038241271431307705 -> 0.03920390299817347 (+2.52e-02)
/// - `get_local_atmosphere_humid/1500m/25C/830hPa/rh100` +25000m.density: 0.037964907373655 -> 0.03990031371536156 (+5.10e-02)
#[rustfmt::skip]
const LOCAL_AIR_HUMID: Pinned = &[
    (
        "get_local_atmosphere_humid/0m/15C/1013.25hPa/rh37.5",
        &[
            0x3fefeff8156b0190, 0x3ff3f05b03835ee2, 0x4075566ef3fa83bd, 0x3ff38fc7f384bdc2,
            0x40754a270eebb447, 0x3ff38afd4bc61414, 0x40754989aec9bf48, 0x3ff2a47a13589311,
            0x40752b55060f3fc5, 0x3fed1098bcca2451, 0x40748e991dc52527, 0x3fd4010a19e63065,
            0x407275d626af357f, 0x3fa4aa111ffe557a, 0x4072ab009d8c5c79,
        ],
    ),
    (
        "get_local_atmosphere_humid/0m/15C/1013.25hPa/rh50",
        &[
            0x3fefe978c0b95f99, 0x3ff3ec33e800b968, 0x407558402c591d60, 0x3ff38bcd293e575b,
            0x40754bf73b88949e, 0x3ff38704b1486637, 0x40754b59cdff6554, 0x3ff2a0e9466d6194,
            0x40752d2292b33d9c, 0x3fed0cde9bece491, 0x40749059512b6a46, 0x3fd403bac25e8ff7,
            0x40727768a32d0621, 0x3fa4b56e5804932a, 0x4072ac97a1347924,
        ],
    ),
    (
        "get_local_atmosphere_humid/0m/15C/1013.25hPa/rh100",
        &[
            0x3fefcf861a6fc551, 0x3ff3db9e895f7ac8, 0x40755f820d3cdb32, 0x3ff37be889cadc0e,
            0x40755334ef1fc968, 0x3ff37728ca71559a, 0x407552974c0fd5aa, 0x3ff292ab688962f0,
            0x40753455cabd120d, 0x3fecfdf8901cfd84, 0x407497573a4bf840, 0x3fd40e7f69cb7981,
            0x40727daffc321698, 0x3fa4e3202f83c233, 0x4072b2f10f67ac52,
        ],
    ),
    (
        "get_local_atmosphere_humid/0m/35C/1013.25hPa/rh37.5",
        &[
            0x3fedb1f18bd4c298, 0x3ff28389cc89a762, 0x40761cf372081682, 0x3ff23030bf3f1d97,
            0x4076110c17cdd4b6, 0x3ff22c0d7eeac451, 0x407610738fd38a3e, 0x3ff164bb3c9ecdaa,
            0x4075f32e5037b4b7, 0x3feb9eddbe5ea7f5, 0x40755b79a6138291, 0x3fd4ca5d8a812fb4,
            0x4073567b116fa677, 0x3fa9c74a5451ecd9, 0x4073897db6843eec,
        ],
    ),
    (
        "get_local_atmosphere_humid/0m/35C/1013.25hPa/rh50",
        &[
            0x3fed9deca93c544a, 0x3ff276c802419ef7, 0x4076231f9086fa7d, 0x3ff223edc141c06d,
            0x40761734e3b44084, 0x3ff21fd0c418a888, 0x4076169c3126a6d5, 0x3ff159a87a618044,
            0x4075f94ec5fb287d, 0x3feb92a294e29913, 0x4075616fc383638e, 0x3fd4d1ca94385b21,
            0x40735be0e05f685d, 0x3fa9f0f77e5eef52, 0x40738ef1c2762ee3,
        ],
    ),
    (
        "get_local_atmosphere_humid/0m/35C/1013.25hPa/rh100",
        &[
            0x3fed4e1f5710466d, 0x3ff243f0625a6ebc, 0x40763bae93c1a68f, 0x3ff1f30ccbb9f7f0,
            0x40762fb6ae909fdf, 0x3ff1ef08a2b1cb1e, 0x40762f1d529c9cdc, 0x3ff12d7e061864c1,
            0x407611af657e1005, 0x3feb61afd7e8ddf9, 0x40757927e74bc623, 0x3fd4ef6b4bbc52da,
            0x4073715ad8864aa0, 0x3faa9a17c30d6138, 0x4073a4a461740e3a,
        ],
    ),
    (
        "get_local_atmosphere_humid/1500m/25C/830hPa/rh37.5",
        &[
            0x3fe93413b93edc70, 0x3fef7270d36c9751, 0x4075bb027ddc6474, 0x3feedfcb5c869ad7,
            0x4075aeed3440992e, 0x3feed88427f88bb7, 0x4075ae525e3a4c80, 0x3fed7a228607d0b3,
            0x4075909b38d608ec, 0x3fe733dbc9c36f59, 0x4074f682861ff1bc, 0x3fd03a8326719097,
            0x40734ec312a3a1a4, 0x3fa3fbf1c85106fa, 0x4073911ea0099e7f,
        ],
    ),
    (
        "get_local_atmosphere_humid/1500m/25C/830hPa/rh50",
        &[
            0x3fe92869b2a1e351, 0x3fef638dea4a13ee, 0x4075bf3090cb53ad, 0x3feed1817ad31cdd,
            0x4075b318f42cf8e0, 0x3feeca41d436de51, 0x4075b27e005e177b, 0x3fed6d4741d1e1e2,
            0x407594c123b39ea0, 0x3fe72cef8073e42a, 0x4074fa8accd3ce53, 0x3fd03edce92d5ff9,
            0x40735279d6c0cda2, 0x3fa41288b2816364, 0x407394e227cbbf31,
        ],
    ),
    (
        "get_local_atmosphere_humid/1500m/25C/830hPa/rh100",
        &[
            0x3fe8f9de745ca092, 0x3fef282964ed2244, 0x4075cfd940e49d0e, 0x3fee987d4e8b1180,
            0x4075c3b860e9c586, 0x3fee915bb0ad0649, 0x4075c31cf667c74b, 0x3fed39f4e317a3af,
            0x4075a54951fc4f79, 0x3fe7113be06cbd91, 0x40750a9cd923a6e1, 0x3fd0503f48c1597f,
            0x4073614709114f94, 0x3fa46dd05d03e687, 0x4073a3e23906c284,
        ],
    ),
];

/// NWS pressure altitude, the FAA-rule and density-matched density altitudes, the density-altitude
/// entry inverse, and the `atmosphere.density_altitude` service.
///
/// IMMUTABLE: nothing in 0.50.0 may change it, the aloft vapour model included.
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

/// Muzzle Miller Sg (the published dry T/P correction, kept deliberately; a humid correction was
/// rejected), Litz spin drift, the fast path's aerodynamic-jump offset, and `TrajectorySolver`'s
/// aerodynamic-jump components, which depend on station T/P only and so are pinned at every RH.
///
/// IMMUTABLE: nothing in 0.50.0 may change it, the aloft vapour model included.
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
/// `fast_integrate_with_segments`. 0.50.0 migrates the tuple to the uniform `KernelAtmosphere`; the
/// values do not move.
///
/// IMMUTABLE: nothing in 0.50.0 may change it, the aloft vapour model included.
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
/// Carlo, at both stations.
///
/// IMMUTABLE: nothing in 0.50.0 may change it, the aloft vapour model included.
///
/// Compared to a relative 1e-15 until 0.50.0, which made it bits like every other group: a
/// last-ulp movement was thought possible here when 0.50.0 moved the kernels onto `air_at`, but
/// none happened.
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
/// IMMUTABLE: nothing in 0.50.0 may change it, the aloft vapour model included.
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
/// MAY CHANGE ONLY WITH THE ALOFT VAPOUR MODEL: byte-identical through every other change in
/// 0.50.0.
///
/// Re-pinned by the aloft vapour model (0.50.0): the unified aloft vapour model holds the
/// governing air's (the station's or the zone's) water-vapour mole fraction for the lapse, the
/// density and the sound speed, with a virtual-temperature lapse. 202 of 638 values moved, all in
/// cases with a humidity above 0 and all off the station's own height; every RH-0 and
/// station-height value is unchanged. The largest float shift of each moved case,
/// pinned before the aloft vapour model -> pinned now:
/// - `solve_v1/rk4/icao-sea-level/rh0.5` samples[6].windage_m: -0.1662096260304339 -> -0.1662097665282316 (+8.45e-07)
/// - `solve_v1/rk45/1500m-25C-830hPa/rh0.9` samples[1].windage_m: -0.0029060136027916706 -> -0.0029060236587567126 (+3.46e-06)
/// - `card.come_ups/icao-sea-level/rh50` rows[0].drop_adj: 4.782293688882372e-06 -> 4.782677659831119e-06 (+8.03e-05)
/// - `card.come_ups/icao-sea-level/rh100` rows[0].drop_adj: -8.413634726556538e-05 -> -8.41355824696185e-05 (-9.09e-06)
/// - `card.come_ups/1500m-77F-24.51inHg/rh50` rows[7].energy: 1040.0857480660966 -> 1040.0876549435054 (+1.83e-06)
/// - `card.come_ups/1500m-77F-24.51inHg/rh100` rows[7].energy: 1046.994806293369 -> 1046.9985981368106 (+3.62e-06)
#[rustfmt::skip]
const FRONT_ENDS_HUMID: Pinned = &[
    (
        "solve_v1/rk4/icao-sea-level/rh0.5",
        &[
            0x4082c00000000000, 0x3ff43727c932ce97, 0x3fee41cf404ddc15, 0x407e9f35f04051ba,
            0x4095445f20dceda9, 0x40045740e6bb5acd, 0x3fbb74e1c1bee04b, 0x3fdb6fe38aae2234,
            0x3f54b33333333334, 0x000000000000000d, 0x4049000000000000, 0x3fafb76091fb121e,
            0x4088bcfeda5cc7a5, 0xbf4b462c63800800, 0xbf708aa603d48622, 0x400295e85ccea387,
            0x4072c00000000000, 0x3fda5d20f08efe99, 0x4084274f17b2d2f8, 0x3fd9b90a43838afa,
            0xbfc5465c93c1dcc8, 0x3ffe4842b2b4d6e1, 0x4082c00000000000, 0x3fee41cf404ddc15,
            0x407e9f35f04051ba, 0x400779d3688145ea, 0xbfe8476b4d1fa747, 0x3ff7018291bf9f9e,
            0x5b3ae58909d56458,
        ],
    ),
    (
        "solve_v1/rk45/1500m-25C-830hPa/rh0.9",
        &[
            0x4082c00000000000, 0x3ff43d3638c64072, 0x3fec805495a91df6, 0x408141d31c8131c9,
            0x409b0466160d00ab, 0x4009b189fc2b839c, 0x3fbd02e03f1dc358, 0x3fdf44c41c7ffa25,
            0x3f54933333333334, 0x000000000000000d, 0x4049000000000000, 0x3faf965ad85b383e,
            0x4088f0e3cd09662a, 0xbf572c7db87b8000, 0xbf67ce5f926ca527, 0x4002589714d18fdc,
            0x4072c00000000000, 0x3fd9ab846d256bc5, 0x40854070103dc12a, 0x3fd7c8b0ad61b8f2,
            0xbfbeb676d0c4de68, 0x3fff43be05d0bae9, 0x4082c00000000000, 0x3fec805495a91df6,
            0x408141d31c8131c9, 0x40050f5515bc3174, 0xbfe10fc0bd8968e0, 0x3ff963512bd3710c,
            0x87de56d04dbb94f3,
        ],
    ),
    (
        "card.come_ups/icao-sea-level/rh50",
        &[
            0x0000000000000008, 0x0000000000000001, 0x4059000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x3ed40f5c6cba9831, 0x0000000000000001,
            0x3ed40f5c6cba9831, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a3a2be8dda1d70, 0x0000000000000001,
            0x40a32c7d28f874c9, 0x0000000000000001, 0x3fbd7bdce124fd68, 0x0000000000000001,
            0x4069000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x3fe1d04849521e4f, 0x0000000000000001, 0x3fe1d03e41a3e7f2, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40a23c28133d327d, 0x0000000000000001, 0x40a08929d8e99274, 0x0000000000000001,
            0x3fce98b93f459dd0, 0x0000000000000001, 0x4072c00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x3ff557d4d0f4d432, 0x0000000000000001,
            0x3fe8df6158978a15, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a0e4f75470e59e, 0x0000000000000001,
            0x409c633e0075a125, 0x0000000000000001, 0x3fd7d83712f72046, 0x0000000000000001,
            0x4079000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4001ec19609b1fcb, 0x0000000000000001, 0x3fed00bbe082d6c8, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x409f3a15a87769a6, 0x0000000000000001, 0x40983effd683dcba, 0x0000000000000001,
            0x3fe08a5eed4831c4, 0x0000000000000001, 0x407f400000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x400a25a9c65cccc7, 0x0000000000000001,
            0x3ff07320cb8359f8, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x409cc4e5f351da04, 0x0000000000000001,
            0x409494633caf4807, 0x0000000000000001, 0x3fe58b866ec1a8cd, 0x0000000000000001,
            0x4082c00000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4011b8377270b67e, 0x0000000000000001, 0x3ff2958a3d09406a, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x409a676fe2d1cf02, 0x0000000000000001, 0x409155c6650542e7, 0x0000000000000001,
            0x3feafccfabeebf39, 0x0000000000000001, 0x4085e00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x4016fa03449cf2a5, 0x0000000000000001,
            0x3ff5072f48b0f09c, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x409820309537fbfe, 0x0000000000000001,
            0x408cf1f4169e9a1d, 0x0000000000000001, 0x3ff0772bc4fbb68e, 0x0000000000000001,
            0x4089000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x401cf38238315d0d, 0x0000000000000001, 0x3ff7e5fbce51a9a0, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4095efdb1cc6f39e, 0x0000000000000001, 0x4087ee69373f8da3, 0x0000000000000001,
            0x3ff3b9dd360adc2a,
        ],
    ),
    (
        "card.come_ups/icao-sea-level/rh100",
        &[
            0x0000000000000008, 0x0000000000000001, 0x4059000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0xbf160e3e4cf0b181, 0x0000000000000001,
            0xbf160e3e4cf0b181, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a3a3aef981d7de, 0x0000000000000001,
            0x40a32e52bba920de, 0x0000000000000001, 0x3fbd7b280cc17344, 0x0000000000000001,
            0x4069000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x3fe1cd3231189393, 0x0000000000000001, 0x3fe1cde2a30afb19, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40a23df852caec43, 0x0000000000000001, 0x40a08c73fcbb2f1f, 0x0000000000000001,
            0x3fce9734a419d5ba, 0x0000000000000001, 0x4072c00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x3ff553e0632bbbb5, 0x0000000000000001,
            0x3fe8da8e953ee3d7, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a0e7960bb36ab8, 0x0000000000000001,
            0x409c6c0cb1566c45, 0x0000000000000001, 0x3fd7d65cf6c3c6e2, 0x0000000000000001,
            0x4079000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4001e847c5df9116, 0x0000000000000001, 0x3fecf95e5126ccee, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x409f40c35edd75ee, 0x0000000000000001, 0x4098495ff4ce2275, 0x0000000000000001,
            0x3fe088984d818959, 0x0000000000000001, 0x407f400000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x400a1efccb29daa8, 0x0000000000000001,
            0x3ff06d6a0a949324, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x409cccdabb0f5b17, 0x0000000000000001,
            0x40949fc6c80bba70, 0x0000000000000001, 0x3fe588846bcb99b0, 0x0000000000000001,
            0x4082c00000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4011b2e05fd35490, 0x0000000000000001, 0x3ff28d87e8f99cf0, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x409a708bdc89ff4f, 0x0000000000000001, 0x409161be7b9827d6, 0x0000000000000001,
            0x3feaf81c5bebaf5d, 0x0000000000000001, 0x4085e00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x4016f1f7a2b74dd0, 0x0000000000000001,
            0x3ff4fc5d0b8fe500, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40982a578f6a6caf, 0x0000000000000001,
            0x408d0a5584cbb9ac, 0x0000000000000001, 0x3ff073b11ff1259e, 0x0000000000000001,
            0x4089000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x401ce79874d1db40, 0x0000000000000001, 0x3ff7d683486a35c0, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4095faec28824b11, 0x0000000000000001, 0x40880694851dec98, 0x0000000000000001,
            0x3ff3b4de7f67c462,
        ],
    ),
    (
        "card.come_ups/1500m-77F-24.51inHg/rh50",
        &[
            0x0000000000000008, 0x0000000000000001, 0x4059000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x3f41f711ff136095, 0x0000000000000001,
            0x3f41f711ff136095, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a3edff4092017b, 0x0000000000000001,
            0x40a3c08dc25e7ecd, 0x0000000000000001, 0x3fbd43ce229f0678, 0x0000000000000001,
            0x4069000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x3fe11a3db9b26f02, 0x0000000000000001, 0x3fe115bff532aa2a, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40a2cd4f41a40c77, 0x0000000000000001, 0x40a19482c4f081bc, 0x0000000000000001,
            0x3fce2172122d06f5, 0x0000000000000001, 0x4072c00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x3ff44bc2dd84c322, 0x0000000000000001,
            0x3fe77d4801571742, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a1b6498316011d, 0x0000000000000001,
            0x409f33b597641e35, 0x0000000000000001, 0x3fd74900dbceb5fd, 0x0000000000000001,
            0x4079000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4000d891d8a1bd03, 0x0000000000000001, 0x3feacac1a77d6dc8, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40a0a935830e8cb0, 0x0000000000000001, 0x409b9bcf4d8c76a6, 0x0000000000000001,
            0x3fe002386ad29454, 0x0000000000000001, 0x407f400000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x4008445386c61b0d, 0x0000000000000001,
            0x3fedaf06b8917828, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x409f4ab3fb9f2a7f, 0x0000000000000001,
            0x409858d53ce2b6b0, 0x0000000000000001, 0x3fe4a777f8ae5439, 0x0000000000000001,
            0x4082c00000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x401036b3da4e7e40, 0x0000000000000001, 0x3ff052285badc2e6, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x409d5351d0c9876f, 0x0000000000000001, 0x4095621dce06302f, 0x0000000000000001,
            0x3fe99b2129bef7ec, 0x0000000000000001, 0x4085e00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x4014b29445c924c6, 0x0000000000000001,
            0x3ff1ef81adea9a18, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x409b6afc0f107439, 0x0000000000000001,
            0x4092b11ea0a2995d, 0x0000000000000001, 0x3feee57117519106, 0x0000000000000001,
            0x4089000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4019a17130c3afcb, 0x0000000000000001, 0x3ff3bb73abea2c14, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x409990bf34d2b4c2, 0x0000000000000001, 0x40904059c237aec1, 0x0000000000000001,
            0x3ff247f8288fe5ef,
        ],
    ),
    (
        "card.come_ups/1500m-77F-24.51inHg/rh100",
        &[
            0x0000000000000008, 0x0000000000000001, 0x4059000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x3f396a57413995d3, 0x0000000000000001,
            0x3f396a57413995d3, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a3efb7bf81a48b, 0x0000000000000001,
            0x40a3c3f70d1c5944, 0x0000000000000001, 0x3fbd428a323fb35b, 0x0000000000000001,
            0x4069000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x3fe114b6296d1713, 0x0000000000000001, 0x3fe11188de84efe0, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40a2d0a8b57fcfe7, 0x0000000000000001, 0x40a19ac6c98204a6, 0x0000000000000001,
            0x3fce1ec35582ee07, 0x0000000000000001, 0x4072c00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x3ff44505e5b9b63d, 0x0000000000000001,
            0x3fe77555a2065567, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40a1bb2bedddb8de, 0x0000000000000001,
            0x409f44ed7e74c4b2, 0x0000000000000001, 0x3fd745cb6d56be65, 0x0000000000000001,
            0x4079000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4000d23e63d6d563, 0x0000000000000001, 0x3feabeedc3e7e912, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40a0af82c35676be, 0x0000000000000001, 0x409bb0b5d6825a58, 0x0000000000000001,
            0x3fdffe643bd1c466, 0x0000000000000001, 0x407f400000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x40083976bbfb9756, 0x0000000000000001,
            0x3fed9ce1609307cc, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x409f59e0fcf35e2b, 0x0000000000000001,
            0x4098707885ee3564, 0x0000000000000001, 0x3fe4a26f04663219, 0x0000000000000001,
            0x4082c00000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40102e4e26d7db83, 0x0000000000000001, 0x3ff0464b23683f60, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x409d64d96c01c482, 0x0000000000000001, 0x40957bb5e804bb53, 0x0000000000000001,
            0x3fe9936e9e53309b, 0x0000000000000001, 0x4085e00000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x4014a60a1e4407ee, 0x0000000000000001,
            0x3ff1deefddb0b1ac, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000001, 0x409b7eaec0cc1ccd, 0x0000000000000001,
            0x4092cc03e2466fe1, 0x0000000000000001, 0x3feeda40d018992b, 0x0000000000000001,
            0x4089000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x40198fad7dc8d672, 0x0000000000000001, 0x3ff3a68d7e133a10, 0x0000000000000000,
            0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001,
            0x4099a6745627a76c, 0x0000000000000001, 0x40905bfe90828dcb, 0x0000000000000001,
            0x3ff2402b43d8e4fb,
        ],
    ),
];

/// The bridge `solve` command (`bridge_call`) on the `FRONT_ENDS_DRY` requests.
///
/// IMMUTABLE: nothing in 0.50.0 may change it, the aloft vapour model included.
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
/// MAY CHANGE ONLY WITH THE ALOFT VAPOUR MODEL: byte-identical through every other change in
/// 0.50.0.
///
/// Re-pinned by the aloft vapour model (0.50.0): the unified aloft vapour model holds the
/// governing air's (the station's or the zone's) water-vapour mole fraction for the lapse, the
/// density and the sound speed, with a virtual-temperature lapse. 42 of 58 values moved, all in
/// cases with a humidity above 0 and all off the station's own height; every RH-0 and
/// station-height value is unchanged. The largest float shift of each moved case,
/// pinned before the aloft vapour model -> pinned now:
/// - `bridge-solve/rk4/icao-sea-level/rh0.5` samples[6].windage_m: -0.1662096260304339 -> -0.1662097665282316 (+8.45e-07)
/// - `bridge-solve/rk45/1500m-25C-830hPa/rh0.9` samples[1].windage_m: -0.0029060136027916706 -> -0.0029060236587567126 (+3.46e-06)
#[cfg(feature = "bridge")]
#[rustfmt::skip]
const BRIDGE_SOLVE_HUMID: Pinned = &[
    (
        "bridge-solve/rk4/icao-sea-level/rh0.5",
        &[
            0x4082c00000000000, 0x3ff43727c932ce97, 0x3fee41cf404ddc15, 0x407e9f35f04051ba,
            0x4095445f20dceda9, 0x40045740e6bb5acd, 0x3fbb74e1c1bee04b, 0x3fdb6fe38aae2234,
            0x3f54b33333333334, 0x000000000000000d, 0x4049000000000000, 0x3fafb76091fb121e,
            0x4088bcfeda5cc7a5, 0xbf4b462c63800800, 0xbf708aa603d48622, 0x400295e85ccea387,
            0x4072c00000000000, 0x3fda5d20f08efe99, 0x4084274f17b2d2f8, 0x3fd9b90a43838afa,
            0xbfc5465c93c1dcc8, 0x3ffe4842b2b4d6e1, 0x4082c00000000000, 0x3fee41cf404ddc15,
            0x407e9f35f04051ba, 0x400779d3688145ea, 0xbfe8476b4d1fa747, 0x3ff7018291bf9f9e,
            0x5b3ae58909d56458,
        ],
    ),
    (
        "bridge-solve/rk45/1500m-25C-830hPa/rh0.9",
        &[
            0x4082c00000000000, 0x3ff43d3638c64072, 0x3fec805495a91df6, 0x408141d31c8131c9,
            0x409b0466160d00ab, 0x4009b189fc2b839c, 0x3fbd02e03f1dc358, 0x3fdf44c41c7ffa25,
            0x3f54933333333334, 0x000000000000000d, 0x4049000000000000, 0x3faf965ad85b383e,
            0x4088f0e3cd09662a, 0xbf572c7db87b8000, 0xbf67ce5f926ca527, 0x4002589714d18fdc,
            0x4072c00000000000, 0x3fd9ab846d256bc5, 0x40854070103dc12a, 0x3fd7c8b0ad61b8f2,
            0xbfbeb676d0c4de68, 0x3fff43be05d0bae9, 0x4082c00000000000, 0x3fec805495a91df6,
            0x408141d31c8131c9, 0x40050f5515bc3174, 0xbfe10fc0bd8968e0, 0x3ff963512bd3710c,
            0x87de56d04dbb94f3,
        ],
    ),
];

/// The FFI single solve `ballistics_calculate_trajectory` (RK4, sampling and spin drift on) at both
/// stations, RH 0.
///
/// IMMUTABLE: nothing in 0.50.0 may change it, the aloft vapour model included.
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
/// MAY CHANGE ONLY WITH THE ALOFT VAPOUR MODEL: byte-identical through every other change in
/// 0.50.0.
///
/// Re-pinned by the aloft vapour model (0.50.0): the unified aloft vapour model holds the
/// governing air's (the station's or the zone's) water-vapour mole fraction for the lapse, the
/// density and the sound speed, with a virtual-temperature lapse. 48 of 68 values moved, all in
/// cases with a humidity above 0 and all off the station's own height; every RH-0 and
/// station-height value is unchanged. The largest float shift of each moved case,
/// pinned before the aloft vapour model -> pinned now:
/// - `ffi-trajectory/rk4/icao-sea-level/rh50` impact_energy: 1360.9028469099487 -> 1360.9033619221536 (+3.78e-07)
/// - `ffi-trajectory/rk4/icao-sea-level/rh100` impact_energy: 1365.0156018568925 -> 1365.0166285531138 (+7.52e-07)
/// - `ffi-trajectory/rk4/1500m-25C-830hPa/rh50` impact_energy: 1721.766284074095 -> 1721.7671717118649 (+5.16e-07)
/// - `ffi-trajectory/rk4/1500m-25C-830hPa/rh100` impact_energy: 1730.7141932216796 -> 1730.7159555847395 (+1.02e-06)
#[cfg(feature = "ffi")]
#[rustfmt::skip]
const FFI_TRAJECTORY_HUMID: Pinned = &[
    (
        "ffi-trajectory/rk4/icao-sea-level/rh50",
        &[
            0x4082c00000000000, 0x3fc0888f190e8a73, 0x3fee421041989671, 0x407e9eaa343db311,
            0x4095439d0ae86068, 0x7ff8000000000000, 0x00000000000003b3, 0x3fee421041989671,
            0x4082c00000000000, 0xc00423ec34f14817, 0xbfebb70533b0d645, 0x407e9eaa343db311,
            0xad877748067969bb, 0x000000000000000d, 0x04d5dd856d8048b5, 0x40842724cab5a084,
            0x3ffe48032350300b,
        ],
    ),
    (
        "ffi-trajectory/rk4/icao-sea-level/rh100",
        &[
            0x4082c00000000000, 0x3fc0894039dc54ee, 0x3fee3c2d4aa317cf, 0x407eaa807de0731d,
            0x4095541107134f37, 0x7ff8000000000000, 0x00000000000003b2, 0x3fee3c2d4aa317cf,
            0x4082c00000000000, 0xc0041bd3e9b58544, 0xbfeb9f7957dadc37, 0x407eaa807de0731d,
            0xa02b09170e06e842, 0x000000000000000d, 0xc97b49820f9c8e06, 0x40842a988dbb5cca,
            0x3ffe42e90fcbfd0c,
        ],
    ),
    (
        "ffi-trajectory/rk4/1500m-25C-830hPa/rh50",
        &[
            0x4082c00000000000, 0x3fc0c0a6a6cb2acd, 0x3fec882049f3bb49, 0x408138728ca11337,
            0x409ae71195761381, 0x7ff8000000000000, 0x000000000000037d, 0x3fec882049f3bb49,
            0x4082c00000000000, 0xc001cf6349c3e26c, 0xbfe4cf45551d6a1b, 0x408138728ca11337,
            0x57dd55a6a2010197, 0x000000000000000d, 0xa27d7c359e7e9931, 0x40853b23474ca883,
            0x3fff4f1901d03624,
        ],
    ),
    (
        "ffi-trajectory/rk4/1500m-25C-830hPa/rh100",
        &[
            0x4082c00000000000, 0x3fc0c1f233818054, 0x3fec7e93eed79e7a, 0x408143e391368f42,
            0x409b0add2375f761, 0x7ff8000000000000, 0x000000000000037c, 0x3fec7e93eed79e7a,
            0x4082c00000000000, 0xc001c2b0f59b2908, 0xbfe4a913e8acf6bd, 0x408143e391368f42,
            0x0a8536724ab18b02, 0x000000000000000d, 0xe85be4f9652fdf36, 0x408541984bf62dd9,
            0x3fff40add86fb5aa,
        ],
    ),
];
