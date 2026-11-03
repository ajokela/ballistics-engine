//! The public docs say what the humidity code does: no stale comment, and no doc that
//! contradicts the code.
//!
//! Every number the docs give for the humidity change is recomputed here by running the engine,
//! so a doc can neither quote an estimate nor drift from the code. Each is checked IN PLACE: a
//! table must hold exactly the rows the engine gives, whole (so each cell is held to its own
//! quantity, not to the same digits found anywhere in the section), and a figure in prose is
//! matched with the words around it. "0.49.0" is the released 0.49.0 crate itself, run beside
//! this one (the `ballistics_engine_049` dev-dependency):
//!
//! - **Sound speed.** The kernels' station sound speed in 0.49.0 (dry, 0 %) and now (the
//!   station's own humidity), for every row of the table, and the change.
//! - **Trajectories.** `try_integrate_trajectory` RK4 at 1 ms: a 2 s flat shot and a 30-degree,
//!   30 s lob, on 0.49.0 and 0.50.0; the `fast_integrate` flat shot the prose quotes; and the lean
//!   Monte Carlo at 500 m and 1000 m. [`legacy_kernel_air_reproduces_0_49_0`] holds
//!   `KernelAtmosphere::from_legacy_params` with `RelativeHumidity::DRY` to 0.49.0's kernel bit
//!   for bit, and [`the_legacy_recipe_is_exact_only_without_humid_zones`] shows that holds only
//!   without humid `AtmoSock` zones, as both documents now say.
//! - **The aloft vapour model.** The air at +1 km and +3 km, and three `TrajectorySolver`
//!   shots (flat to 1000 m, a 30-degree lob, and both from a 1500 m station), 0.49.0 against
//!   0.50.0. [`the_air_before_0_50_0_is_the_old_aloft_model`] also rebuilds the old air from
//!   the old formulas. The aloft change's own share, which only the development build before it
//!   could measure, is labelled as a one-time measurement, as is the CLI golden-corpus count.
//! - **The Miller exception.** The factor a humid density would put on the muzzle Miller Sg,
//!   which 0.50.0 deliberately does not apply, in every document that quotes it, and the spin
//!   drift and aerodynamic jump it would move. A pin shows the solve's Sg does not move with
//!   humidity.
//! - **0.49.0's behaviour.** The deprecated `get_local_atmosphere`'s sound-speed shift and the
//!   steam a NaN humidity flew as.
//! - **Constants.** The engine's own dry sound speed at 15 °C, which `SPEED_OF_SOUND_MPS`'s doc
//!   states.
//! - **The model's own limits.** The held vapour's supersaturation aloft, and the gaps between
//!   the lapse's, the density's and the sound speed's gas constants, at every site that states
//!   them.
//! - **Disclosures.** A kernel caller's new height reference, the breaking list and the new
//!   items in both documents, the default and its exceptions, the sensor band on every public
//!   humidity field, and no private app's internals in the public release notes.
//!
//! It also checks the CHANGELOG's 0.50.0 section carries the breaking API list, the sensor band
//! from its constant and the downstream follow-ups, and that the stale wording is gone. The section
//! is found under `## [Unreleased]` before the release cut and under `## [0.50.0] - <date>` after
//! it, and a test applies `scripts/release/bump.sh`'s cut in memory, so cutting 0.50.0 cannot turn
//! this suite red.

use std::fs;
use std::path::Path;

use ballistics_engine::atmosphere::{calculate_air_density_cimp, moist_speed_of_sound};
use ballistics_engine::trajectory_integration::{try_integrate_trajectory, TrajectoryParams};
use ballistics_engine::{
    air_at, decode_solve_request_v1, solve_v1, BallisticInputs, DragModel, KernelAtmosphere,
    LocalAir, RelativeHumidity, StationAir, HUMIDITY_SENSOR_TOLERANCE_PERCENT,
};
/// The released 0.49.0 engine, run beside this one: the "before" of every 0.49.0 -> 0.50.0
/// figure.
#[cfg(not(target_arch = "wasm32"))]
use ballistics_engine_049 as old;

fn read(relative: &str) -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join(relative);
    fs::read_to_string(&path).unwrap_or_else(|error| panic!("read {}: {error}", path.display()))
}

/// The release whose CHANGELOG section holds the humidity notes.
const HUMIDITY_RELEASE: &str = "0.50.0";
/// The release just before it.
const PREVIOUS_RELEASE: &str = "0.49.0";

/// The CHANGELOG section that holds the 0.50.0 humidity notes: its heading line, then its body up
/// to the next `## [` heading, split the way `scripts/release/extract-notes.sh` splits them.
///
/// Before the release cut the section is `## [Unreleased]`. `scripts/release/bump.sh 0.50.0`
/// renames it `## [0.50.0] - <date>`, and work for the next release later opens a new
/// `## [Unreleased]` above it. So the 0.50.0 heading, with any date, wins when it exists, and the
/// `## [Unreleased]` section is the fallback. Either way the section must sit directly above
/// `## [0.49.0]`.
fn release_section() -> String {
    release_section_in(&read("CHANGELOG.md"))
}

/// [`release_section`] on a CHANGELOG text.
fn release_section_in(changelog: &str) -> String {
    let starts: Vec<usize> = std::iter::once(0)
        .chain(changelog.match_indices('\n').map(|(at, _)| at + 1))
        .filter(|&at| changelog[at..].starts_with("## ["))
        .collect();
    let heading = |at: usize| changelog[at..].lines().next().unwrap_or_default();
    let released = format!("## [{HUMIDITY_RELEASE}]");
    let index = starts
        .iter()
        .position(|&at| heading(at).starts_with(&released))
        .or_else(|| {
            starts
                .iter()
                .position(|&at| heading(at).starts_with("## [Unreleased]"))
        })
        .unwrap_or_else(|| {
            panic!("CHANGELOG.md has neither a `{released}` nor a `## [Unreleased]` section")
        });
    let previous = format!("## [{PREVIOUS_RELEASE}]");
    let next = starts.get(index + 1).copied();
    assert!(
        next.is_some_and(|at| heading(at).starts_with(&previous)),
        "CHANGELOG.md: `{}` must sit directly above `{previous}`, not above `{}`",
        heading(starts[index]),
        next.map_or("the end of the file", heading)
    );
    changelog[starts[index]..next.unwrap_or(changelog.len())].to_string()
}

/// A section without its heading line.
fn body_of(section: &str) -> &str {
    section.split_once('\n').map_or("", |(_, body)| body)
}

/// `scripts/release/bump.sh`'s CHANGELOG cut, in Rust: the first `## [Unreleased]` becomes
/// `## [<new>] - <date>`. The test below checks the script still makes exactly this replace.
fn bump_cut(changelog: &str, new: &str, date: &str) -> String {
    changelog.replacen("## [Unreleased]", &format!("## [{new}] - {date}"), 1)
}

/// Whether the humidity notes still wait for the release cut, under `## [Unreleased]`.
fn awaits_the_cut(changelog: &str) -> bool {
    release_section_in(changelog).starts_with("## [Unreleased]")
}

/// The CHANGELOG once 0.50.0 is cut: `bump.sh`'s cut while the notes await it, the text as it
/// stands once they sit under their dated heading.
fn as_cut(changelog: &str, date: &str) -> String {
    if awaits_the_cut(changelog) {
        bump_cut(changelog, HUMIDITY_RELEASE, date)
    } else {
        changelog.to_string()
    }
}

/// The release cut must not lose the humidity notes: on the CHANGELOG that `bump.sh 0.50.0`
/// writes, the lookup finds the same section body, under its dated `## [0.50.0]` heading. Once
/// the cut has happened, the notes must sit under that dated heading.
#[test]
fn the_release_cut_keeps_the_0_50_0_section() {
    let script = read("scripts/release/bump.sh");
    for line in [
        "if '## [Unreleased]' not in s:",
        "s = s.replace('## [Unreleased]', f'## [{new}] - {date}', 1)",
        "DATE=$(date +%Y-%m-%d)",
    ] {
        assert_contains(&script, line, "scripts/release/bump.sh's CHANGELOG cut");
    }
    let changelog = read("CHANGELOG.md");
    let before = release_section_in(&changelog);
    if !awaits_the_cut(&changelog) {
        assert!(
            before.starts_with(&format!("## [{HUMIDITY_RELEASE}] - ")),
            "the cut notes sit under `{}`",
            before.lines().next().unwrap_or_default()
        );
        return;
    }
    for date in ["2026-11-02", "2027-01-31"] {
        let after = release_section_in(&bump_cut(&changelog, HUMIDITY_RELEASE, date));
        assert_eq!(
            after.lines().next(),
            Some(format!("## [{HUMIDITY_RELEASE}] - {date}").as_str()),
            "the cut section's heading"
        );
        assert_eq!(body_of(&after), body_of(&before), "the cut section's body");
    }
}

/// After 0.50.0 ships, work for the next release opens a new `## [Unreleased]` above it; the
/// lookup must still find the 0.50.0 notes, not the new section. This holds before the cut (on
/// the CHANGELOG `bump.sh` would write) and after it (on the CHANGELOG as it stands).
#[test]
fn a_later_unreleased_section_does_not_hide_the_0_50_0_notes() {
    let cut = as_cut(&read("CHANGELOG.md"), "2026-11-02");
    let notes = release_section_in(&cut);
    assert!(
        notes.starts_with(&format!("## [{HUMIDITY_RELEASE}] - ")),
        "the cut notes sit under `{}`",
        notes.lines().next().unwrap_or_default()
    );
    let at = cut.find(&notes).expect("the notes sit in the CHANGELOG");
    let reopened = format!(
        "{}## [Unreleased]\n\n### Added\n- Something for the next release.\n\n{}",
        &cut[..at],
        &cut[at..]
    );
    assert_eq!(
        release_section_in(&reopened),
        notes,
        "the lookup must return the 0.50.0 notes, not the new section"
    );
}

/// The section becomes the GitHub release body as it stands (`scripts/release/extract-notes.sh`),
/// and the cut renames only its heading, so its prose must already read as released.
#[test]
fn the_0_50_0_section_reads_as_released() {
    let section = release_section();
    let body = body_of(&section);
    for unreleased_wording in [
        "To be released",
        "to be released",
        "Unreleased",
        "this branch",
        "the train",
    ] {
        assert!(
            !body.contains(unreleased_wording),
            "CHANGELOG 0.50.0 section: `{unreleased_wording}` reads wrong once the release is cut"
        );
    }
}

fn assert_contains(haystack: &str, needle: &str, what: &str) {
    assert!(haystack.contains(needle), "{what}: missing `{needle}`");
}

// =================================================================================================
// The CHANGELOG's 0.50.0 section
// =================================================================================================

/// Before the cut the notes sit under `## [Unreleased]`, after it under `## [0.50.0] - <date>`,
/// either way directly above 0.49.0 (`release_section` asserts the position).
#[test]
fn the_changelog_has_the_0_50_0_section_directly_above_0_49_0() {
    let section = release_section();
    let heading = section.lines().next().unwrap_or_default();
    assert!(
        heading == "## [Unreleased]" || heading.starts_with(&format!("## [{HUMIDITY_RELEASE}]")),
        "the humidity notes sit under `{heading}`"
    );
}

#[test]
fn the_0_50_0_section_lists_every_breaking_rust_api_change() {
    let section = release_section();
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
        // Typed air causes added while preserving the originating station reading.
        "CardServiceError::Atmosphere",
        "DsfServiceErrorV1::Atmosphere",
        "UncertaintyTruingErrorV1::Atmosphere",
        "WindScenarioError::Atmosphere",
        "TruingPlanErrorV1::atmosphere_error",
    ] {
        assert_contains(&section, name, "CHANGELOG 0.50.0 breaking list");
    }
}

#[test]
fn the_0_50_0_section_states_the_sensor_band_from_its_constant() {
    let section = release_section();
    let band = HUMIDITY_SENSOR_TOLERANCE_PERCENT;
    assert_contains(
        &section,
        "HUMIDITY_SENSOR_TOLERANCE_PERCENT",
        "the band constant",
    );
    assert_prose(
        &section,
        &format!(
            "`HUMIDITY_SENSOR_TOLERANCE_PERCENT` = {band} points around 0-100 %, so a percent \
             from {} to {}, or a fraction from {} to {}.",
            -band,
            100.0 + band,
            -band / 100.0,
            1.0 + band / 100.0
        ),
        "the sensor band",
    );
}

/// The follow-ups name every public repository that pins the engine, each as its own bullet, and
/// one neutral line for the apps and services that embed it: their own task lists live with
/// them, not in a public release note.
#[test]
fn the_0_50_0_section_lists_the_downstream_follow_ups() {
    let section = release_section();
    for repo in [
        "rockets-engine",
        "winds-engine",
        "ballistics-engine-py",
        "ballistics-engine-rb",
        "handbook",
    ] {
        assert_contains(&section, &format!("- **{repo}**"), "a downstream follow-up");
    }
    assert_contains(
        &section,
        "- **Apps and services that embed the engine**",
        "the embedders' follow-up",
    );
}

/// The `# Humidity Consistency (0.50.0)` section of MCCOY_MIGRATION.md, to the next `# ` heading.
fn migration_section() -> String {
    let migration = read("MCCOY_MIGRATION.md");
    let start = migration
        .find("\n# Humidity Consistency (0.50.0)")
        .expect("MCCOY_MIGRATION.md has the 0.50.0 humidity notes")
        + 1;
    let body = &migration[start..];
    let end = body.find("\n# ").unwrap_or(body.len());
    body[..end].to_string()
}

/// The CHANGELOG section ships as the GitHub release body and in the crates.io package, and the
/// migration notes ship in the package too, so neither may carry a private app's or service's
/// internals: its environment variables, deploy details, private models, defects or chores.
/// Those follow-ups live with the apps; the notes keep one neutral line for embedders.
#[test]
fn the_release_notes_carry_no_private_app_internals() {
    for (what, text) in [
        ("the CHANGELOG's 0.50.0 section", release_section()),
        ("MCCOY_MIGRATION.md's 0.50.0 section", migration_section()),
    ] {
        for private in [
            "Flask",
            "flask",
            "BALLISTICS_",
            "Dockerfile",
            "/v1/calculate",
            "ballistics_rust",
            "hardcoded",
            "the live apps",
            "neither app",
            "iOS",
            "Android",
        ] {
            assert!(
                !text.contains(private),
                "{what}: `{private}` is a private app's material, not a public release note's"
            );
        }
    }
}

// =================================================================================================
// Numbers the CHANGELOG quotes, from the engine, each in its own row and column
// =================================================================================================
//
// Every quoted number is checked in place: a table must hold exactly the rows the engine gives,
// each row whole (so a cell is checked against its own quantity, not against the same digits
// anywhere in the section), and a figure in prose is checked with the words around it. "0.49.0"
// is the released 0.49.0 crate itself, run beside this one (the `ballistics_engine_049`
// dev-dependency; native only, so these tests are too).

/// The rows of the Markdown table whose header line is `header`, in order (not the header or the
/// `|---|` rule).
fn table_rows(section: &str, header: &str) -> Vec<String> {
    let mut lines = section.lines().skip_while(|line| *line != header);
    assert_eq!(lines.next(), Some(header), "no table headed `{header}`");
    lines.next();
    lines
        .take_while(|line| line.starts_with('|'))
        .map(str::to_string)
        .collect()
}

/// The table headed `header` holds exactly `expected`, row for row.
fn assert_table(section: &str, header: &str, expected: &[String]) {
    let rows = table_rows(section, header);
    assert_eq!(
        rows, expected,
        "the table headed `{header}` is not the engine's:\n  doc:    {rows:#?}\n  engine: {expected:#?}"
    );
}

/// `text` with every run of whitespace, line breaks included, as one space, and each line's
/// leading `///`, `//!` or `//` dropped, so prose can be matched across a line wrap in a
/// Markdown file or a doc comment.
fn normalized(text: &str) -> String {
    text.lines()
        .map(|line| {
            let line = line.trim_start();
            line.strip_prefix("///")
                .or_else(|| line.strip_prefix("//!"))
                .or_else(|| line.strip_prefix("//"))
                .unwrap_or(line)
        })
        .flat_map(str::split_whitespace)
        .collect::<Vec<_>>()
        .join(" ")
}

/// `prose` (one sentence or clause, single-spaced) appears in `text` as written, line wraps aside.
fn assert_prose(text: &str, prose: &str, what: &str) {
    assert!(
        normalized(text).contains(prose),
        "{what}: missing `{prose}`"
    );
}

const SEA_LEVEL_HPA: f64 = 1013.25;
/// The humid conditions most of the CHANGELOG's figures quote, at sea level: (temperature C, RH
/// percent).
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

/// What 0.49.0's kernels and lean Monte Carlo were given for a sea-level station: its humid
/// CIPM density ratio (0.49.0's `calculate_atmosphere`), in the legacy standard-mode tuple.
#[cfg(not(target_arch = "wasm32"))]
fn released_tuple(temperature_c: f64, humidity_percent: f64) -> (f64, f64, f64, f64) {
    let (density, _) = old::atmosphere::calculate_atmosphere(
        0.0,
        Some(temperature_c),
        Some(SEA_LEVEL_HPA),
        humidity_percent,
    );
    (0.0, temperature_c, SEA_LEVEL_HPA, density / 1.225)
}

/// The kernel sound-speed table's stations, in its order: (temperature C, RH percent).
#[cfg(not(target_arch = "wasm32"))]
const SOUND_SPEED_ROWS: [(f64, f64); 5] = [
    (-10.0, 50.0),
    (15.0, 50.0),
    (15.0, 100.0),
    (30.0, 50.0),
    (35.0, 90.0),
];

#[cfg(not(target_arch = "wasm32"))]
#[test]
fn the_0_50_0_sound_speed_shifts_are_engine_output() {
    let expected: Vec<String> = SOUND_SPEED_ROWS
        .iter()
        .map(|&(temperature_c, humidity_percent)| {
            // 0.49.0's kernels flew the station's air with the sound speed at 0 %: its own
            // per-step air at the station's height, humidity argument 0.
            let (.., ratio) = released_tuple(temperature_c, humidity_percent);
            let (_, before) = old::atmosphere::get_local_atmosphere_humid(
                0.0,
                0.0,
                temperature_c,
                SEA_LEVEL_HPA,
                ratio,
                0.0,
            );
            let after = station(temperature_c, humidity_percent).speed_of_sound_mps();
            format!(
                "| {temperature_c} °C / {humidity_percent} % | {before:.3} m/s | {after:.3} m/s | \
                 {:+.3}% |",
                (after / before - 1.0) * 100.0
            )
        })
        .collect();
    assert_table(
        &release_section(),
        "| Station, 1013.25 hPa | Kernel sound speed, 0.49.0 | 0.50.0 | Change |",
        &expected,
    );
}

const MASS_KG: f64 = 0.011_34; // .308 175 gr
const DIAMETER_M: f64 = 0.007_82;
const LENGTH_M: f64 = 0.031_5; // 1.24 in
const MUZZLE_VELOCITY_MPS: f64 = 792.48; // 2600 fps
const BC_G7: f64 = 0.243;

/// `TrajectoryParams` for `engine` (this crate or 0.49.0), with its own atmosphere field.
macro_rules! kernel_params {
    ($engine:ident, $($atmosphere:tt)+) => {
        $engine::trajectory_integration::TrajectoryParams {
            mass_kg: MASS_KG,
            bc: BC_G7,
            drag_model: $engine::DragModel::G7,
            wind_segments: vec![],
            $($atmosphere)+,
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
    };
}

fn params(atmosphere: KernelAtmosphere) -> TrajectoryParams {
    kernel_params!(ballistics_engine, atmosphere)
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

fn launch(angle_rad: f64) -> [f64; 6] {
    [
        0.0,
        0.0,
        0.0,
        MUZZLE_VELOCITY_MPS * angle_rad.cos(),
        MUZZLE_VELOCITY_MPS * angle_rad.sin(),
        0.0,
    ]
}

/// `(downrange m, height m, speed m/s)` of a kernel state `[x, y, z, vx, vy, vz]`.
fn endpoint(state: &impl std::ops::Index<usize, Output = f64>) -> (f64, f64, f64) {
    let speed = (state[3] * state[3] + state[4] * state[4] + state[5] * state[5]).sqrt();
    (state[0], state[1], speed)
}

/// `(downrange m, height m, speed m/s)` at the end of an RK4 1 ms flight.
fn fly(atmosphere: KernelAtmosphere, angle_rad: f64, seconds: f64) -> (f64, f64, f64) {
    let steps = try_integrate_trajectory(
        launch(angle_rad),
        (0.0, seconds),
        params(atmosphere),
        "RK4",
        1e-6,
        0.001,
    )
    .expect("valid kernel input");
    let (time, state) = steps.last().unwrap();
    assert_eq!(*time, seconds);
    endpoint(state)
}

/// The same flight on the released 0.49.0 kernel, given what 0.49.0's callers gave it.
#[cfg(not(target_arch = "wasm32"))]
fn fly_released(
    temperature_c: f64,
    humidity_percent: f64,
    angle_rad: f64,
    seconds: f64,
) -> (f64, f64, f64) {
    let steps = old::trajectory_integration::try_integrate_trajectory(
        launch(angle_rad),
        (0.0, seconds),
        kernel_params!(old, atmos_params: released_tuple(temperature_c, humidity_percent)),
        "RK4",
        1e-6,
        0.001,
    )
    .expect("valid kernel input");
    let (time, state) = steps.last().unwrap();
    assert_eq!(*time, seconds);
    endpoint(state)
}

const FLAT: (f64, f64) = (0.002, 2.0);
const LOB: (f64, f64) = (std::f64::consts::FRAC_PI_6, 30.0);

#[cfg(not(target_arch = "wasm32"))]
fn bits((x, y, v): (f64, f64, f64)) -> [u64; 3] {
    [x.to_bits(), y.to_bits(), v.to_bits()]
}

/// The rebuilt "before" is 0.49.0's kernel, bit for bit: the released 0.49.0 crate flies the
/// same shots beside it, on whatever platform runs this. The constants are what the 0.49.0 build
/// itself (main 261c1eb) flew at 15 C / 50 % on macOS, measured when these docs were written;
/// 1e-12 relative leaves room for another libm.
#[test]
fn legacy_kernel_air_reproduces_0_49_0() {
    #[cfg(not(target_arch = "wasm32"))]
    for (temperature_c, humidity_percent) in QUOTED.into_iter().chain([(15.0, 0.0)]) {
        for (shot, (angle_rad, seconds)) in [("flat", FLAT), ("lob", LOB)] {
            assert_eq!(
                bits(fly(
                    legacy_kernel_air(temperature_c, humidity_percent),
                    angle_rad,
                    seconds
                )),
                bits(fly_released(
                    temperature_c,
                    humidity_percent,
                    angle_rad,
                    seconds
                )),
                "{shot} at {temperature_c} C / {humidity_percent} %: the rebuild is not 0.49.0"
            );
        }
    }
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

/// The kernel shift table: each cell is `try_integrate_trajectory` RK4 at 1 ms on 0.50.0 minus the
/// same on 0.49.0, for its own shot, quantity and station.
#[cfg(not(target_arch = "wasm32"))]
#[test]
fn the_0_50_0_trajectory_shifts_are_engine_output() {
    /// `(downrange m, height mm, speed m/s)` of the flat shot and `(downrange m, height m)` of the
    /// lob, 0.50.0 minus 0.49.0.
    type Shift = (f64, f64, f64, f64, f64);
    let shifts: Vec<Shift> = QUOTED
        .iter()
        .map(|&(temperature_c, humidity_percent)| {
            let after = KernelAtmosphere::Station(station(temperature_c, humidity_percent));
            let (x0, y0, v0) = fly_released(temperature_c, humidity_percent, FLAT.0, FLAT.1);
            let (x1, y1, v1) = fly(after, FLAT.0, FLAT.1);
            let (lx0, ly0, _) = fly_released(temperature_c, humidity_percent, LOB.0, LOB.1);
            let (lx1, ly1, _) = fly(after, LOB.0, LOB.1);
            (x1 - x0, (y1 - y0) * 1000.0, v1 - v0, lx1 - lx0, ly1 - ly0)
        })
        .collect();
    let row = |label: String, cell: &dyn Fn(&Shift) -> String| {
        let cells: Vec<String> = shifts.iter().map(cell).collect();
        format!("| {label} | {} |", cells.join(" | "))
    };
    let flat_deg = FLAT.0.to_degrees();
    let expected = [
        row(
            format!("2 s flat ({flat_deg:.2}°): downrange at 2 s"),
            &|s| format!("{:+.3} m", s.0),
        ),
        row("2 s flat: height at 2 s".into(), &|s| {
            format!("{:+.1} mm", s.1)
        }),
        row("2 s flat: speed at 2 s".into(), &|s| {
            format!("{:+.3} m/s", s.2)
        }),
        row("30° lob: downrange at 30 s".into(), &|s| {
            format!("{:+.3} m", s.3)
        }),
        row("30° lob: height at 30 s".into(), &|s| {
            format!("{:+.3} m", s.4)
        }),
    ];
    let header = format!(
        "| Shot | {} |",
        QUOTED
            .iter()
            .map(|(t, rh)| format!("{t} °C / {rh} %"))
            .collect::<Vec<_>>()
            .join(" | ")
    );
    assert_table(&release_section(), &header, &expected);

    // "`fast_integrate` moves the same on the flat shot": its own 2 s flat shot, 0.50.0 minus
    // 0.49.0, at the two conditions the sentence names, equal to the table's at its precision.
    let fast_shift = |temperature_c: f64, humidity_percent: f64| {
        let mut inputs = BallisticInputs {
            bc_value: BC_G7,
            bc_type: DragModel::G7,
            bullet_mass: MASS_KG,
            muzzle_velocity: MUZZLE_VELOCITY_MPS,
            bullet_diameter: DIAMETER_M,
            bullet_length: LENGTH_M,
            twist_rate: 10.0,
            caliber_inches: DIAMETER_M / 0.0254,
            weight_grains: MASS_KG / 0.000_064_798_91,
            temperature: temperature_c,
            pressure: SEA_LEVEL_HPA,
            ground_threshold: -1.0e5,
            ..BallisticInputs::default()
        };
        inputs.humidity_percent = humidity_percent;
        let after = ballistics_engine::fast_trajectory::fast_integrate(
            &inputs,
            &ballistics_engine::wind::WindSock::new(vec![]),
            ballistics_engine::fast_trajectory::FastIntegrationParams {
                horiz: 1.0e6,
                vert: 0.0,
                initial_state: launch(FLAT.0),
                t_span: (0.0, FLAT.1),
                atmosphere: KernelAtmosphere::Station(station(temperature_c, humidity_percent)),
                atmo_sock: None,
            },
        );
        let old_inputs = old::BallisticInputs {
            bc_value: BC_G7,
            bc_type: old::DragModel::G7,
            bullet_mass: MASS_KG,
            muzzle_velocity: MUZZLE_VELOCITY_MPS,
            bullet_diameter: DIAMETER_M,
            bullet_length: LENGTH_M,
            twist_rate: 10.0,
            caliber_inches: DIAMETER_M / 0.0254,
            weight_grains: MASS_KG / 0.000_064_798_91,
            temperature: temperature_c,
            pressure: SEA_LEVEL_HPA,
            ground_threshold: -1.0e5,
            humidity: humidity_percent / 100.0,
            ..old::BallisticInputs::default()
        };
        let before = old::fast_trajectory::fast_integrate(
            &old_inputs,
            &old::wind::WindSock::new(vec![]),
            old::fast_trajectory::FastIntegrationParams {
                horiz: 1.0e6,
                vert: 0.0,
                initial_state: launch(FLAT.0),
                t_span: (0.0, FLAT.1),
                atmo_params: released_tuple(temperature_c, humidity_percent),
                atmo_sock: None,
            },
        );
        let last = |t: &[f64]| t.len() - 1;
        assert_eq!(after.t[last(&after.t)], FLAT.1);
        assert_eq!(before.t[last(&before.t)], FLAT.1);
        after.y[0][last(&after.t)] - before.y[0][last(&before.t)]
    };
    let (f15, f35) = (fast_shift(15.0, 50.0), fast_shift(35.0, 90.0));
    assert_eq!(format!("{f15:+.3}"), format!("{:+.3}", shifts[0].0));
    assert_eq!(format!("{f35:+.3}"), format!("{:+.3}", shifts[2].0));
    assert_prose(
        &release_section(),
        &format!(
            "`fast_integrate` moves the same on the flat shot: {f15:+.3} m at 15 °C / 50 % and \
             {f35:+.3} m at 35 °C / 90 %."
        ),
        "the fast_integrate shift",
    );
}

/// The lean Monte Carlo lines: `solve_trajectory_for_monte_carlo` at its default humidity (50 %),
/// the same bullet at 15 C and sea level, on 0.49.0 and on 0.50.0.
#[cfg(not(target_arch = "wasm32"))]
#[test]
fn the_0_50_0_lean_monte_carlo_shifts_are_engine_output() {
    macro_rules! lean_mc {
        ($engine:ident, $target:expr) => {{
            let inputs = $engine::BallisticInputs {
                bc_value: BC_G7,
                bc_type: $engine::DragModel::G7,
                bullet_mass: MASS_KG,
                muzzle_velocity: MUZZLE_VELOCITY_MPS,
                bullet_diameter: DIAMETER_M,
                bullet_length: LENGTH_M,
                twist_rate: 10.0,
                caliber_inches: DIAMETER_M / 0.0254,
                weight_grains: MASS_KG / 0.000_064_798_91,
                sight_height: 0.0,
                muzzle_height: 0.0,
                altitude: 0.0,
                temperature: 15.0,
                pressure: SEA_LEVEL_HPA,
                ground_threshold: -1.0e5,
                target_distance: $target,
                ..$engine::BallisticInputs::default()
            };
            let out = $engine::monte_carlo::solve_trajectory_for_monte_carlo(&inputs)
                .expect("lean Monte Carlo");
            (out.drop, out.velocity, out.mach)
        }};
    }
    assert_eq!(BallisticInputs::default().humidity_percent, 50.0);
    assert_eq!(old::BallisticInputs::default().humidity, 0.5);
    let section = release_section();
    let (d0, v0, m0) = lean_mc!(old, 500.0);
    let (d1, v1, m1) = lean_mc!(ballistics_engine, 500.0);
    assert_prose(
        &section,
        &format!(
            "- at 500 m, drop {:+.3} mm ({d0:.5} to {d1:.5} m), speed {:+.3} m/s and Mach \
             {m0:.5} to {m1:.5};",
            (d1 - d0) * 1000.0,
            v1 - v0
        ),
        "the lean Monte Carlo at 500 m",
    );
    let (d0, v0, m0) = lean_mc!(old, 1000.0);
    let (d1, v1, m1) = lean_mc!(ballistics_engine, 1000.0);
    assert_prose(
        &section,
        &format!(
            "- at 1000 m, drop {:+.2} mm ({d0:.4} to {d1:.4} m), speed {:+.3} m/s and Mach \
             {m0:.5} to {m1:.5}.",
            (d1 - d0) * 1000.0,
            v1 - v0
        ),
        "the lean Monte Carlo at 1000 m",
    );
}

/// The figures no engine at hand can recompute are labelled as measurements: the aloft change's
/// own share of the kernel and lean Monte Carlo shifts (measured against the development build
/// just before it, which no released crate is), and the CLI golden corpus count.
#[test]
fn the_0_50_0_measured_figures_say_so() {
    let section = release_section();
    assert_prose(
        &section,
        "On the 30° lob that change alone was measured at -0.385 m of the downrange at \
         15 °C / 50 %, -1.020 m at 30 °C / 50 % and -2.515 m at 35 °C / 90 %, against the \
         development build just before it; no released engine reproduces that build, so these \
         are not recomputed.",
        "the kernels' aloft-alone share",
    );
    assert_prose(
        &section,
        "on its own it was measured, the same way, at -0.13 µm of drop and +0.0001 m/s of \
         speed at 500 m, and -0.010 mm and -0.0007 m/s at 1000 m.",
        "the lean Monte Carlo's aloft-alone share",
    );
    assert_prose(
        &section,
        "A golden diff of the CLI corpus against the 0.49.0 binary, run once, measured this: \
         70 of its 418 runs move,",
        "the corpus count",
    );
}

// =================================================================================================
// One vapour model aloft, on every surface
// =================================================================================================

const SEA_LEVEL_STATION_HPA: f64 = 1013.25;

/// The air over a sea-level station at 1013.25 hPa before the aloft change, as the 0.49.0 build
/// (main 261c1eb) and the 4281ccf build flew it, bit for bit alike: `(T C, RH %, height m, density
/// kg/m3, sound speed m/s)`. The rows of the aloft table, in its order.
const AIR_BEFORE_0_50_0: [(f64, f64, f64, f64, f64); 4] = [
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
/// temperature and pressure. They are also what the released 0.49.0 crate gives, run beside
/// this one.
#[test]
fn the_air_before_0_50_0_is_the_old_aloft_model() {
    for (temperature_c, humidity_percent, height_m, density, sound_speed) in AIR_BEFORE_0_50_0 {
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
        #[cfg(not(target_arch = "wasm32"))]
        {
            let (released_density, released_sound_speed) =
                released_air_aloft(temperature_c, humidity_percent, height_m);
            assert!(
                close(released_density, density),
                "0.49.0 density {what}: {released_density:?}"
            );
            assert!(
                close(released_sound_speed, sound_speed),
                "0.49.0 c {what}: {released_sound_speed:?}"
            );
        }
    }
}

/// 0.49.0's per-step air `height_m` over a sea-level station: `TrajectorySolver`'s lookup, the
/// station's CIPM ratio and relative humidity carried up by `get_local_atmosphere_humid`.
#[cfg(not(target_arch = "wasm32"))]
fn released_air_aloft(temperature_c: f64, humidity_percent: f64, height_m: f64) -> (f64, f64) {
    let (.., ratio) = released_tuple(temperature_c, humidity_percent);
    old::atmosphere::get_local_atmosphere_humid(
        height_m,
        0.0,
        temperature_c,
        SEA_LEVEL_STATION_HPA,
        ratio,
        humidity_percent,
    )
}

/// The aloft table, row by row: 0.49.0's air (the constants, checked against the released crate
/// above) against 0.50.0's.
#[test]
fn the_aloft_air_table_is_engine_output() {
    let expected: Vec<String> = AIR_BEFORE_0_50_0
        .iter()
        .map(
            |&(temperature_c, humidity_percent, height_m, density, sound_speed)| {
                let air = air_over_sea_level(temperature_c, humidity_percent, height_m);
                format!(
                    "| {temperature_c} °C / {humidity_percent} % | +{} km | {sound_speed:.3} → \
                     {:.3} m/s ({:+.3}%) | {density:.4} → {:.4} kg/m³ ({:+.3}%) |",
                    height_m / 1000.0,
                    air.speed_of_sound_mps,
                    (air.speed_of_sound_mps / sound_speed - 1.0) * 100.0,
                    air.density_kg_m3,
                    (air.density_kg_m3 / density - 1.0) * 100.0
                )
            },
        )
        .collect();
    assert_table(
        &release_section(),
        "| Station | Height | Sound speed | Density |",
        &expected,
    );
}

#[cfg(not(target_arch = "wasm32"))]
const SIGHT_HEIGHT_M: f64 = 0.0381; // 1.5 in

/// `(drop at 1000 m, speed at 1000 m, lob range, lob apex, lob impact speed)` from
/// `TrajectorySolver` (RK4) on `engine` (this crate or 0.49.0): the flat shot zeroed at 100 m, the
/// lob at 30 degrees until it falls 100 m below the muzzle.
#[cfg(not(target_arch = "wasm32"))]
macro_rules! solver_shots {
    ($engine:ident, $altitude_m:expr, $pressure_hpa:expr, $temperature_c:expr, $humidity:expr) => {{
        use $engine::{
            calculate_zero_angle_with_station_presence, AtmosphericConditions, BallisticInputs,
            DragModel, StationAtmospherePresence, TrajectorySolver, WindConditions,
        };
        let (altitude_m, pressure_hpa, temperature_c) =
            ($altitude_m, $pressure_hpa, $temperature_c);
        let shot_inputs = |angle_rad: f64| BallisticInputs {
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
        };
        let air = AtmosphericConditions {
            temperature: temperature_c,
            pressure: pressure_hpa,
            humidity: $humidity,
            altitude: altitude_m,
        };
        let presence = StationAtmospherePresence::BOTH;
        let zero = calculate_zero_angle_with_station_presence(
            shot_inputs(0.0),
            100.0,
            SIGHT_HEIGHT_M,
            WindConditions::default(),
            air.clone(),
            presence,
        )
        .unwrap();
        let mut flat = shot_inputs(zero);
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
            shot_inputs(std::f64::consts::FRAC_PI_6),
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
    }};
}

/// The shot table, row by row: each cell is `TrajectorySolver` (RK4) on 0.50.0 minus the same on
/// the released 0.49.0, for its own station, shot, quantity and humidity.
#[cfg(not(target_arch = "wasm32"))]
#[test]
fn the_aloft_shot_table_is_engine_output() {
    type Shots = (f64, f64, f64, f64, f64);
    let shots = |altitude_m: f64, pressure_hpa: f64, (t, rh): (f64, f64)| -> (Shots, Shots) {
        (
            solver_shots!(old, altitude_m, pressure_hpa, t, rh),
            solver_shots!(ballistics_engine, altitude_m, pressure_hpa, t, rh),
        )
    };
    let conditions = [(15.0, 50.0), (35.0, 90.0)];
    let sea: Vec<_> = conditions
        .iter()
        .map(|&c| shots(0.0, SEA_LEVEL_STATION_HPA, c))
        .collect();
    let high: Vec<_> = conditions
        .iter()
        .map(|&c| shots(1500.0, 845.6, c))
        .collect();
    let row = |label: &str, runs: &[(Shots, Shots)], cell: &dyn Fn(&Shots, &Shots) -> String| {
        let cells: Vec<String> = runs
            .iter()
            .map(|(before, after)| cell(before, after))
            .collect();
        format!("| {label} | {} |", cells.join(" | "))
    };
    let drop = |b: &Shots, a: &Shots| format!("{:+.3} mm", (a.0 - b.0) * 1000.0);
    let speed = |b: &Shots, a: &Shots| format!("{:+.4} m/s", a.1 - b.1);
    let range = |b: &Shots, a: &Shots| format!("{:+.3} m ({:.2} → {:.2} m)", a.2 - b.2, b.2, a.2);
    let apex = |b: &Shots, a: &Shots| format!("{:+.3} m", a.3 - b.3);
    let impact = |b: &Shots, a: &Shots| format!("{:+.3} m/s", a.4 - b.4);
    let expected = [
        row("Sea level, flat: drop at 1000 m", &sea, &drop),
        row("Sea level, flat: speed at 1000 m", &sea, &speed),
        row("Sea level, 30° lob: range", &sea, &range),
        row("Sea level, 30° lob: apex", &sea, &apex),
        row("Sea level, 30° lob: impact speed", &sea, &impact),
        row(
            "1500 m station (845.6 hPa), flat: drop at 1000 m",
            &high,
            &drop,
        ),
        row("1500 m station, flat: speed at 1000 m", &high, &speed),
        row("1500 m station, 30° lob: range", &high, &range),
    ];
    assert_table(
        &release_section(),
        "| Shot | 15 °C / 50 % | 35 °C / 90 % |",
        &expected,
    );
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

/// The Miller-factor list each document gives, as written: the clause after "multiply Sg by" up to
/// the sentence's end, split into `(factor, "T °C / RH %")` items.
fn miller_factors(text: &str, what: &str) -> Vec<(String, String)> {
    let text = normalized(text);
    let start = text
        .find("multiply Sg by ")
        .unwrap_or_else(|| panic!("{what}: no Miller factor list"))
        + "multiply Sg by ".len();
    let list = &text[start..start + text[start..].find("%.").expect("the list ends") + 1];
    list.split(", ")
        .flat_map(|item| item.split(" and "))
        .map(|item| {
            let (factor, conditions) = item
                .split_once(" at ")
                .unwrap_or_else(|| panic!("{what}: `{item}` is not `F at T °C / RH %`"));
            (factor.to_string(), conditions.to_string())
        })
        .collect()
}

#[test]
fn the_miller_exception_is_documented_with_engine_numbers() {
    for (file, text) in [
        ("src/stability.rs", read("src/stability.rs")),
        ("CHANGELOG 0.50.0", release_section()),
        ("CLI_USAGE.md", read("CLI_USAGE.md")),
    ] {
        let factors = miller_factors(&text, file);
        for (temperature_c, humidity_percent) in QUOTED {
            let conditions = format!("{temperature_c} °C / {humidity_percent} %");
            assert!(
                factors.iter().any(|(_, at)| *at == conditions),
                "{file}: no Miller factor at {conditions}"
            );
        }
        for (factor, conditions) in factors {
            let (temperature, humidity) = conditions
                .strip_suffix(" %")
                .and_then(|c| c.split_once(" °C / "))
                .unwrap_or_else(|| panic!("{file}: `{conditions}`"));
            let engine = format!(
                "{:.4}",
                humid_miller_factor(temperature.parse().unwrap(), humidity.parse().unwrap())
            );
            assert_eq!(
                factor, engine,
                "{file}: the Miller factor at {conditions} is the engine's {engine}"
            );
        }
    }
}

/// What the dry Miller Sg keeps: Litz spin drift at 1000 yd and the aerodynamic jump in a 10 mph
/// crosswind, with the published (dry) Sg and with the Sg a humid correction would give, each in
/// the words each document quotes them with.
#[test]
fn the_miller_exception_quotes_engine_drift_and_jump() {
    use ballistics_engine::aerodynamic_jump::litz_crosswind_jump_moa;
    use ballistics_engine::spin_drift::{effective_sg_from_inputs, litz_drift_inches};
    use ballistics_engine::{
        AtmosphericConditions, StationAtmospherePresence, TrajectorySolver, WindConditions,
    };
    let at = |temperature_c: f64, humidity_percent: f64| {
        let inputs = BallisticInputs {
            bc_value: BC_G7,
            bc_type: DragModel::G7,
            bullet_mass: MASS_KG,
            muzzle_velocity: MUZZLE_VELOCITY_MPS,
            bullet_diameter: DIAMETER_M,
            bullet_length: LENGTH_M,
            twist_rate: 10.0,
            caliber_inches: DIAMETER_M / 0.0254,
            weight_grains: MASS_KG / 0.000_064_798_91,
            temperature: temperature_c,
            pressure: SEA_LEVEL_HPA,
            ground_threshold: -1.0e5,
            ..BallisticInputs::default()
        };
        let sg = effective_sg_from_inputs(&inputs, temperature_c, SEA_LEVEL_HPA);
        let humid_sg = sg * humid_miller_factor(temperature_c, humidity_percent);
        let mut solver = TrajectorySolver::new_with_station_presence(
            inputs,
            WindConditions::default(),
            AtmosphericConditions {
                temperature: temperature_c,
                pressure: SEA_LEVEL_HPA,
                humidity: humidity_percent,
                altitude: 0.0,
            },
            StationAtmospherePresence::BOTH,
        );
        solver.set_max_range(914.4); // 1000 yd
        let time_of_flight = solver.solve().unwrap().time_of_flight;
        let calibers = LENGTH_M / DIAMETER_M;
        (
            litz_drift_inches(sg, time_of_flight),
            litz_drift_inches(humid_sg, time_of_flight),
            litz_crosswind_jump_moa(sg, calibers, 10.0, true),
            litz_crosswind_jump_moa(humid_sg, calibers, 10.0, true),
        )
    };
    let (d15, h15, j15, k15) = at(15.0, 50.0);
    let (d35, h35, j35, k35) = at(35.0, 90.0);
    assert_prose(
        &release_section(),
        &format!(
            "For the bullet above, spin drift at 1000 yd would go from {d15:.2} to {h15:.2} in at \
             15 °C / 50 %, and from {d35:.2} to {h35:.2} in at 35 °C / 90 %. - The aerodynamic \
             jump in a 10 mph crosswind would go from {j15:.4} to {k15:.4} MOA, and from \
             {j35:.4} to {k35:.4} MOA."
        ),
        "CHANGELOG: the Miller exception's drift and jump",
    );
    assert_prose(
        &read("src/stability.rs"),
        &format!(
            "drift at 1000 yd from {d15:.2} to {h15:.2} in at 15 °C / 50 % and from {d35:.2} to \
             {h35:.2} in at 35 °C / 90 %, and the aerodynamic jump in a 10 mph crosswind from \
             {j15:.4} to {k15:.4} MOA and from {j35:.4} to {k35:.4} MOA."
        ),
        "src/stability.rs: the Miller exception's drift and jump",
    );
}

/// Two 0.49.0 behaviours the section quotes numbers for, recomputed on the released crate: the
/// deprecated dry `get_local_atmosphere`'s sound-speed constant, and the steam a NaN humidity
/// flew as.
#[cfg(not(target_arch = "wasm32"))]
#[test]
fn the_0_50_0_section_quotes_0_49_0_behaviour_from_the_released_crate() {
    let section = release_section();
    #[allow(deprecated)]
    let (_, now) =
        ballistics_engine::atmosphere::get_local_atmosphere(0.0, 0.0, 15.0, SEA_LEVEL_HPA, 1.0);
    let (_, then) = old::atmosphere::get_local_atmosphere(0.0, 0.0, 15.0, SEA_LEVEL_HPA, 1.0);
    assert_prose(
        &section,
        &format!(
            "That is {:+.2e} relative, {:+.5} m/s at 15 °C.",
            now / then - 1.0,
            now - then
        ),
        "the deprecated get_local_atmosphere's sound-speed shift",
    );
    let (steam_density, steam_sound_speed) =
        old::atmosphere::calculate_atmosphere(0.0, Some(15.0), Some(SEA_LEVEL_HPA), f64::NAN);
    assert_prose(
        &section,
        &format!(
            "0.49.0 capped the NaN water-vapour mole fraction at 1 and flew pure steam: \
             {steam_density:.3} kg/m³ and {steam_sound_speed:.1} m/s at 15 °C and 1013.25 hPa."
        ),
        "0.49.0's NaN-humidity steam",
    );
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
// Model descriptions and stale wording
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
        // The aloft model before 0.50.0, which re-derived the vapour from the RH aloft. The base
        // air's water-vapour mole fraction is held at every height for the lapse, the density
        // and the sound speed alike.
        (
            "src/derivatives.rs",
            "the sound speed is the moist value at the base RH",
            "the base air's vapour mole fraction is held aloft, not its RH",
        ),
        (
            "src/derivatives.rs",
            "RH, for density AND sound speed) before the altitude lapse",
            "the zone's vapour mole fraction sets the lapse too",
        ),
        (
            "src/trajectory_integration.rs",
            "RH, for density and sound speed) before the altitude lapse",
            "the zone's vapour mole fraction sets the lapse too",
        ),
        (
            "src/fast_trajectory.rs",
            "RH, for\n    /// density and sound speed) as the base before the altitude lapse",
            "the zone's vapour mole fraction sets the lapse too",
        ),
        (
            "src/fast_trajectory.rs",
            "the base humidity for density and sound speed alike",
            "the base air's vapour mole fraction sets the lapse too",
        ),
        (
            "src/monte_carlo.rs",
            "and its RH the flight's sound speed",
            "the station's vapour mole fraction sets the lapse and the sound speed",
        ),
        (
            "src/atmosphere.rs",
            "the relative humidity it gave the sound speed",
            "the trace records the base humidity whose vapour mole fraction is held aloft",
        ),
        (
            "src/atmosphere.rs",
            "so the\n/// pressure aloft falls as the density the projectile flies says it must",
            "the lapse's gas constant and the density's agree only to about 1e-3",
        ),
    ] {
        assert!(!read(path).contains(stale), "{path}: stale comment ({why})");
    }
}

// =================================================================================================
// Disclosures checked against the engine: re-baseline recipes, a kernel caller's height
// reference, the breaking list in both documents, the new items, and the model's own limits
// =================================================================================================

/// [`fly`] with downrange `AtmoSock` zones `(T C, P hPa, RH %, until m)`.
#[cfg(not(target_arch = "wasm32"))]
fn fly_zoned(
    atmosphere: KernelAtmosphere,
    zones: Option<Vec<(f64, f64, f64, f64)>>,
    angle_rad: f64,
    seconds: f64,
) -> (f64, f64, f64) {
    let mut params = params(atmosphere);
    params.atmo_sock = zones.map(ballistics_engine::atmosphere::AtmoSock::new);
    let steps = try_integrate_trajectory(
        launch(angle_rad),
        (0.0, seconds),
        params,
        "RK4",
        1e-6,
        0.001,
    )
    .expect("valid kernel input");
    let (time, state) = steps.last().unwrap();
    assert_eq!(*time, seconds);
    endpoint(state)
}

/// [`fly_released`] with the same zones, on the released 0.49.0 kernel.
#[cfg(not(target_arch = "wasm32"))]
fn fly_released_zoned(
    temperature_c: f64,
    humidity_percent: f64,
    zones: Option<Vec<(f64, f64, f64, f64)>>,
    angle_rad: f64,
    seconds: f64,
) -> (f64, f64, f64) {
    let mut params =
        kernel_params!(old, atmos_params: released_tuple(temperature_c, humidity_percent));
    params.atmo_sock = zones.map(old::atmosphere::AtmoSock::new);
    let steps = old::trajectory_integration::try_integrate_trajectory(
        launch(angle_rad),
        (0.0, seconds),
        params,
        "RK4",
        1e-6,
        0.001,
    )
    .expect("valid kernel input");
    let (time, state) = steps.last().unwrap();
    assert_eq!(*time, seconds);
    endpoint(state)
}

/// The 2 s flat shot's downrange on `fast_integrate`, this crate or 0.49.0, at a sea-level
/// station with the given zones.
#[cfg(not(target_arch = "wasm32"))]
macro_rules! fast_flat_downrange {
    (
        $engine:ident,
        $temperature_c:expr,
        $humidity:ident: $humidity_value:expr,
        $atmosphere:ident: $air:expr,
        $zones:expr
    ) => {{
        let mut inputs = $engine::BallisticInputs {
            temperature: $temperature_c,
            bc_value: BC_G7,
            bc_type: $engine::DragModel::G7,
            bullet_mass: MASS_KG,
            muzzle_velocity: MUZZLE_VELOCITY_MPS,
            bullet_diameter: DIAMETER_M,
            bullet_length: LENGTH_M,
            twist_rate: 10.0,
            caliber_inches: DIAMETER_M / 0.0254,
            weight_grains: MASS_KG / 0.000_064_798_91,
            pressure: SEA_LEVEL_HPA,
            ground_threshold: -1.0e5,
            ..$engine::BallisticInputs::default()
        };
        inputs.$humidity = $humidity_value;
        let solution = $engine::fast_trajectory::fast_integrate(
            &inputs,
            &$engine::wind::WindSock::new(vec![]),
            $engine::fast_trajectory::FastIntegrationParams {
                horiz: 1.0e6,
                vert: 0.0,
                initial_state: launch(FLAT.0),
                t_span: (0.0, FLAT.1),
                $atmosphere: $air,
                atmo_sock: $zones.map($engine::atmosphere::AtmoSock::new),
            },
        );
        let last = solution.t.len() - 1;
        assert_eq!(solution.t[last], FLAT.1);
        solution.y[0][last]
    }};
}

/// `KernelAtmosphere::from_legacy_params(.., RelativeHumidity::DRY)` rebuilds 0.49.0's kernel air
/// bit for bit only for a flight with no `AtmoSock` zones or only dry ones. A zone with a
/// humidity above 0 flies its own moist sound speed and lapse, which 0.49.0's kernels never did,
/// so no recipe rebuilds 0.49.0 there; both documents say so, with the engine's own number.
#[cfg(not(target_arch = "wasm32"))]
#[test]
fn the_legacy_recipe_is_exact_only_without_humid_zones() {
    let humid_zone_shift = |temperature_c: f64, humidity_percent: f64, angle: f64, secs: f64| {
        let zone = Some(vec![(
            temperature_c,
            SEA_LEVEL_HPA,
            humidity_percent,
            1.0e6,
        )]);
        let recipe = fly_zoned(
            legacy_kernel_air(temperature_c, humidity_percent),
            zone.clone(),
            angle,
            secs,
        );
        let released = fly_released_zoned(temperature_c, humidity_percent, zone, angle, secs);
        recipe.0 - released.0
    };
    for (temperature_c, humidity_percent) in QUOTED.into_iter().chain([(15.0, 0.0)]) {
        let dry_zone = vec![(temperature_c, SEA_LEVEL_HPA, 0.0, 1.0e6)];
        let humid_zone = vec![(temperature_c, SEA_LEVEL_HPA, 50.0, 1.0e6)];
        for (shot, (angle, secs)) in [("flat", FLAT), ("lob", LOB)] {
            for zones in [None, Some(dry_zone.clone())] {
                assert_eq!(
                    bits(fly_zoned(
                        legacy_kernel_air(temperature_c, humidity_percent),
                        zones.clone(),
                        angle,
                        secs
                    )),
                    bits(fly_released_zoned(
                        temperature_c,
                        humidity_percent,
                        zones.clone(),
                        angle,
                        secs
                    )),
                    "{shot} at {temperature_c} C / {humidity_percent} %, zones {zones:?}: the \
                     recipe is not 0.49.0"
                );
            }
            // A humid zone, over a humid station or a dry one: not 0.49.0.
            assert_ne!(
                bits(fly_zoned(
                    legacy_kernel_air(temperature_c, humidity_percent),
                    Some(humid_zone.clone()),
                    angle,
                    secs
                )),
                bits(fly_released_zoned(
                    temperature_c,
                    humidity_percent,
                    Some(humid_zone.clone()),
                    angle,
                    secs
                )),
                "{shot} at {temperature_c} C / {humidity_percent} % with a 50 % zone"
            );
        }
        // fast_integrate, which the migration notes also name: the same rule.
        if humidity_percent > 0.0 {
            for zones in [None, Some(dry_zone.clone())] {
                let ratio_tuple = released_tuple(temperature_c, humidity_percent);
                let recipe = fast_flat_downrange!(
                    ballistics_engine,
                    temperature_c,
                    humidity_percent: 0.0,
                    atmosphere: KernelAtmosphere::from_legacy_params(
                        ratio_tuple,
                        RelativeHumidity::DRY
                    )
                    .unwrap(),
                    zones.clone()
                );
                let released = fast_flat_downrange!(
                    old,
                    temperature_c,
                    humidity: humidity_percent / 100.0,
                    atmo_params: ratio_tuple,
                    zones.clone()
                );
                assert_eq!(
                    recipe.to_bits(),
                    released.to_bits(),
                    "fast_integrate at {temperature_c} C / {humidity_percent} %, zones {zones:?}"
                );
            }
        }
    }
    // A zone at the station's own air carries the station's whole humidity shift: the table's.
    let shift = humid_zone_shift(15.0, 50.0, FLAT.0, FLAT.1);
    let station_shift = {
        let (x0, ..) = fly_released(15.0, 50.0, FLAT.0, FLAT.1);
        let (x1, ..) = fly(
            KernelAtmosphere::Station(station(15.0, 50.0)),
            FLAT.0,
            FLAT.1,
        );
        x1 - x0
    };
    assert_eq!(format!("{shift:+.3}"), format!("{station_shift:+.3}"));

    let section = release_section();
    for prose in [
        "With the station and every `AtmoSock` zone at 0 %, every kernel is bit-identical to \
         0.49.0."
            .to_string(),
        "`KernelAtmosphere::from_legacy_params((alt, t, p, ratio), RelativeHumidity::DRY)` \
         rebuilds 0.49.0's kernel air, the humid density with the 0 % sound speed, for a flight \
         with no `AtmoSock` zones or only dry ones; that was checked bit for bit against the \
         0.49.0 build."
            .to_string(),
        format!(
            "A zone with a humidity above 0 flies its own moist sound speed and lapse, at its own \
             height too, which 0.49.0's kernels never did, so no API rebuilds 0.49.0 there: under \
             this recipe a zone at the station's own 15 °C / 50 % still moves the 2 s flat shot \
             by {shift:+.3} m, the shift in the table above."
        ),
        "A zoned caller migrates with `StationAir::new` at the station's humidity, so the station \
         and its zones fly one model."
            .to_string(),
        "(inside an `AtmoSock` zone the zone's own, at the zone's height too)".to_string(),
    ] {
        assert_prose(
            &section,
            &prose,
            "the CHANGELOG's kernel re-baseline recipe",
        );
    }
    let migration = migration_section();
    for prose in [
        "It was checked bit for bit against the 0.49.0 build on `try_integrate_trajectory` and \
         `fast_integrate`, for a flight with no `AtmoSock` zones or only dry ones.",
        "A zone with a humidity above 0 now flies its own moist sound speed and lapse, at its own \
         height too, and no API reproduces 0.49.0 there; migrate a zoned caller with \
         `StationAir::new` at the station's humidity, so the station and its zones fly one model.",
        "with the station and every zone at 0 % nothing moves, bit for bit",
    ] {
        assert_prose(&migration, prose, "MCCOY_MIGRATION.md's re-baseline recipe");
    }
    assert!(
        !normalized(&section).contains("At 0 % every kernel is bit-identical to 0.49.0."),
        "the CHANGELOG still says a kernel at 0 % is 0.49.0 whatever its zones"
    );
}

/// `compute_derivatives` measures a point's height from its `KernelAtmosphere` station, where
/// 0.49.0 measured it from `inputs.altitude`. A caller that passes an absolute altitude in
/// `pos.y` with `inputs.altitude` at 0 (rockets-engine's shape) subtracts the station's altitude
/// and gets 0.49.0's acceleration to the bit; without the shift it reads the air one station
/// altitude too high. The rockets follow-up, the migration notes and the rustdoc say so.
#[cfg(not(target_arch = "wasm32"))]
#[test]
fn an_absolute_height_kernel_caller_subtracts_the_station_altitude() {
    use nalgebra::Vector3;
    let pad_m = 1500.0;
    let (temperature_c, pressure_hpa) = (5.25, 845.6);
    let mut too_little_drag_at_the_pad = f64::NAN;
    for humidity_percent in [0.0, 30.0] {
        let (density, _) = old::atmosphere::calculate_atmosphere(
            pad_m,
            Some(temperature_c),
            Some(pressure_hpa),
            humidity_percent,
        );
        let tuple = (pad_m, temperature_c, pressure_hpa, density / 1.225);
        let kernel = KernelAtmosphere::from_legacy_params(tuple, RelativeHumidity::DRY).unwrap();
        let old_inputs = old::BallisticInputs {
            altitude: 0.0,
            ..Default::default()
        };
        let inputs = BallisticInputs {
            altitude: 0.0,
            ..Default::default()
        };
        let velocity = Vector3::new(150.0, 90.0, 3.0);
        for absolute_m in [1500.0, 1600.1, 2500.0, 5000.0] {
            let then = old::derivatives::compute_derivatives(
                Vector3::new(100.0, absolute_m, 0.0),
                velocity,
                &old_inputs,
                Vector3::zeros(),
                tuple,
                0.3,
                None,
                0.0,
                None,
            );
            let at = |height_m: f64| {
                ballistics_engine::derivatives::compute_derivatives(
                    Vector3::new(100.0, height_m, 0.0),
                    velocity,
                    &inputs,
                    Vector3::zeros(),
                    &kernel,
                    0.3,
                    None,
                    0.0,
                    None,
                )
            };
            let shifted = at(absolute_m - pad_m);
            assert_eq!(
                shifted.map(f64::to_bits),
                then.map(f64::to_bits),
                "{humidity_percent} % at {absolute_m} m: the shifted height is not 0.49.0"
            );
            let unshifted = at(absolute_m);
            assert_ne!(
                unshifted[3], then[3],
                "the unshifted height reads the same air"
            );
            if humidity_percent == 0.0 && absolute_m == pad_m {
                too_little_drag_at_the_pad = (1.0 - unshifted[3] / then[3]) * 100.0;
            }
        }
    }
    let section = release_section();
    assert_prose(
        &section,
        &format!(
            "`compute_derivatives` now measures a point's height from the `KernelAtmosphere` \
             station, so `pos.y` is the height above the station's altitude. rockets-engine \
             passes the absolute altitude in `pos.y` (with `inputs.altitude` at 0): subtract \
             `StationAir::altitude_m()` from it, and a dry flight's acceleration is 0.49.0's to \
             the bit. Without the shift the air is read one station altitude too high: at a \
             1500 m pad, {too_little_drag_at_the_pad:.0} % too little drag."
        ),
        "the rockets-engine follow-up",
    );
    assert_prose(
        &section,
        "`air_at`, which `Air::at` moves to, still takes the absolute altitude.",
        "the rockets-engine follow-up's second height convention",
    );
    assert!(
        !section.contains("Set `inputs.humidity_percent`"),
        "compute_derivatives does not read the inputs' humidity, so setting it is a no-op"
    );
    assert_prose(
        &migration_section(),
        "A caller that passes an absolute altitude in `pos.y` (with `inputs.altitude` at 0) \
         subtracts the station's altitude, `StationAir::altitude_m()`, from it; then a dry \
         flight's acceleration is 0.49.0's to the bit.",
        "MCCOY_MIGRATION.md's height reference",
    );
    assert_prose(
        &read("src/derivatives.rs"),
        "A caller holding an absolute altitude passes its height above the station, `altitude - \
         station.altitude_m()`, as `pos.y`.",
        "compute_derivatives' rustdoc",
    );
}

/// Every Rust-source break and deprecation is in both the CHANGELOG's breaking list and the
/// migration notes, so a binding maintainer who reads either one meets each of them.
#[test]
fn every_breaking_change_is_in_the_changelog_and_the_migration_notes() {
    let (section, migration) = (release_section(), migration_section());
    for name in [
        "atmos_params",
        "atmo_params",
        "KernelAtmosphere",
        "compute_derivatives",
        "try_solve_trajectory_rust",
        "solve_trajectory_rust",
        "aerodynamic_jump_launch_offset_rad",
        "try_integrate_trajectory",
        "KernelInputError",
        "humidity_percent",
        "humidity_percent()",
        "get_local_atmosphere",
        "from_legacy_params",
        "ToleranceAxisV1",
        "near_at_physical_boundary",
        "far_at_physical_boundary",
        "compute_wez",
        "compute_wez_with_humidity",
        "WindScenarioError",
        "InvalidHumidity",
        "ZeroFailed",
        "set_humidity",
        "InputAxis",
        "inputs.altitude",
    ] {
        assert_contains(&section, name, "the CHANGELOG's breaking list");
        assert_contains(&migration, name, "MCCOY_MIGRATION.md's 0.50.0 notes");
    }
    for deprecated in [
        "`BallisticInputs::humidity_percent()` accessor is `#[deprecated]`",
        "`wez::compute_wez`, the 24-argument form, is `#[deprecated]`",
    ] {
        assert_prose(&section, deprecated, "the CHANGELOG's deprecation shims");
    }
}

/// Every new public item outside a new type is named in the CHANGELOG, and the one deprecation
/// note on the humidity work says how to migrate, constructor included.
#[test]
fn the_0_50_0_section_names_every_new_public_item() {
    let section = release_section();
    for name in [
        "ProfileData::validate_humidity",
        "atmosphere::reading_text",
        "flown_humidity_note",
        "stored_station",
        "air_refusal_in",
        "a7p_air_refusal",
        "physical_domain",
        "StationPressureOrigin",
        "AtmoSock::validate",
        "compute_wez_with_humidity",
    ] {
        assert_contains(&section, name, "the CHANGELOG's new items");
    }
    assert_contains(
        &read("src/atmosphere.rs"),
        "note = \"dry air by contract; use air_at(&KernelAtmosphere::Station(StationAir::new(..)?), \
         None, 0.0, altitude_m), or StationAir::with_density_ratio for a caller-held density \
         ratio\"",
        "get_local_atmosphere's deprecation note",
    );
}

/// The water-vapour mole fraction a first-order moist sound speed implies at `temp_k`: the
/// engine's `c = sqrt(1.4 (1 - 0.062 x_v) · 287.0531 (1 + 0.378 x_v) · T)` solved for `x_v`, in
/// the cancellation-free form of the root near 0.
fn x_v_from_sound_speed(sound_speed_mps: f64, temp_k: f64) -> f64 {
    let q = sound_speed_mps * sound_speed_mps / (1.4 * 287.0531 * temp_k) - 1.0;
    let (b, a) = (0.378 - 0.062, 0.062 * 0.378);
    2.0 * q / (b + (b * b - 4.0 * a * q).sqrt())
}

/// Over a sea-level station at 1013.25 hPa, the air `height_m` up holds the station's
/// water-vapour mole fraction: this many times the most the air there could hold, its saturation
/// mole fraction at the local temperature and pressure (the engine's own Wagner-Pruss pressure
/// and enhancement factor, read through its sound speed at 100 %).
fn times_saturation(temperature_c: f64, humidity_percent: f64, height_m: f64) -> f64 {
    let station = station(temperature_c, humidity_percent);
    let held = x_v_from_sound_speed(station.speed_of_sound_mps(), temperature_c + 273.15);
    let air = air_at(&KernelAtmosphere::Station(station), None, 0.0, height_m);
    let held_aloft = x_v_from_sound_speed(air.speed_of_sound_mps, air.temperature_k);
    assert!((held_aloft / held - 1.0).abs() < 1e-9, "x_v is held aloft");
    let saturated = x_v_from_sound_speed(
        moist_speed_of_sound(air.temperature_k, air.pressure_pa, 100.0),
        air.temperature_k,
    );
    held / saturated
}

/// The height above a sea-level station where its held vapour reaches saturation.
fn saturation_height_m(temperature_c: f64, humidity_percent: f64) -> f64 {
    let (mut below, mut above) = (0.0_f64, 12_000.0_f64);
    for _ in 0..80 {
        let mid = 0.5 * (below + above);
        if times_saturation(temperature_c, humidity_percent, mid) < 1.0 {
            below = mid;
        } else {
            above = mid;
        }
    }
    0.5 * (below + above)
}

/// The held vapour's supersaturation aloft, quantified from the engine at every disclosure site:
/// every height above a saturated station is supersaturated, and the crossing heights and the
/// multiples at +3 km and +10 km are the engine's.
#[test]
fn the_supersaturation_aloft_is_quantified_from_the_engine() {
    for temperature_c in [0.0, 15.0, 35.0, 45.0] {
        for height_m in [1.0, 10.0, 100.0, 1000.0, 3000.0, 10_000.0] {
            let times = times_saturation(temperature_c, 100.0, height_m);
            assert!(
                times > 1.0,
                "{temperature_c} C / 100 % at +{height_m} m: {times}"
            );
        }
    }
    let sentence = format!(
        "Held aloft, `x_v` exceeds saturation in the colder air: at every height above a \
         saturated station, and from about {:.1} km above a 15 °C / 50 % station and {:.1} km \
         above a 35 °C / 90 % one. At +3 km it is {:.2} times saturation over the first and {:.2} \
         times over the second, and at +10 km {:.0} and {:.0} times. Condensation is not \
         modelled.",
        saturation_height_m(15.0, 50.0) / 1000.0,
        saturation_height_m(35.0, 90.0) / 1000.0,
        times_saturation(15.0, 50.0, 3000.0),
        times_saturation(35.0, 90.0, 3000.0),
        times_saturation(15.0, 50.0, 10_000.0),
        times_saturation(35.0, 90.0, 10_000.0),
    );
    assert_prose(&release_section(), &sentence, "the CHANGELOG's aloft model");
    assert_prose(
        &read("CLI_USAGE.md"),
        &sentence,
        "CLI_USAGE.md's aloft model",
    );
    // The atmosphere module's doc and air_at's own.
    assert_eq!(
        normalized(&read("src/atmosphere.rs"))
            .matches(&sentence)
            .count(),
        2,
        "src/atmosphere.rs states it in the module doc and on air_at: {sentence}"
    );
}

/// The engine's gas constants for one held `x_v`, checked against their own docs:
/// - the lapse's ideal-gas `R_AIR / (M_m / M_a)` against the `Z R / M_m` the CIPM-2007 density
///   implies (its compressibility, and its `R / M_a` against `R_AIR`), over a grid of stations;
/// - the sound speed's first-order `R_AIR (1 + 0.378 x_v)` against the lapse's exact one.
#[test]
fn the_gas_constant_gaps_are_quantified_from_the_engine() {
    const R_AIR: f64 = 287.0531;
    const M_V_OVER_M_A: f64 = 18.015_28e-3 / 28.965_46e-3;
    let gaps = |temperature_c: f64, pressure_hpa: f64, humidity_percent: f64| {
        let station = StationAir::new(
            0.0,
            temperature_c,
            pressure_hpa,
            RelativeHumidity::try_percent(humidity_percent).unwrap(),
        )
        .unwrap();
        let temp_k = temperature_c + 273.15;
        let x_v = x_v_from_sound_speed(station.speed_of_sound_mps(), temp_k);
        let lapse_r = R_AIR / (1.0 - x_v * (1.0 - M_V_OVER_M_A));
        let density_r = pressure_hpa * 100.0 / (station.density_kg_m3() * temp_k);
        let sound_r = R_AIR * (1.0 + 0.378 * x_v);
        let sound_speed = station.speed_of_sound_mps();
        let with_lapse_r = (1.4 * (1.0 - 0.062 * x_v) * lapse_r * temp_k).sqrt();
        (
            lapse_r / density_r - 1.0,
            sound_r / lapse_r - 1.0,
            with_lapse_r - sound_speed,
            (with_lapse_r / sound_speed - 1.0) * 100.0,
        )
    };
    let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
    for temperature_c in [-20.0, 0.0, 15.0, 30.0, 45.0] {
        for pressure_hpa in [700.0, 850.0, 1013.25] {
            for humidity_percent in [0.0, 50.0, 100.0] {
                let (lapse_vs_density, ..) = gaps(temperature_c, pressure_hpa, humidity_percent);
                low = low.min(lapse_vs_density);
                high = high.max(lapse_vs_density);
            }
        }
    }
    let atmosphere = read("src/atmosphere.rs");
    assert_prose(
        &atmosphere,
        &format!(
            "The two gas constants still differ a little, dry air included: the lapse's \
             ideal-gas `R_AIR / (M_m / M_a)` and the `Z R / M_m` the CIPM-2007 density implies \
             (its compressibility `Z`, and its `R / M_a = 287.0478` against `R_AIR = 287.0531`) \
             differ by {low:.1e} to {high:.1e} from -20 to 45 °C, 700 to 1013.25 hPa and 0 to \
             100 %."
        ),
        "the lapse's gas constant against the density's",
    );
    let (_, mild, ..) = gaps(15.0, 1013.25, 50.0);
    let (_, hot, low_speed, low_percent) = gaps(45.0, 1013.25, 100.0);
    let (.., thin_speed, _) = gaps(45.0, 700.0, 100.0);
    assert_prose(
        &atmosphere,
        &format!(
            "Its gas constant, `287.0531 (1 + 0.378 x_v)`, is first order in `x_v`: it differs \
             from the exact `R_AIR / (M_m / M_a)` the lapse uses by {mild:.1e} at 15 °C / 50 % \
             and {hot:.1e} at 45 °C / 100 % (1013.25 hPa), where the sound speed is \
             {low_speed:.2} m/s ({low_percent:.3} %) below what the lapse's gas constant gives. \
             The gap grows as the pressure falls: {thin_speed:.2} m/s at 45 °C / 100 % and \
             700 hPa."
        ),
        "the sound speed's gas constant against the lapse's",
    );
}

/// The default sentence names its exceptions, in all three documents that state it.
#[test]
fn the_default_is_stated_with_its_exceptions() {
    let exceptions = "unless another reading supplies it (a saved profile, the shot day for \
                      `--zero-humidity`), with two exceptions: an `.a7p` file's absent humidity \
                      is 0 % (the format's own meaning), and the bridge's `true.*` commands \
                      require `humidity_pct`.";
    for (what, text, blanket) in [
        (
            "the CHANGELOG",
            release_section(),
            "An omitted humidity is 50 % everywhere.",
        ),
        (
            "src/atmosphere.rs",
            read("src/atmosphere.rs"),
            "An omitted humidity is [`DEFAULT_HUMIDITY_PERCENT`], 50 %, on every surface.",
        ),
        (
            "CLI_USAGE.md",
            read("CLI_USAGE.md"),
            "(0.50.0 made this one rule everywhere)",
        ),
    ] {
        assert_prose(&text, exceptions, what);
        assert!(
            !normalized(&text).contains(blanket),
            "{what}: the blanket default `{blanket}` contradicts its documented exceptions"
        );
    }
}

/// Every public field that carries a humidity reading says in its rustdoc that a reading within
/// the sensor band flies as the end it reads as, so a binding author neither refuses 101 nor
/// clamps it by hand. A field that is not a reading (an echo, a test trace) is listed with why.
#[test]
fn every_public_humidity_field_doc_states_the_sensor_band() {
    const NOT_A_READING: [(&str, &str, &str); 3] = [
        (
            "src/atmosphere.rs",
            "humidity_sent",
            "the value a refusal message quotes, as sent",
        ),
        (
            "src/atmosphere.rs",
            "humidity_percent",
            "the test-only air trace's record of the base humidity",
        ),
        (
            "src/truing.rs",
            "humidity",
            "the pub(crate) forward model's copy of the validated request humidity",
        ),
    ];
    let mut files = Vec::new();
    let mut pending = vec![Path::new(env!("CARGO_MANIFEST_DIR")).join("src")];
    while let Some(dir) = pending.pop() {
        for entry in fs::read_dir(&dir).unwrap().flatten() {
            let path = entry.path();
            if path.is_dir() {
                pending.push(path);
            } else if path.extension().is_some_and(|ext| ext == "rs")
                && !path.to_string_lossy().ends_with("_tests.rs")
            {
                files.push(path);
            }
        }
    }
    let mut checked = 0;
    for path in files {
        let relative = path
            .strip_prefix(env!("CARGO_MANIFEST_DIR"))
            .unwrap()
            .to_string_lossy()
            .replace('\\', "/");
        let text = fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        for (index, line) in lines.iter().enumerate() {
            let Some(rest) = line.trim_start().strip_prefix("pub ") else {
                continue;
            };
            let Some((name, ty)) = rest.split_once(':') else {
                continue;
            };
            let ty = ty.trim().trim_end_matches(',');
            if !name.contains("humid") || !matches!(ty, "f64" | "Option<f64>" | "c_double") {
                continue;
            }
            if NOT_A_READING
                .iter()
                .any(|&(file, field, _)| file == relative && field == name)
            {
                continue;
            }
            // The field's doc comment: the `///` lines above it, past any attributes (a
            // `#[serde(...)]` may span several lines).
            let mut doc = Vec::new();
            let mut in_attribute = false;
            for above in lines[..index].iter().rev().map(|line| line.trim()) {
                if in_attribute {
                    in_attribute = !above.starts_with("#[");
                } else if above.starts_with("///") {
                    doc.push(above);
                } else if above.starts_with("#[") {
                    in_attribute = !above.ends_with(']');
                } else if above.ends_with(")]") {
                    in_attribute = true;
                } else {
                    break;
                }
            }
            let doc = doc.join(" ");
            assert!(
                doc.contains("HUMIDITY_SENSOR_TOLERANCE_PERCENT")
                    || doc.contains("5 points")
                    || doc.contains("within the band"),
                "{relative}:{}: `pub {name}` does not say how a sensor reading near 0-100 flies",
                index + 1
            );
            checked += 1;
        }
    }
    assert!(checked >= 14, "only {checked} humidity fields were checked");
}

/// The commands that take a saved profile (`--saved-profile <NAME>` or `--profile <NAME>`; not
/// `trajectory`'s `--profile <FILE>`, a CSV) and fly a humidity: those that take `--humidity`,
/// and `dsf`, which takes no station flag and flies the stored air whole.
#[cfg(feature = "cli")]
fn saved_profile_humidity_commands() -> std::collections::BTreeSet<String> {
    let mut found = std::collections::BTreeSet::new();
    let mut pending: Vec<Vec<String>> = vec![Vec::new()];
    while let Some(path) = pending.pop() {
        let help = help_of(&path);
        let takes_profile = help.lines().map(str::trim_start).any(|line| {
            line.starts_with("--saved-profile <NAME>") || line.starts_with("--profile <NAME>")
        });
        if !path.is_empty()
            && takes_profile
            && (declares_option(&help, "humidity") || path == ["dsf"])
        {
            found.insert(path.join(" "));
        }
        for name in subcommands_in(&help) {
            let mut child = path.clone();
            child.push(name);
            pending.push(child);
        }
    }
    found
}

/// The humidity table's saved-profile row names exactly the commands that fly a saved profile's
/// stored humidity (0.50.0: every command that takes a saved profile and flies a humidity; before
/// it, the row named `--saved-profile` and `--profile` while twelve `--profile` commands flew 50 %).
#[cfg(feature = "cli")]
#[test]
fn the_cli_usage_humidity_table_names_every_command_that_flies_a_stored_humidity() {
    let flown = saved_profile_humidity_commands();
    for known in [
        "trajectory",
        "dsf",
        "come-ups",
        "compare",
        "reticle hold",
        "adaptive-card",
    ] {
        assert!(
            flown.contains(known),
            "the help walk missed `{known}`: {flown:?}"
        );
    }
    assert!(
        !flown.contains("stability") && !flown.contains("dial-plan"),
        "stability's Miller correction is dry and dial-plan flies no air: {flown:?}"
    );
    assert_eq!(
        humidity_table_commands("`humidity`"),
        flown,
        "CLI_USAGE.md's humidity table must name exactly the commands that fly a saved \
         profile's stored humidity"
    );
}

/// `come-ups` on one saved profile at 1000 yd, by `ballistics` itself: `(velocity fps, come-up
/// mil)`, with `extra` typed after the profile.
#[cfg(feature = "cli")]
fn come_ups_on_a_saved_profile(home: &Path, extra: &[&str]) -> (f64, f64) {
    let output = std::process::Command::new(env!("CARGO_BIN_EXE_ballistics"))
        .env("HOME", home)
        .args([
            "come-ups",
            "--profile",
            "hot",
            "--zero-distance",
            "100",
            "--start",
            "1000",
            "--end",
            "1000",
            "-o",
            "json",
        ])
        .args(extra)
        .output()
        .expect("spawn ballistics");
    assert!(
        output.status.success(),
        "come-ups failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let doc: serde_json::Value = serde_json::from_slice(&output.stdout).expect("come-ups JSON");
    let row = &doc["data"][0];
    (
        row["velocity"].as_f64().expect("velocity"),
        row["come_up"].as_f64().expect("come_up"),
    )
}

/// The CHANGELOG's measured example of the saved-profile weather change is engine output. The
/// 0.49.0 row is the profile flown with the flags' defaults typed (59 °F, 50 %), the air 0.49.0's
/// `come-ups --profile` flew for every profile; the released 0.49.0 CLI gives 1199.5585 fps and
/// 10.8345 mil there, the same at the table's precision (the aloft change moves the fourth
/// decimal).
#[cfg(feature = "cli")]
#[test]
fn the_saved_profile_weather_example_is_engine_output() {
    let home =
        std::env::temp_dir().join(format!("bx-humidity-docs-profile-{}", std::process::id()));
    let _ = fs::remove_dir_all(&home);
    fs::create_dir_all(&home).unwrap();
    let saved = std::process::Command::new(env!("CARGO_BIN_EXE_ballistics"))
        .env("HOME", &home)
        .args([
            "profile",
            "save",
            "hot",
            "-v",
            "2700",
            "-b",
            "0.475",
            "-m",
            "168",
            "-d",
            "0.308",
            "--zero-distance",
            "100",
            "--temperature",
            "95",
            "--humidity",
            "100",
        ])
        .output()
        .expect("spawn ballistics");
    assert!(saved.status.success(), "profile save failed");
    let (old_v, old_c) =
        come_ups_on_a_saved_profile(&home, &["--temperature", "59", "--humidity", "50"]);
    let (new_v, new_c) = come_ups_on_a_saved_profile(&home, &[]);
    let _ = fs::remove_dir_all(&home);
    // The table sits inside its list item, so its lines are indented.
    let section = release_section();
    let header = "| `come-ups --profile hot`, 1000 yd | velocity | come-up |";
    let rows: Vec<&str> = section
        .lines()
        .map(str::trim_start)
        .skip_while(|line| *line != header)
        .skip(2)
        .take_while(|line| line.starts_with('|'))
        .collect();
    assert_eq!(
        rows,
        [
            format!(
                "| 0.49.0: the flags' defaults, 59 °F and 50 % | {old_v:.2} fps | {old_c:.2} mil |"
            ),
            format!(
                "| 0.50.0: the stored weather, 95 °F and 100 % | {new_v:.2} fps | {new_c:.2} mil |"
            ),
        ],
        "the table headed `{header}` is not the engine's"
    );
}

/// `ballistics` run with `HOME` at `home`, which must succeed; its stdout as JSON.
#[cfg(feature = "cli")]
fn cli_json(home: &Path, args: &[&str]) -> serde_json::Value {
    let output = std::process::Command::new(env!("CARGO_BIN_EXE_ballistics"))
        .env("HOME", home)
        .args(args)
        .output()
        .expect("spawn ballistics");
    assert!(
        output.status.success(),
        "{args:?} failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).expect("JSON output")
}

/// The CHANGELOG's two other measured figures for the saved-profile change are engine output.
/// What 0.49.0 flew is reproduced by typing it: the converted 29.92 inHg a weatherless imperial
/// profile carried into a metric `trajectory` run (at sea level its presence changes nothing),
/// and a metric profile's numbers read as imperial ones by `reticle hold` (the profile is not
/// converted). Both were checked against the build before the change, bit for bit.
#[cfg(feature = "cli")]
#[test]
fn the_saved_profile_weather_figures_are_engine_output() {
    let home = std::env::temp_dir().join(format!("bx-humidity-docs-units-{}", std::process::id()));
    let _ = fs::remove_dir_all(&home);
    fs::create_dir_all(&home).unwrap();
    let ok = |args: &[&str]| {
        let output = std::process::Command::new(env!("CARGO_BIN_EXE_ballistics"))
            .env("HOME", &home)
            .args(args)
            .output()
            .expect("spawn ballistics");
        assert!(output.status.success(), "{args:?} failed");
        output.stdout
    };
    ok(&[
        "profile",
        "save",
        "plain",
        "-v",
        "2700",
        "-b",
        "0.475",
        "-m",
        "168",
        "-d",
        "0.308",
        "--zero-distance",
        "100",
    ]);
    ok(&[
        "-u",
        "metric",
        "profile",
        "save",
        "m823",
        "-v",
        "823",
        "-b",
        "0.475",
        "-m",
        "10.9",
        "-d",
        "7.82",
        "--zero-distance",
        "91.44",
    ]);
    let reticle = home.join("reticle.json");
    fs::write(
        &reticle,
        ok(&[
            "reticle",
            "generate",
            "tree",
            "--rows",
            "8",
            "--row-spacing",
            "1.0",
            "--spread-step",
            "0.5",
            "-o",
            "json",
        ]),
    )
    .unwrap();
    let reticle = reticle.to_str().unwrap();

    let converted_default = format!("{:?}", 29.92_f64 * 33.8639);
    let trajectory = |extra: &[&str]| {
        let mut args = vec![
            "-u",
            "metric",
            "trajectory",
            "--max-range",
            "600",
            "--auto-zero",
            "100",
            "-o",
            "json",
            "--saved-profile",
            "plain",
        ];
        args.extend_from_slice(extra);
        cli_json(&home, &args)["impact_velocity"]
            .as_f64()
            .expect("impact_velocity")
    };
    let before = trajectory(&["--pressure", &converted_default]);
    let after = trajectory(&[]);
    let hold = |extra: &[&str]| {
        let mut args = vec![
            "reticle",
            "hold",
            "--reticle-json",
            reticle,
            "--range",
            "300",
            "-o",
            "json",
        ];
        args.extend_from_slice(extra);
        cli_json(&home, &args)["hold"]["down_mil"]
            .as_f64()
            .expect("hold.down_mil")
    };
    let read_as_imperial = hold(&[
        "-v",
        "823",
        "-b",
        "0.475",
        "-m",
        "10.9",
        "-d",
        "7.82",
        "--zero-distance",
        "91.44",
    ]);
    let converted = hold(&["--profile", "m823"]);
    let _ = fs::remove_dir_all(&home);
    let section = release_section();
    assert_prose(
        &section,
        &format!("moves from {before:.4} to {after:.4} m/s at 600 m"),
        "the weatherless cross-unit trajectory figure",
    );
    assert_prose(
        &section,
        &format!("`reticle hold --range 300` held {read_as_imperial:.2} mil; it now holds {converted:.2} mil"),
        "the unconverted-profile reticle hold figure",
    );
}
