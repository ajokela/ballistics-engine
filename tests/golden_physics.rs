//! MBA-1318 — Golden-reference physics validation harness.
//!
//! A versioned dataset of reference trajectories (`data/validation/*.json`) with per-source
//! tolerance bands, replayed through the real [`TrajectorySolver`]. The point is to catch physics
//! regressions in CI instead of in the field: any change that shifts drop / drift / time-of-flight /
//! velocity beyond the justified tolerance for a case makes this test fail with a rich, aligned
//! delta table naming the exact case and observable.
//!
//! Gated behind the non-default `validation` feature so the everyday `cargo test` build stays fast
//! and dependency-free. Run it with:
//!
//! ```text
//! cargo test --features validation --test golden_physics
//! ```
//!
//! The harness itself pulls in **zero** new dependencies: `serde` / `serde_json` are already core
//! dependencies of the crate. See `data/validation/README.md` for the case-file schema and the
//! provenance rules for each `source.kind`.
#![cfg(feature = "validation")]

use ballistics_engine::drag::DragTable;
use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, DragModel, TrajectoryResult, TrajectorySolver,
    WindConditions,
};
use serde::Deserialize;
use std::fs;
use std::path::{Path, PathBuf};

/// Minimum number of case files we expect to find. A guard so that a data directory that has been
/// emptied / mis-pathed fails loudly instead of silently "passing" with nothing to check.
const MIN_EXPECTED_CASES: usize = 10;

// --------------------------------------------------------------------------------------------
// Case-file schema (mirrors data/validation/README.md). serde defaults keep the JSON terse: any
// field left out takes the documented engine default.
// --------------------------------------------------------------------------------------------

#[derive(Debug, Deserialize)]
struct Case {
    id: String,
    #[allow(dead_code)] // present for humans / provenance; not used by the runner
    description: String,
    source: Source,
    inputs: Inputs,
    expectations: Vec<Expectation>,
    /// A case we keep on purpose even though the engine disagrees with the reference: a REAL,
    /// tracked discrepancy. Its expectations are reported informationally and excluded from the
    /// hard pass/fail so the suite stays green while the discrepancy is worked (see tracking_note).
    #[serde(default)]
    known_discrepancy: bool,
    #[serde(default)]
    tracking_note: Option<String>,
}

#[derive(Debug, Deserialize)]
struct Source {
    /// "analytic" | "cross-implementation" | "published"
    kind: String,
    #[allow(dead_code)]
    citation: String,
    #[allow(dead_code)]
    #[serde(default)]
    retrieved: Option<String>,
    #[allow(dead_code)]
    #[serde(default)]
    generator_version: Option<String>,
}

fn default_true() -> bool {
    true
}

#[derive(Debug, Deserialize)]
struct Inputs {
    muzzle_velocity_mps: f64,
    bullet_mass_kg: f64,
    bullet_diameter_m: f64,
    bullet_length_m: f64,
    bc_value: f64,
    bc_type: String,
    #[serde(default)]
    muzzle_angle_rad: f64,
    #[serde(default)]
    shooting_angle_rad: f64,
    #[serde(default)]
    twist_rate_in: f64,
    #[serde(default = "default_true")]
    is_twist_right: bool,
    #[serde(default)]
    sight_height_m: f64,
    #[serde(default)]
    muzzle_height_m: f64,
    #[serde(default)]
    altitude_m: f64,
    temperature_c: f64,
    pressure_hpa: f64,
    #[serde(default)]
    humidity_frac: f64,
    #[serde(default)]
    wind_speed_mps: f64,
    #[serde(default)]
    wind_angle_rad: f64,
    #[serde(default)]
    vertical_wind_mps: f64,
    #[serde(default = "default_true")]
    use_rk4: bool,
    #[serde(default = "default_true")]
    use_adaptive_rk45: bool,
    #[serde(default)]
    use_enhanced_spin_drift: bool,
    #[serde(default)]
    enable_advanced_effects: bool,
    #[serde(default)]
    enable_magnus: bool,
    #[serde(default)]
    enable_coriolis: bool,
    #[serde(default)]
    latitude_deg: Option<f64>,
    #[serde(default)]
    shot_azimuth_rad: f64,
    /// Harness directive: install an all-zero custom drag table so the only force acting is
    /// gravity. Turns the solve into an exact vacuum trajectory (the analytic anchors rely on it).
    #[serde(default)]
    zero_drag: bool,
    /// Override the ground-stop threshold. Defaults to a very negative value so the projectile
    /// always reaches the furthest observation range rather than being stopped by ground impact.
    #[serde(default)]
    ground_threshold_m: Option<f64>,
}

#[derive(Debug, Deserialize)]
struct Expectation {
    /// "drop_m" | "drift_m" | "tof_s" | "velocity_mps"
    observable: String,
    range_m: f64,
    value: f64,
    #[serde(default)]
    tol_abs: Option<f64>,
    #[serde(default)]
    tol_rel: Option<f64>,
    #[allow(dead_code)]
    tolerance_justification: String,
}

// --------------------------------------------------------------------------------------------
// Runner
// --------------------------------------------------------------------------------------------

fn validation_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("data/validation")
}

/// Build the engine's `BallisticInputs` from a case's `inputs` block. Imperial mirror fields
/// (`caliber_inches` / `weight_grains`) are derived from the SI values so the stability / spin-drift
/// path stays self-consistent. Wind and atmosphere are carried on the constructor args, not here —
/// the solver integrates from `self.wind` / `self.atmosphere`, not these fields.
fn build_inputs(inp: &Inputs) -> BallisticInputs {
    let bc_type = DragModel::from_str(&inp.bc_type)
        .unwrap_or_else(|| panic!("unknown bc_type/drag model: {}", inp.bc_type));

    let mut b = BallisticInputs {
        muzzle_velocity: inp.muzzle_velocity_mps,
        bullet_mass: inp.bullet_mass_kg,
        bullet_diameter: inp.bullet_diameter_m,
        bullet_length: inp.bullet_length_m,
        bc_value: inp.bc_value,
        bc_type,
        muzzle_angle: inp.muzzle_angle_rad,
        shooting_angle: inp.shooting_angle_rad,
        twist_rate: inp.twist_rate_in,
        is_twist_right: inp.is_twist_right,
        sight_height: inp.sight_height_m,
        muzzle_height: inp.muzzle_height_m,
        altitude: inp.altitude_m,
        temperature: inp.temperature_c,
        pressure: inp.pressure_hpa,
        humidity: inp.humidity_frac,
        wind_speed: inp.wind_speed_mps,
        wind_angle: inp.wind_angle_rad,
        use_rk4: inp.use_rk4,
        use_adaptive_rk45: inp.use_adaptive_rk45,
        use_enhanced_spin_drift: inp.use_enhanced_spin_drift,
        enable_advanced_effects: inp.enable_advanced_effects,
        enable_magnus: inp.enable_magnus,
        enable_coriolis: inp.enable_coriolis,
        latitude: inp.latitude_deg,
        shot_azimuth: inp.shot_azimuth_rad,
        caliber_inches: inp.bullet_diameter_m / 0.0254,
        weight_grains: inp.bullet_mass_kg / 0.00006479891,
        ground_threshold: inp.ground_threshold_m.unwrap_or(-1.0e6),
        ..BallisticInputs::default()
    };

    if inp.zero_drag {
        // All-zero Cd across a wide Mach span => drag acceleration is exactly zero, leaving gravity
        // as the sole force: an exact vacuum parabola. interpolate() clamps beyond the endpoints.
        b.custom_drag_table = Some(DragTable::new(
            vec![0.0, 1.0, 2.0, 3.0, 5.0, 10.0],
            vec![0.0; 6],
        ));
    }
    b
}

fn build_wind(inp: &Inputs) -> WindConditions {
    WindConditions {
        speed: inp.wind_speed_mps,
        direction: inp.wind_angle_rad,
        vertical_speed: inp.vertical_wind_mps,
    }
}

fn build_atmosphere(inp: &Inputs) -> AtmosphericConditions {
    AtmosphericConditions {
        temperature: inp.temperature_c,
        pressure: inp.pressure_hpa,
        // BallisticInputs.humidity is a 0-1 fraction; AtmosphericConditions.humidity is a percent.
        humidity: (inp.humidity_frac * 100.0).clamp(0.0, 100.0),
        altitude: inp.altitude_m,
    }
}

/// Interpolate one observable from the solved trajectory at a downrange distance `range_m`.
/// Points are stored in ascending downrange (x) order; we bracket `range_m` and interpolate
/// linearly in x. Returns an error string if the trajectory does not reach `range_m`.
fn observe(result: &TrajectoryResult, observable: &str, range_m: f64) -> Result<f64, String> {
    let pts = &result.points;
    if pts.len() < 2 {
        return Err(format!("trajectory has only {} point(s)", pts.len()));
    }
    let last_x = pts.last().unwrap().position.x;
    if range_m > last_x + 1e-6 {
        return Err(format!(
            "trajectory reached only {last_x:.2} m; cannot observe at {range_m:.2} m"
        ));
    }

    // Scalar of interest at a given trajectory point.
    let scalar = |p: &ballistics_engine::TrajectoryPoint| -> Result<f64, String> {
        match observable {
            "drop_m" => Ok(p.position.y),
            "drift_m" => Ok(p.position.z),
            "tof_s" => Ok(p.time),
            "velocity_mps" => Ok(p.velocity_magnitude),
            other => Err(format!("unknown observable: {other}")),
        }
    };

    // Find the bracketing segment [p0, p1] with p0.x <= range_m <= p1.x.
    for w in pts.windows(2) {
        let (p0, p1) = (&w[0], &w[1]);
        if p0.position.x <= range_m && range_m <= p1.position.x {
            let dx = p1.position.x - p0.position.x;
            let (v0, v1) = (scalar(p0)?, scalar(p1)?);
            if dx.abs() < 1e-12 {
                return Ok(v1);
            }
            let t = (range_m - p0.position.x) / dx;
            return Ok(v0 + t * (v1 - v0));
        }
    }
    // range_m within [0, last_x] but no bracket found only if it is below the first point's x.
    Err(format!(
        "no trajectory segment brackets {range_m:.2} m (first point x = {:.2} m)",
        pts[0].position.x
    ))
}

/// Effective tolerance = tol_abs + tol_rel*|expected|. At least one must be present. This is the
/// numpy `allclose` convention: tol_abs is the floor near zero, tol_rel scales with magnitude.
fn effective_tol(exp: &Expectation) -> f64 {
    let a = exp.tol_abs.unwrap_or(0.0);
    let r = exp.tol_rel.map(|r| r * exp.value.abs()).unwrap_or(0.0);
    a + r
}

struct Row {
    case_id: String,
    kind: String,
    observable: String,
    range_m: f64,
    actual: f64,
    expected: f64,
    delta: f64,
    tol: f64,
    known_discrepancy: bool,
}

fn load_case_files() -> Vec<PathBuf> {
    let dir = validation_dir();
    let mut files: Vec<PathBuf> = fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", dir.display()))
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("json"))
        .collect();
    files.sort();
    files
}

#[test]
fn golden_physics_validation() {
    let files = load_case_files();
    assert!(
        files.len() >= MIN_EXPECTED_CASES,
        "expected >= {} validation case files in {}, found {}",
        MIN_EXPECTED_CASES,
        validation_dir().display(),
        files.len(),
    );

    let mut hard_failures: Vec<Row> = Vec::new();
    let mut known_rows: Vec<Row> = Vec::new();
    let mut total_expectations = 0usize;

    // Set GOLDEN_DUMP=1 to print the delta table for EVERY expectation (pass or fail) — handy for
    // calibrating tolerances and for reporting head-room. Never affects pass/fail.
    let dump = std::env::var("GOLDEN_DUMP").is_ok();
    let mut all_rows: Vec<Row> = Vec::new();

    for file in &files {
        let text = fs::read_to_string(file)
            .unwrap_or_else(|e| panic!("cannot read {}: {e}", file.display()));
        let case: Case = serde_json::from_str(&text)
            .unwrap_or_else(|e| panic!("cannot parse {}: {e}", file.display()));

        // Basic dataset hygiene: every expectation must carry at least one tolerance.
        for exp in &case.expectations {
            assert!(
                exp.tol_abs.is_some() || exp.tol_rel.is_some(),
                "case '{}' observable '{}' @ {} m has no tolerance (need tol_abs and/or tol_rel)",
                case.id,
                exp.observable,
                exp.range_m,
            );
            assert!(
                !exp.tolerance_justification.trim().is_empty(),
                "case '{}' observable '{}' @ {} m is missing tolerance_justification",
                case.id,
                exp.observable,
                exp.range_m,
            );
        }

        let inputs = build_inputs(&case.inputs);
        let wind = build_wind(&case.inputs);
        let atmosphere = build_atmosphere(&case.inputs);

        let max_obs_range = case
            .expectations
            .iter()
            .map(|e| e.range_m)
            .fold(0.0_f64, f64::max);

        let mut solver = TrajectorySolver::new(inputs, wind, atmosphere);
        // Integrate comfortably past the furthest observation so interpolation always brackets it.
        solver.set_max_range(max_obs_range * 1.02 + 5.0);

        let result = solver
            .solve()
            .unwrap_or_else(|e| panic!("case '{}' failed to solve: {e:?}", case.id));

        for exp in &case.expectations {
            total_expectations += 1;
            let actual = observe(&result, &exp.observable, exp.range_m)
                .unwrap_or_else(|e| panic!("case '{}': {e}", case.id));
            let delta = actual - exp.value;
            let tol = effective_tol(exp);
            let within = delta.abs() <= tol;
            if dump {
                all_rows.push(Row {
                    case_id: case.id.clone(),
                    kind: case.source.kind.clone(),
                    observable: exp.observable.clone(),
                    range_m: exp.range_m,
                    actual,
                    expected: exp.value,
                    delta,
                    tol,
                    known_discrepancy: case.known_discrepancy,
                });
            }
            if !within || case.known_discrepancy {
                let row = Row {
                    case_id: case.id.clone(),
                    kind: case.source.kind.clone(),
                    observable: exp.observable.clone(),
                    range_m: exp.range_m,
                    actual,
                    expected: exp.value,
                    delta,
                    tol,
                    known_discrepancy: case.known_discrepancy,
                };
                if case.known_discrepancy {
                    known_rows.push(row);
                } else {
                    hard_failures.push(row);
                }
            }
        }

        if case.known_discrepancy && case.tracking_note.is_none() {
            panic!(
                "case '{}' is marked known_discrepancy but has no tracking_note",
                case.id
            );
        }
    }

    if dump {
        eprintln!(
            "\n=== GOLDEN_DUMP: all {} expectations across {} case files ===",
            total_expectations,
            files.len()
        );
        eprintln!("{}", render_table(&all_rows));
    }

    // Informational: surface known discrepancies (they never fail the suite).
    if !known_rows.is_empty() {
        eprintln!("\n=== KNOWN DISCREPANCIES (informational, excluded from pass/fail) ===");
        eprintln!("{}", render_table(&known_rows));
    }

    if !hard_failures.is_empty() {
        panic!(
            "\nGOLDEN PHYSICS VALIDATION FAILED: {} of {} expectations outside tolerance across {} case files\n{}\n\
             Investigate before adjusting tolerances: is the reference bad, the tolerance unjustified, \
             or is this a REAL physics discrepancy? A real one should be kept and marked \
             {{\"known_discrepancy\": true, \"tracking_note\": \"...\"}}.",
            hard_failures.len(),
            total_expectations,
            files.len(),
            render_table(&hard_failures),
        );
    }
}

/// Render an aligned delta table. Columns: case, kind, observable, range, actual, expected, delta,
/// tol, and how many multiples of tol the delta is (the "n×tol" over-run).
fn render_table(rows: &[Row]) -> String {
    let mut out = String::new();
    let header = format!(
        "{:<34} {:<20} {:<13} {:>9} {:>14} {:>14} {:>14} {:>12} {:>8}",
        "case", "kind", "observable", "range_m", "actual", "expected", "delta", "tol", "delta/tol",
    );
    out.push_str(&header);
    out.push('\n');
    out.push_str(&"-".repeat(header.len()));
    out.push('\n');
    for r in rows {
        let ratio = if r.tol > 0.0 {
            r.delta.abs() / r.tol
        } else {
            f64::INFINITY
        };
        let flag = if r.known_discrepancy { " [known]" } else { "" };
        out.push_str(&format!(
            "{:<34} {:<20} {:<13} {:>9.1} {:>14.5} {:>14.5} {:>+14.5} {:>12.5} {:>8.2}{}\n",
            truncate(&r.case_id, 34),
            truncate(&r.kind, 20),
            r.observable,
            r.range_m,
            r.actual,
            r.expected,
            r.delta,
            r.tol,
            ratio,
            flag,
        ));
    }
    out
}

fn truncate(s: &str, n: usize) -> String {
    if s.len() <= n {
        s.to_string()
    } else {
        format!("{}…", &s[..n.saturating_sub(1)])
    }
}
