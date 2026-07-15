//! Integration tests for `trajectory --target-speed` (MBA-1325): the mover-ring hold
//! technique. Every trajectory point already knows its own time of flight, so the ring
//! radius a shooter watches for a mover — `target_speed x ToF(range)` around the hold
//! point — falls out of an already-solved trajectory as pure post-processing.
//!
//! Conventions locked by the ticket (see CLI_USAGE.md "Mover Ring"):
//! - `mover_ring_m` is ALWAYS meters (units in the name, MBA-1315 hygiene), regardless
//!   of `--units`.
//! - `mover_ring_mil = mover_ring_m / downrange_m * 1000`, omitted at the muzzle
//!   (downrange = 0, ratio undefined) — the table prints "-" there instead.
//! - The flag is additive: absent (or 0), table/JSON/CSV output is byte-identical to
//!   the pre-MBA-1325 contract (see `no_flag_*` tests below).

use serde_json::Value;
use std::path::PathBuf;
use std::process::Command;

fn get_cli_binary() -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("target");
    path.push("debug");
    path.push("ballistics");
    if !path.exists() {
        path.pop();
        path.pop();
        path.push("release");
        path.push("ballistics");
    }
    path
}

fn run(args: &[&str]) -> String {
    let out = Command::new(get_cli_binary())
        .args(args)
        .output()
        .expect("run");
    assert!(
        out.status.success(),
        "command failed ({:?}): {}",
        args,
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).to_string()
}

fn run_json(args: &[&str]) -> Value {
    serde_json::from_str(&run(args)).expect("json")
}

const BASE_METRIC: &[&str] = &[
    "--units", "metric", "trajectory", "-v", "823", "-b", "0.475", "-m", "10.9", "-d", "7.82",
    "--max-range", "450",
];

// ---------------------------------------------------------------------------------
// JSON: exact ring math (metric run, so the JSON `z` downrange field is already in
// meters and needs no unit round-trip against the target-speed conversion).
// ---------------------------------------------------------------------------------

/// `mover_ring_m` is EXACTLY `target_speed_mps * point_time` at every printed point.
#[test]
fn json_ring_m_equals_target_speed_times_time() {
    let target_speed_mps = 3.0_f64; // metric units: --target-speed is already m/s
    let mut args: Vec<&str> = BASE_METRIC.to_vec();
    args.extend(["--full", "-o", "json", "--target-speed", "3"]);
    let json = run_json(&args);
    let points = json["trajectory"].as_array().expect("trajectory array");
    assert!(points.len() > 5, "expected several points, got {}", points.len());

    let mut checked = 0;
    for p in points {
        let time = p["time"].as_f64().expect("time");
        let ring_m = p["mover_ring_m"].as_f64().expect("mover_ring_m present");
        let expected = target_speed_mps * time;
        assert!(
            (ring_m - expected).abs() < 1e-9,
            "ring_m {ring_m} != target_speed_mps*time {expected} at t={time}"
        );
        checked += 1;
    }
    assert!(checked > 5);
}

/// `mover_ring_mil = mover_ring_m / downrange_m * 1000`, and both fields are OMITTED
/// entirely at the muzzle (downrange = 0), not present-as-zero/NaN.
#[test]
fn json_ring_mil_matches_formula_and_omits_at_muzzle() {
    let mut args: Vec<&str> = BASE_METRIC.to_vec();
    args.extend(["--full", "-o", "json", "--target-speed", "3"]);
    let json = run_json(&args);
    let points = json["trajectory"].as_array().expect("trajectory array");

    // First point is the muzzle: z (downrange, metric = meters) must be 0.
    let muzzle = &points[0];
    assert_eq!(muzzle["z"].as_f64().unwrap(), 0.0, "first point should be the muzzle");
    assert_eq!(
        muzzle["mover_ring_m"].as_f64().unwrap(),
        0.0,
        "mover_ring_m at the muzzle is 0 (target_speed * t=0), not omitted"
    );
    assert!(
        muzzle.get("mover_ring_mil").is_none(),
        "mover_ring_mil must be OMITTED at the muzzle (0/0 is degenerate), got {:?}",
        muzzle.get("mover_ring_mil")
    );

    let mut checked = 0;
    for p in points.iter().skip(1) {
        let downrange_m = p["z"].as_f64().expect("z (downrange, metric = meters)");
        if downrange_m <= 0.0 {
            continue;
        }
        let ring_m = p["mover_ring_m"].as_f64().expect("mover_ring_m");
        let ring_mil = p["mover_ring_mil"]
            .as_f64()
            .expect("mover_ring_mil present past the muzzle");
        let expected = ring_m / downrange_m * 1000.0;
        assert!(
            (ring_mil - expected).abs() < 1e-6,
            "ring_mil {ring_mil} != ring_m/downrange_m*1000 {expected} at downrange={downrange_m}"
        );
        checked += 1;
    }
    assert!(checked > 5, "expected several non-muzzle points checked, got {checked}");
}

/// Imperial (default units): --target-speed is mph, converted to m/s before the ring
/// math — same formula, different input unit.
#[test]
fn json_ring_m_imperial_uses_mph_conversion() {
    let target_speed_mph = 5.0_f64;
    let target_speed_mps = target_speed_mph * 0.44704;
    let out = run_json(&[
        "trajectory", "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308", "--max-range",
        "500", "--full", "-o", "json", "--target-speed", "5",
    ]);
    let points = out["trajectory"].as_array().expect("trajectory array");
    let last = points.last().expect("at least one point");
    let time = last["time"].as_f64().expect("time");
    let ring_m = last["mover_ring_m"].as_f64().expect("mover_ring_m");
    let expected = target_speed_mps * time;
    assert!(
        (ring_m - expected).abs() < 1e-9,
        "ring_m {ring_m} != target_speed_mps*time {expected}"
    );
}

// ---------------------------------------------------------------------------------
// Table: header, muzzle "-", and cross-checked ring math (metric so X is meters).
// ---------------------------------------------------------------------------------

/// Parse the "Trajectory Points:" table's data rows into (time_s, x_m, ring_mil).
/// `ring_mil` is `None` when the printed cell is the literal "-" placeholder.
fn parse_table_ring_rows(table: &str) -> Vec<(f64, f64, Option<f64>)> {
    let mut rows = Vec::new();
    let mut in_section = false;
    for line in table.lines() {
        if line.trim_end() == "Trajectory Points:" {
            in_section = true;
            continue;
        }
        if !in_section {
            continue;
        }
        if line.starts_with('└') {
            break;
        }
        if !line.starts_with('│') {
            continue; // border (┌/├) or blank lines
        }
        let cells: Vec<&str> = line
            .split('│')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .collect();
        if cells.len() != 6 {
            continue;
        }
        // The header row's cells ("Time (s)", "X (m)", ...) don't parse as f64, so
        // this naturally skips it without a separate check.
        if let (Ok(time), Ok(x)) = (cells[0].parse::<f64>(), cells[1].parse::<f64>()) {
            let ring = if cells[5] == "-" {
                None
            } else {
                Some(cells[5].parse::<f64>().expect("ring cell should be numeric or '-'"))
            };
            rows.push((time, x, ring));
        }
    }
    rows
}

#[test]
fn table_has_ring_column_with_correct_header_and_muzzle_dash() {
    let mut args: Vec<&str> = BASE_METRIC.to_vec();
    args.extend(["--full", "-o", "table", "--target-speed", "3"]);
    let out = run(&args);
    assert!(out.contains("Ring(mil)"), "table should carry a Ring(mil) header:\n{out}");

    let rows = parse_table_ring_rows(&out);
    assert!(rows.len() > 3, "expected several parsed rows, got {}", rows.len());
    let (t0, x0, ring0) = rows[0];
    assert!(t0.abs() < 1e-9 && x0.abs() < 1e-9, "first row should be the muzzle (t=0, x=0)");
    assert_eq!(ring0, None, "muzzle row must print '-' for Ring, not a number");

    // The table prints time to 3dp and X to 2dp, so re-deriving target_speed*t/x*1000
    // from the DISPLAYED text (not the full-precision values the engine actually used)
    // carries its own rounding error — amplified at small x, where ring_mil ~ 1/x. The
    // bit-exact formula is already covered by the JSON tests above; this is a coarser
    // cross-check that the printed column is in the right ballpark, restricted to points
    // far enough downrange that text rounding is negligible.
    let mut checked = 0;
    for (t, x, ring) in rows.iter().skip(1) {
        if *x < 50.0 {
            continue;
        }
        let ring = ring.expect("non-muzzle row should have a numeric Ring value");
        let expected = 3.0 * t / x * 1000.0;
        let tol = (expected.abs() * 0.01).max(0.02);
        assert!(
            (ring - expected).abs() < tol,
            "table ring {ring} != target_speed*t/x*1000 {expected} (tol {tol}) at t={t} x={x}"
        );
        checked += 1;
    }
    assert!(checked > 2, "expected several far-downrange rows checked, got {checked}");
}

// ---------------------------------------------------------------------------------
// CSV: extra column, header carries the unit, correct value.
// ---------------------------------------------------------------------------------

#[test]
fn csv_raw_points_gains_ring_mil_column_when_flag_set() {
    let mut args: Vec<&str> = BASE_METRIC.to_vec();
    args.extend(["--full", "-o", "csv", "--target-speed", "3"]);
    let out = run(&args);
    let mut lines = out.lines();
    let header = lines.next().expect("header line");
    assert!(header.ends_with(",ring_mil"), "CSV header should end with ,ring_mil: {header}");

    // As in the table test above: CSV prints x/z to 2dp, so re-deriving the ratio from
    // that rounded text has its own error, amplified at small downrange. Restrict the
    // formula cross-check to far-enough-downrange rows; the bit-exact check already
    // lives in the JSON tests.
    let mut checked = 0;
    for line in lines {
        let cols: Vec<&str> = line.split(',').collect();
        assert_eq!(cols.len(), 7, "expected 7 columns (6 base + ring_mil): {line}");
        let time: f64 = cols[0].parse().expect("time");
        let downrange_m: f64 = cols[3].parse().expect("z column (downrange, meters)");
        if downrange_m <= 0.0 {
            assert_eq!(cols[6], "", "muzzle row's ring_mil field should be empty: {line}");
            continue;
        }
        let ring_mil: f64 = cols[6].parse().expect("ring_mil");
        if downrange_m < 50.0 {
            continue;
        }
        let expected = 3.0 * time / downrange_m * 1000.0;
        let tol = (expected.abs() * 0.01).max(0.02);
        assert!(
            (ring_mil - expected).abs() < tol,
            "csv ring_mil {ring_mil} != formula {expected} (tol {tol}) on line: {line}"
        );
        checked += 1;
    }
    assert!(checked > 3, "expected several non-muzzle rows checked, got {checked}");
}

#[test]
fn csv_sampled_points_also_gains_ring_mil_column() {
    let mut args: Vec<&str> = BASE_METRIC.to_vec();
    args.extend(["--full", "--sample-trajectory", "-o", "csv", "--target-speed", "3"]);
    let out = run(&args);
    let mut lines = out.lines();
    let header = lines.next().expect("header line");
    assert!(
        header.ends_with(",ring_mil"),
        "sampled CSV header should also end with ,ring_mil: {header}"
    );
    let row = lines.next().expect("at least one sampled row");
    assert_eq!(row.split(',').count(), 7, "expected 7 columns in sampled CSV row: {row}");
}

// ---------------------------------------------------------------------------------
// Additive / no-regression: absent (or zero) --target-speed leaves every format
// byte-for-byte identical to the pre-MBA-1325 contract.
// ---------------------------------------------------------------------------------

#[test]
fn no_flag_table_has_no_ring_column() {
    let out = run(&[
        "trajectory", "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308", "--full", "-o",
        "table",
    ]);
    assert!(!out.contains("Ring"), "Ring column must not appear without --target-speed");
    // The original 5-column border must be exactly preserved (byte-identical width).
    assert!(out.contains("┌──────────┬──────────┬──────────┬──────────┬──────────┐"));
    assert!(!out.contains("┌──────────┬──────────┬──────────┬──────────┬──────────┬──────────┐"));
}

#[test]
fn no_flag_json_has_no_mover_ring_fields() {
    let out = run(&[
        "trajectory", "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308", "--full", "-o",
        "json",
    ]);
    assert!(!out.contains("mover_ring"), "mover_ring fields must not appear without --target-speed");
}

#[test]
fn zero_target_speed_behaves_identically_to_omitted() {
    let with_zero = run(&[
        "trajectory", "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308", "--full", "-o",
        "json", "--target-speed", "0",
    ]);
    let omitted = run(&[
        "trajectory", "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308", "--full", "-o",
        "json",
    ]);
    assert_eq!(with_zero, omitted, "--target-speed 0 must match omitting the flag entirely");
}

#[test]
fn no_flag_csv_header_is_unchanged() {
    let out = run(&[
        "trajectory", "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308", "--full", "-o",
        "csv",
    ]);
    let header = out.lines().next().expect("header");
    assert_eq!(header, "time_s,x_yd,y_yd,z_yd,velocity_fps,energy_ft-lb");
}

// ---------------------------------------------------------------------------------
// Review fix M1: --target-speed shares lead's f64_range(0.0, 300.0) clap bound.
// ---------------------------------------------------------------------------------

#[test]
fn trajectory_rejects_out_of_range_target_speed() {
    // "=" syntax so the negative case reaches the value parser rather than being
    // eaten earlier as an unexpected hyphen argument (also a rejection, but this
    // asserts the f64_range bound specifically).
    for bad in ["--target-speed=301", "--target-speed=1000000000", "--target-speed=-1"] {
        let out = Command::new(get_cli_binary())
            .args([
                "trajectory", "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308", bad,
            ])
            .output()
            .expect("run");
        assert!(
            !out.status.success(),
            "{bad} must be rejected, but the command succeeded"
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains("is not in range"),
            "expected the f64_range clap error for {bad}, got:\n{stderr}"
        );
    }
}

// ---------------------------------------------------------------------------------
// Review fix M3: the ring TABLE column honors --adjustment-unit; the MOA figure uses
// the crate's locked printed-table dial constant (MBA-724), exactly 3.438 x mil.
// ---------------------------------------------------------------------------------

#[test]
fn table_ring_column_honors_moa_adjustment_unit() {
    // --target-speed 300 (the bound) makes ring_mil large enough (~190 at the last
    // point) that the dial constant is DISCRIMINABLE from the exact-angle 3437.7467
    // variant at 2 printed decimals: the two differ by mil x 6.57e-5 > 0.012, more
    // than the 0.01 print quantum, so a wrong constant cannot round to the same cell.
    let base = [
        "trajectory", "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308", "--max-range",
        "400", "--target-speed", "300", "--full",
    ];

    // Full-precision mil from JSON (JSON is adjustment-unit-invariant by contract).
    let mut json_args = base.to_vec();
    json_args.extend(["-o", "json"]);
    let json = run_json(&json_args);
    let last = json["trajectory"]
        .as_array()
        .expect("points")
        .last()
        .expect("last")
        .clone();
    let mil = last["mover_ring_mil"].as_f64().expect("mover_ring_mil");
    assert!(
        mil > 160.0,
        "sensitivity precondition: need ring_mil > 160 to pin the constant, got {mil}"
    );

    // MOA table: header flips to Ring(moa); the last row's cell must equal the
    // 2dp-rounded mil x 3.438 exactly (string match pins the constant).
    let mut moa_args = base.to_vec();
    moa_args.extend(["--adjustment-unit", "moa", "-o", "table"]);
    let moa_table = run(&moa_args);
    assert!(moa_table.contains("Ring(moa)"), "moa run must show Ring(moa):\n{moa_table}");
    assert!(!moa_table.contains("Ring(mil)"), "moa run must not show Ring(mil)");
    let moa_rows = parse_table_ring_rows(&moa_table);
    let (_, _, moa_cell) = moa_rows.last().expect("moa table rows");
    let moa_cell = moa_cell.expect("last row has a numeric ring cell");
    assert_eq!(
        format!("{:.2}", moa_cell),
        format!("{:.2}", mil * 3.438),
        "Ring(moa) cell must be exactly mil x 3.438 (the MBA-724 dial constant), 2dp-rounded"
    );

    // MIL table (default): same last row prints the mil value itself.
    let mut mil_args = base.to_vec();
    mil_args.extend(["-o", "table"]);
    let mil_table = run(&mil_args);
    assert!(mil_table.contains("Ring(mil)"));
    let mil_rows = parse_table_ring_rows(&mil_table);
    let (_, _, mil_cell) = mil_rows.last().expect("mil table rows");
    assert_eq!(
        format!("{:.2}", mil_cell.expect("numeric ring cell")),
        format!("{:.2}", mil),
        "Ring(mil) cell must match the JSON mover_ring_mil, 2dp-rounded"
    );
}
