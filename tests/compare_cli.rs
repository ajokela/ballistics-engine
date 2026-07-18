//! MBA-735: `compare` — side-by-side multi-load comparison. These tests run the real
//! binary and pin the CLI contract: load-spec parsing/validation, per-load independent
//! zeroing at shared conditions, physical ordering (higher BC ⇒ less drop/drift at
//! distance), JSON delta semantics (baseline deltas are zero), CSV shape and header
//! sanitization, and the PDF rejection.
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

fn out(args: &[&str]) -> Output {
    Command::new(BIN).args(args).output().expect("spawn ballistics")
}

fn ok_stdout(o: &Output) -> String {
    assert!(
        o.status.success(),
        "command failed: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    String::from_utf8_lossy(&o.stdout).into_owned()
}

fn err_text(o: &Output) -> String {
    assert!(!o.status.success(), "command unexpectedly succeeded");
    format!(
        "{}{}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    )
}

const LOAD_A: &str = "175 SMK:g7:0.243:175:2650";
const LOAD_B: &str = "168 ELD-M:g7:0.523:168:2700";

fn base_args<'a>(extra: &[&'a str]) -> Vec<&'a str> {
    let mut v = vec![
        "compare",
        "--load",
        LOAD_A,
        "--load",
        LOAD_B,
        "--zero-distance",
        "100",
        "--end",
        "500",
        "--step",
        "100",
    ];
    v.extend_from_slice(extra);
    v
}

#[test]
fn table_lists_both_loads_and_all_ranges() {
    let stdout = ok_stdout(&out(&base_args(&[])));
    assert!(stdout.contains("Load Comparison"), "missing title:\n{stdout}");
    assert!(stdout.contains("175 SMK"), "missing load A name:\n{stdout}");
    assert!(stdout.contains("168 ELD-M"), "missing load B name:\n{stdout}");
    for range in ["100", "200", "300", "400", "500"] {
        assert!(
            stdout.lines().any(|l| l.trim_start().starts_with(range)),
            "missing range row {range}:\n{stdout}"
        );
    }
}

#[test]
fn json_has_deltas_and_physical_ordering() {
    let stdout = ok_stdout(&out(&base_args(&["-o", "json"])));
    let v: serde_json::Value = serde_json::from_str(&stdout).expect("valid JSON");
    let cmp = &v["compare"];
    assert_eq!(cmp["loads"].as_array().unwrap().len(), 2);
    let rows = cmp["rows"].as_array().unwrap();
    assert_eq!(rows.len(), 5);

    let last = rows.last().unwrap()["loads"].as_array().unwrap().clone();
    let (a, b) = (&last[0], &last[1]);
    // Baseline (load #1) deltas are exactly zero
    assert_eq!(a["delta_drop"].as_f64().unwrap(), 0.0);
    assert_eq!(a["delta_energy"].as_f64().unwrap(), 0.0);
    // Higher-BC load drops and drifts less at 500 and keeps more velocity
    let (a_drop, b_drop) = (a["drop_adj"].as_f64().unwrap(), b["drop_adj"].as_f64().unwrap());
    let (a_drift, b_drift) =
        (a["drift_adj"].as_f64().unwrap(), b["drift_adj"].as_f64().unwrap());
    assert!(b_drop.abs() < a_drop.abs(), "high BC should drop less: {b_drop} vs {a_drop}");
    assert!(
        b_drift.abs() < a_drift.abs(),
        "high BC should drift less: {b_drift} vs {a_drift}"
    );
    assert!(b["velocity"].as_f64().unwrap() > a["velocity"].as_f64().unwrap());
    // Every numeric field is finite
    for load in &last {
        for key in
            ["drop", "drop_adj", "drift", "drift_adj", "velocity", "energy", "time", "delta_drop"]
        {
            assert!(load[key].as_f64().unwrap().is_finite(), "non-finite {key}");
        }
    }
}

#[test]
fn csv_has_sanitized_headers_and_row_per_range() {
    let o = out(&[
        "compare",
        "--load",
        "A:g7:0.243:175:2650",
        "--load",
        "B,x:g7:0.523:168:2700", // comma in the name must be sanitized in headers
        "--zero-distance",
        "100",
        "--end",
        "300",
        "--step",
        "100",
        "-o",
        "csv",
    ]);
    let stdout = ok_stdout(&o);
    let mut lines = stdout.lines();
    let header = lines.next().unwrap();
    assert!(header.starts_with("range_yd,"), "header: {header}");
    assert!(header.contains("A_drop_in"), "header: {header}");
    assert!(header.contains("B x_drop_in"), "comma not sanitized: {header}");
    // 1 range column + 2 loads x 7 fields
    assert_eq!(header.split(',').count(), 1 + 2 * 7, "header: {header}");
    assert_eq!(lines.count(), 3, "expected 3 data rows");
}

#[test]
fn metric_units_flow_through() {
    let stdout = ok_stdout(&out(&[
        "--units",
        "metric",
        "compare",
        "--load",
        "A:g7:0.243:11.34:808", // grams, m/s
        "--load",
        "B:g7:0.523:10.89:823",
        "--zero-distance",
        "100",
        "--end",
        "300",
        "--step",
        "100",
        "-o",
        "csv",
    ]));
    let header = stdout.lines().next().unwrap();
    assert!(header.starts_with("range_m,"), "header: {header}");
    assert!(header.contains("_drop_mm"), "header: {header}");
    assert!(header.contains("_velocity_m/s"), "header: {header}");
}

#[test]
fn rejects_fewer_than_two_loads() {
    let text = err_text(&out(&["compare", "--load", LOAD_A, "--zero-distance", "100"]));
    assert!(text.contains("at least 2 loads"), "got: {text}");
}

#[test]
fn rejects_malformed_specs_with_field_names() {
    // wrong field count
    let text = err_text(&out(&[
        "compare", "--load", "A:g7:0.243", "--load", LOAD_B, "--zero-distance", "100",
    ]));
    assert!(text.contains("NAME:DRAG:BC:MASS:VELOCITY"), "got: {text}");
    // bad drag model
    let text = err_text(&out(&[
        "compare", "--load", "A:g9:0.2:170:2600", "--load", LOAD_B, "--zero-distance", "100",
    ]));
    assert!(text.contains("must be g1 or g7"), "got: {text}");
    // non-numeric BC
    let text = err_text(&out(&[
        "compare", "--load", "A:g7:abc:170:2600", "--load", LOAD_B, "--zero-distance", "100",
    ]));
    assert!(text.contains("BC"), "got: {text}");
    // non-positive velocity
    let text = err_text(&out(&[
        "compare", "--load", "A:g7:0.2:170:0", "--load", LOAD_B, "--zero-distance", "100",
    ]));
    assert!(text.contains("VELOCITY"), "got: {text}");
    // start >= end
    let text = err_text(&out(&[
        "compare", "--load", LOAD_A, "--load", LOAD_B, "--zero-distance", "100", "--start",
        "600", "--end", "500",
    ]));
    assert!(text.contains("--start must be less than --end"), "got: {text}");
}

#[test]
fn rejects_pdf_output() {
    let text = err_text(&out(&base_args(&["-o", "pdf"])));
    assert!(text.contains("PDF output is not supported"), "got: {text}");
}

#[test]
fn optional_diameter_field_is_accepted() {
    // 6-field spec with explicit diameter parses and runs
    let stdout = ok_stdout(&out(&[
        "compare",
        "--load",
        "A:g7:0.243:175:2650:0.308",
        "--load",
        "B:g1:0.505:180:2600:0.338",
        "--zero-distance",
        "100",
        "--end",
        "300",
    ]));
    assert!(stdout.contains("dia 0.338"), "explicit diameter not echoed:\n{stdout}");
}

#[test]
fn profile_bc_segments_flow_into_compare() {
    // MBA-735 follow-up: a saved profile's velocity-BC segments must reach both the
    // zeroing and the sampled runs (the scalar-BC caveat is lifted for compare).
    // Hermetic: profiles live under $HOME/.ballistics, so point HOME at a temp dir,
    // save a scalar profile, then inject bc_segments into its JSON the way an .a7p
    // import would (profile save has no --bc-segment flag).
    use std::fs;
    let home = std::env::temp_dir().join(format!("compare-prof-{}", std::process::id()));
    fs::create_dir_all(&home).unwrap();

    let save = Command::new(BIN)
        .env("HOME", &home)
        .args([
            "profile", "save", "segtest", "-v", "2650", "-b", "0.243", "-m", "175", "-d",
            "0.308", "--drag-model", "g7",
        ])
        .output()
        .expect("spawn profile save");
    assert!(save.status.success(), "{}", String::from_utf8_lossy(&save.stderr));

    // Inject two velocity-BC segments (velocity_mps is always m/s in ProfileData):
    // markedly higher BC than the scalar 0.243 so the trajectories must diverge.
    let path = home.join(".ballistics").join("profiles").join("segtest.json");
    let mut profile: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();
    profile["bc_segments"] = serde_json::json!([
        {"bc": 0.35, "velocity_mps": 807.7},
        {"bc": 0.30, "velocity_mps": 500.0}
    ]);
    fs::write(&path, serde_json::to_string_pretty(&profile).unwrap()).unwrap();

    // Scalar twin via --load (identical numbers, no segments) vs the segmented profile.
    let o = Command::new(BIN)
        .env("HOME", &home)
        .args([
            "compare",
            "--load",
            "twin:g7:0.243:175:2650",
            "--profile",
            "segtest",
            "--zero-distance",
            "100",
            "--end",
            "600",
            "--step",
            "100",
            "-o",
            "json",
        ])
        .output()
        .expect("spawn compare");
    let stdout = ok_stdout(&o);
    let v: serde_json::Value = serde_json::from_str(&stdout).expect("valid JSON");
    let loads = v["compare"]["loads"].as_array().unwrap();
    assert_eq!(loads[0]["bc_segments"], false, "scalar twin mis-flagged");
    assert_eq!(loads[1]["bc_segments"], true, "profile segments not detected");

    // Segments (higher BC everywhere) must produce measurably less drop at 600.
    let last = v["compare"]["rows"].as_array().unwrap().last().unwrap()["loads"]
        .as_array()
        .unwrap()
        .clone();
    let twin_drop = last[0]["drop_adj"].as_f64().unwrap();
    let seg_drop = last[1]["drop_adj"].as_f64().unwrap();
    assert!(
        seg_drop.abs() < twin_drop.abs() * 0.98,
        "segments did not change the physics: twin {twin_drop} vs segmented {seg_drop}"
    );

    // And the table legend tags the segment-backed load.
    let table = ok_stdout(
        &Command::new(BIN)
            .env("HOME", &home)
            .args([
                "compare",
                "--load",
                "twin:g7:0.243:175:2650",
                "--profile",
                "segtest",
                "--zero-distance",
                "100",
                "--end",
                "300",
            ])
            .output()
            .expect("spawn compare table"),
    );
    assert!(table.contains("[BC segments]"), "table legend missing tag:\n{table}");

    fs::remove_dir_all(&home).ok();
}
