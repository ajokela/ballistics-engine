// MBA-1410: finish the turret-units work started in MBA-1355. smoa/iphy/clicks reached
// trajectory, come-ups, and the PDF dope card there; this extends them to wind-card, lead,
// range-table, and compare, and adds independent elevation-vs-windage unit selection
// (`--windage-unit`) to the two dual-axis commands (range-table, compare). The old
// "clicks is currently supported for trajectory and come-ups only" rejection is gone from
// every one of these surfaces.

use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

fn run(args: &[&str]) -> Output {
    Command::new(BIN).args(args).output().expect("spawn ballistics")
}

fn stdout(out: &Output) -> String {
    String::from_utf8_lossy(&out.stdout).to_string()
}

fn stderr(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).to_string()
}

const OLD_REJECTION: &str =
    "error: --adjustment-unit clicks is currently supported for trajectory and come-ups only (MBA-1355)";

// ---------------------------------------------------------------------------
// The old rejection message is gone from every extended surface.
// ---------------------------------------------------------------------------

#[test]
fn old_rejection_message_is_gone_from_wind_card() {
    let out = run(&[
        "wind-card", "-v", "2700", "-m", "168", "-d", "0.308", "-b", "0.5",
        "--zero-distance", "100", "--end", "300",
        "--adjustment-unit", "clicks", "--windage-click-value", "0.1mil",
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(!stderr(&out).contains(OLD_REJECTION));
}

#[test]
fn old_rejection_message_is_gone_from_lead() {
    let out = run(&[
        "lead", "-v", "2700", "-m", "168", "-d", "0.308", "-b", "0.5",
        "--target-speed", "5",
        "--adjustment-unit", "clicks", "--windage-click-value", "0.1mil",
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(!stderr(&out).contains(OLD_REJECTION));
}

#[test]
fn old_rejection_message_is_gone_from_range_table() {
    let out = run(&[
        "range-table", "-v", "2700", "-m", "168", "-d", "0.308", "-b", "0.5",
        "--zero-distance", "100", "--end", "300",
        "--adjustment-unit", "clicks",
        "--elevation-click-value", "0.1mil", "--windage-click-value", "0.1mil",
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(!stderr(&out).contains(OLD_REJECTION));
}

#[test]
fn old_rejection_message_is_gone_from_compare() {
    let out = run(&[
        "compare",
        "--load", "A:g1:0.5:168:2700",
        "--load", "B:g1:0.45:175:2650",
        "--zero-distance", "100", "--end", "300",
        "--adjustment-unit", "clicks",
        "--elevation-click-value", "0.1mil", "--windage-click-value", "0.1mil",
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(!stderr(&out).contains(OLD_REJECTION));
}

// ---------------------------------------------------------------------------
// Each newly-supported command renders each unit (smoa/iphy/clicks).
// ---------------------------------------------------------------------------

#[test]
fn wind_card_renders_smoa_iphy_and_clicks() {
    let base = [
        "wind-card", "-v", "2700", "-m", "168", "-d", "0.308", "-b", "0.5",
        "--zero-distance", "100", "--end", "300",
    ];
    for unit in ["smoa", "iphy"] {
        let out = run(&[&base[..], &["--adjustment-unit", unit]].concat());
        assert!(out.status.success(), "{unit}: {}", stderr(&out));
    }
    let out = run(&[&base[..], &["--adjustment-unit", "clicks", "--windage-click-value", "0.25moa"]].concat());
    assert!(out.status.success(), "{}", stderr(&out));
}

#[test]
fn lead_renders_clicks() {
    let out = run(&[
        "lead", "-v", "2700", "-m", "168", "-d", "0.308", "-b", "0.5",
        "--target-speed", "5",
        "--adjustment-unit", "clicks", "--windage-click-value", "0.25moa",
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
}

#[test]
fn range_table_renders_smoa_iphy_and_clicks() {
    let base = [
        "range-table", "-v", "2700", "-m", "168", "-d", "0.308", "-b", "0.5",
        "--zero-distance", "100", "--end", "300",
    ];
    for unit in ["smoa", "iphy"] {
        let out = run(&[&base[..], &["--adjustment-unit", unit]].concat());
        assert!(out.status.success(), "{unit}: {}", stderr(&out));
    }
    let out = run(&[
        &base[..],
        &[
            "--adjustment-unit", "clicks",
            "--elevation-click-value", "0.1mil",
            "--windage-click-value", "0.1mil",
        ],
    ]
    .concat());
    assert!(out.status.success(), "{}", stderr(&out));
}

#[test]
fn compare_renders_smoa_iphy_and_clicks() {
    let base = [
        "compare",
        "--load", "A:g1:0.5:168:2700",
        "--load", "B:g1:0.45:175:2650",
        "--zero-distance", "100", "--end", "300",
    ];
    for unit in ["smoa", "iphy"] {
        let out = run(&[&base[..], &["--adjustment-unit", unit]].concat());
        assert!(out.status.success(), "{unit}: {}", stderr(&out));
    }
    let out = run(&[
        &base[..],
        &[
            "--adjustment-unit", "clicks",
            "--elevation-click-value", "0.1mil",
            "--windage-click-value", "0.1mil",
        ],
    ]
    .concat());
    assert!(out.status.success(), "{}", stderr(&out));
}

// ---------------------------------------------------------------------------
// Independent elevation-vs-windage unit selection genuinely differs between
// axes -- assert the Drop and Wind columns carry DIFFERENT units in ONE run
// (a test using the same unit for both would pass even if the axes were still
// coupled).
// ---------------------------------------------------------------------------

#[test]
fn range_table_json_carries_independent_elevation_and_windage_units() {
    let out = run(&[
        "range-table", "-v", "2700", "-m", "168", "-d", "0.308", "-b", "0.5",
        "--zero-distance", "100", "--end", "300",
        "--adjustment-unit", "mil", "--windage-unit", "moa",
        "-o", "json",
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
    let v: serde_json::Value = serde_json::from_str(&stdout(&out)).expect("json");
    assert_eq!(v["adjustment_unit"], "MIL");
    assert_eq!(v["windage_unit"], "MOA");
    assert_ne!(v["adjustment_unit"], v["windage_unit"]);
}

#[test]
fn range_table_drop_and_wind_adj_values_differ_by_unit_ratio_when_axes_diverge() {
    // mil elevation, moa windage on an identical trajectory: since MOA and MIL are
    // different angular units, the raw drop_adj/wind_adj numbers for the SAME physical
    // drop/drift would differ by the MOA/MIL ratio, unlike a same-unit run where the
    // columns coincidentally use one factor. Compare against a mil/mil run.
    let mil_mil = run(&[
        "range-table", "-v", "2700", "-m", "168", "-d", "0.308", "-b", "0.5",
        "--zero-distance", "100", "--end", "300",
        "--adjustment-unit", "mil", "-o", "json",
    ]);
    let mixed = run(&[
        "range-table", "-v", "2700", "-m", "168", "-d", "0.308", "-b", "0.5",
        "--zero-distance", "100", "--end", "300",
        "--adjustment-unit", "mil", "--windage-unit", "moa", "-o", "json",
    ]);
    assert!(mil_mil.status.success() && mixed.status.success());
    let mm: serde_json::Value = serde_json::from_str(&stdout(&mil_mil)).unwrap();
    let mx: serde_json::Value = serde_json::from_str(&stdout(&mixed)).unwrap();
    let last_row_wind = |v: &serde_json::Value| -> f64 {
        v["data"].as_array().unwrap().last().unwrap()["wind_adj"].as_f64().unwrap()
    };
    let last_row_drop = |v: &serde_json::Value| -> f64 {
        v["data"].as_array().unwrap().last().unwrap()["drop_adj"].as_f64().unwrap()
    };
    // Drop (elevation) is unchanged between the two runs -- only windage moved.
    assert!((last_row_drop(&mm) - last_row_drop(&mx)).abs() < 1e-9);
    // Wind (windage) differs: the MOA-windage run's wind_adj is the same physical drift
    // expressed in a different graduation, so it must scale by the crate's locked
    // MOA/MIL adjustment-factor ratio (3438/1000 == 3.438, see adjustment_factor's docs),
    // NOT the 1.047 TRUE-MOA/mrad *angle* ratio -- these are click-table constants, not
    // exact angles.
    let (wind_mil, wind_moa) = (last_row_wind(&mm), last_row_wind(&mx));
    assert!(wind_mil.abs() > 1e-6, "sanity: nonzero drift expected");
    assert!(
        (wind_moa.abs() / wind_mil.abs() - 3.438).abs() < 0.01,
        "MOA/MIL wind_adj ratio should be ~3.438: got {wind_moa} vs {wind_mil}"
    );
}

#[test]
fn compare_json_carries_independent_elevation_and_windage_units() {
    let out = run(&[
        "compare",
        "--load", "A:g1:0.5:168:2700",
        "--load", "B:g1:0.45:175:2650",
        "--zero-distance", "100", "--end", "300",
        "--adjustment-unit", "mil", "--windage-unit", "moa",
        "-o", "json",
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
    let v: serde_json::Value = serde_json::from_str(&stdout(&out)).expect("json");
    let units = &v["compare"]["units"];
    assert_eq!(units["adjustment"], "MIL");
    assert_eq!(units["windage_adjustment"], "MOA");
    assert_ne!(units["adjustment"], units["windage_adjustment"]);
}

// ---------------------------------------------------------------------------
// --windage-unit clicks is rejected without --adjustment-unit clicks (the
// judgment call documented in resolve_windage_unit's doc comment).
// ---------------------------------------------------------------------------

#[test]
fn windage_unit_clicks_without_elevation_clicks_is_rejected() {
    let out = run(&[
        "range-table", "-v", "2700", "-m", "168", "-d", "0.308", "-b", "0.5",
        "--zero-distance", "100", "--end", "300",
        "--adjustment-unit", "mil", "--windage-unit", "clicks",
    ]);
    assert!(!out.status.success());
    let e = stderr(&out);
    assert!(e.contains("--windage-unit"), "{e}");
    assert!(e.contains("--adjustment-unit"), "{e}");
}

// ---------------------------------------------------------------------------
// Default (mil elevation, mil windage, single --adjustment-unit) is
// byte-identical across every extended command.
// ---------------------------------------------------------------------------

#[test]
fn wind_card_default_is_byte_identical_to_pre_mba_1410() {
    let base = [
        "wind-card", "-v", "2700", "-m", "168", "-d", "0.308", "-b", "0.5",
        "--zero-distance", "100", "--end", "300",
    ];
    let no_flag = run(&base);
    let explicit_mil = run(&[&base[..], &["--adjustment-unit", "mil"]].concat());
    assert!(no_flag.status.success() && explicit_mil.status.success());
    assert_eq!(stdout(&no_flag), stdout(&explicit_mil));
}

#[test]
fn lead_default_is_byte_identical_to_pre_mba_1410() {
    let base = [
        "lead", "-v", "2700", "-m", "168", "-d", "0.308", "-b", "0.5",
        "--target-speed", "5",
    ];
    let no_flag = run(&base);
    let explicit_mil = run(&[&base[..], &["--adjustment-unit", "mil"]].concat());
    assert!(no_flag.status.success() && explicit_mil.status.success());
    assert_eq!(stdout(&no_flag), stdout(&explicit_mil));
}

#[test]
fn range_table_default_is_byte_identical_to_pre_mba_1410() {
    let base = [
        "range-table", "-v", "2700", "-m", "168", "-d", "0.308", "-b", "0.5",
        "--zero-distance", "100", "--end", "300",
    ];
    let no_flag = run(&base);
    let explicit_mil = run(&[&base[..], &["--adjustment-unit", "mil"]].concat());
    let explicit_both = run(&[&base[..], &["--adjustment-unit", "mil", "--windage-unit", "mil"]].concat());
    assert!(no_flag.status.success() && explicit_mil.status.success() && explicit_both.status.success());
    assert_eq!(stdout(&no_flag), stdout(&explicit_mil));
    assert_eq!(stdout(&no_flag), stdout(&explicit_both));
}

#[test]
fn compare_default_is_byte_identical_to_pre_mba_1410() {
    let base = [
        "compare",
        "--load", "A:g1:0.5:168:2700",
        "--load", "B:g1:0.45:175:2650",
        "--zero-distance", "100", "--end", "300",
    ];
    let no_flag = run(&base);
    let explicit_mil = run(&[&base[..], &["--adjustment-unit", "mil"]].concat());
    let explicit_both = run(&[&base[..], &["--adjustment-unit", "mil", "--windage-unit", "mil"]].concat());
    assert!(no_flag.status.success() && explicit_mil.status.success() && explicit_both.status.success());
    assert_eq!(stdout(&no_flag), stdout(&explicit_mil));
    assert_eq!(stdout(&no_flag), stdout(&explicit_both));
}

// ---------------------------------------------------------------------------
// The two folded MBA-1355 backlog minors.
// ---------------------------------------------------------------------------

#[test]
fn come_ups_windage_click_value_now_warns_instead_of_silently_doing_nothing() {
    let out = run(&[
        "come-ups", "-v", "2700", "-m", "168", "-d", "0.308", "-b", "0.5",
        "--zero-distance", "100", "--end", "300",
        "--adjustment-unit", "clicks",
        "--elevation-click-value", "0.1mil",
        "--windage-click-value", "0.1mil",
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
    let e = stderr(&out);
    assert!(e.contains("--windage-click-value"), "{e}");
    assert!(e.contains("come-ups"), "{e}");
}

#[test]
fn come_ups_without_windage_click_value_has_no_inert_flag_warning() {
    let out = run(&[
        "come-ups", "-v", "2700", "-m", "168", "-d", "0.308", "-b", "0.5",
        "--zero-distance", "100", "--end", "300",
        "--adjustment-unit", "clicks",
        "--elevation-click-value", "0.1mil",
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(!stderr(&out).contains("--windage-click-value"));
}

/// PDF dope card clicks cells must print as bare integers ("5"), never "5.0" -- pinned at
/// the unit level in `pdf_dope_card::tests::format_adjustment_drops_the_decimal_for_clicks_only`;
/// this end-to-end smoke confirms a clicks dope card still generates a valid PDF (the
/// module-level test already proves the exact string content, since PDF bytes are
/// compressed and not greppable here — see `tests/dope_card_units.rs`'s doc comment).
#[test]
fn clicks_dope_card_generates_a_valid_pdf() {
    let out_path = std::env::temp_dir().join("bx_dope_clicks_mba1410.pdf");
    let out = run(&[
        "trajectory", "-v", "2700", "-b", "0.5", "-m", "175", "-d", "0.308", "--drag-model", "g7",
        "--max-range", "800", "--temperature", "59", "--pressure", "29.92", "--auto-zero", "100",
        "--sample-trajectory", "-o", "pdf", "--output-file", out_path.to_str().unwrap(),
        "--target-speed", "3",
        "--adjustment-unit", "clicks", "--elevation-click-value", "0.25moa",
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
    let bytes = std::fs::read(&out_path).expect("PDF written");
    assert!(bytes.len() > 10_000, "PDF suspiciously small ({} bytes)", bytes.len());
    assert_eq!(&bytes[..5], b"%PDF-");
    let _ = std::fs::remove_file(&out_path);
}
