// Plan B Task 2: byte-for-byte stdout pins for the four card surfaces (come-ups,
// range-table, wind-card, compare) in table/csv/json form, taken BEFORE the
// `adjustment_display`/`windage_adjustment_display` dial-boundary change that adds
// residual-returning quantization (src/adjustment.rs's `quantize_angle`). These are the
// safety net for that change and for the later CardRow consolidation (Plan B Tasks 9/10):
// if any of those refactors accidentally move a printed number, one of these 12 tests
// fails with a full-string diff instead of the drift going unnoticed.
//
// Every invocation below is copied from an existing CLI test (not invented here):
// - come-ups: flags from `tests/sweep_step_cli.rs`'s `base_args` (-v/-m/-d/-b/--zero-distance),
//   ranges from this file's own `--start 100 --end 300 --step 100` (the pattern the Step-1
//   brief itself names) for exactly 3 rows.
// - range-table: same base flags plus `--end 300 --step 100` (`tests/adjustment_units_extended.rs`
//   uses `--zero-distance 100 --end 300` for range-table; `--step 100` is added, matching
//   `tests/compare_cli.rs`'s own use of `--step 100`, to land exactly on 3 rows instead of 5).
// - wind-card: `--zero-distance 100 --end 300`, byte-for-byte the base array in
//   `tests/adjustment_units_extended.rs`'s `wind_card_renders_smoa_iphy_and_clicks` (default
//   step is 100, so this is already exactly 3 rows).
// - compare: `--load`/`--zero-distance`/`--end 300 --step 100` from `tests/compare_cli.rs`'s
//   `csv_has_sanitized_headers_and_row_per_range` (which asserts exactly 3 data rows), with
//   `LOAD_A`/`LOAD_B` from that same file's constants.
//
// If a change to any of these four surfaces is INTENTIONAL, regenerate the affected
// golden file(s) with the same invocation (see each test) and say so in the commit
// message — silent drift is exactly what this file exists to catch.
//
// Regenerated in 0.50.0 (the unified aloft vapour model): these runs fly the
// default 50 % RH, and the air off the station's height now holds the station's water-vapour
// mole fraction with a virtual-temperature lapse. Only the four `*_json.txt` files moved (the
// table and CSV forms round the change away), every value by at most 1.6e-6 relative and
// 1.6e-9 absolute in the units printed. Old -> new, the largest of each file:
// - come_ups_json: data[2].drop 1.2305099621071154 -> 1.2305099671254405 (+4.1e-9), data[0].drop
//   0.00042704099043231494 -> 0.0004270414041071023 (+9.7e-7 of a near-zero value);
// - range_table_json: data[2].drop 13.289507590756847 -> 13.289507644954757 (+4.1e-9), data[1]
//   wind_drift -2.8520660388064454 -> -2.852066094782608 (+2.0e-8);
// - wind_card_json: data[1].wind_20 -0.7922624146772185 -> -0.7922624302262841 (+2.0e-8);
// - compare_json: rows[2].loads[0].energy 1743.5998823201753 -> 1743.5999063677366 (+1.4e-8),
//   rows[0].loads[1].drop -0.0003852813726742867 -> -0.0003852807537905148 (-1.6e-6 of a
//   near-zero value).
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_ballistics")
}

fn stdout_of(args: &[&str]) -> String {
    let out = Command::new(bin()).args(args).output().expect("run");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).expect("utf8")
}

const COME_UPS: &[&str] = &[
    "come-ups",
    "-v",
    "2700",
    "-m",
    "168",
    "-d",
    "0.308",
    "-b",
    "0.5",
    "--zero-distance",
    "100",
    "--start",
    "100",
    "--end",
    "300",
    "--step",
    "100",
];

const RANGE_TABLE: &[&str] = &[
    "range-table",
    "-v",
    "2700",
    "-m",
    "168",
    "-d",
    "0.308",
    "-b",
    "0.5",
    "--zero-distance",
    "100",
    "--end",
    "300",
    "--step",
    "100",
];

const WIND_CARD: &[&str] = &[
    "wind-card",
    "-v",
    "2700",
    "-m",
    "168",
    "-d",
    "0.308",
    "-b",
    "0.5",
    "--zero-distance",
    "100",
    "--end",
    "300",
];

const COMPARE: &[&str] = &[
    "compare",
    "--load",
    "175 SMK:g7:0.243:175:2650",
    "--load",
    "168 ELD-M:g7:0.523:168:2700",
    "--zero-distance",
    "100",
    "--end",
    "300",
    "--step",
    "100",
];

#[test]
fn come_ups_table_is_byte_identical() {
    let got = stdout_of(COME_UPS);
    let want = include_str!("golden/come_ups_table.txt");
    assert_eq!(got, want);
}

#[test]
fn come_ups_csv_is_byte_identical() {
    let got = stdout_of(&[COME_UPS, &["-o", "csv"]].concat());
    let want = include_str!("golden/come_ups_csv.txt");
    assert_eq!(got, want);
}

#[test]
fn come_ups_json_is_byte_identical() {
    let got = stdout_of(&[COME_UPS, &["-o", "json"]].concat());
    let want = include_str!("golden/come_ups_json.txt");
    assert_eq!(got, want);
}

#[test]
fn range_table_table_is_byte_identical() {
    let got = stdout_of(RANGE_TABLE);
    let want = include_str!("golden/range_table_table.txt");
    assert_eq!(got, want);
}

#[test]
fn range_table_csv_is_byte_identical() {
    let got = stdout_of(&[RANGE_TABLE, &["-o", "csv"]].concat());
    let want = include_str!("golden/range_table_csv.txt");
    assert_eq!(got, want);
}

#[test]
fn range_table_json_is_byte_identical() {
    let got = stdout_of(&[RANGE_TABLE, &["-o", "json"]].concat());
    let want = include_str!("golden/range_table_json.txt");
    assert_eq!(got, want);
}

#[test]
fn wind_card_table_is_byte_identical() {
    let got = stdout_of(WIND_CARD);
    let want = include_str!("golden/wind_card_table.txt");
    assert_eq!(got, want);
}

#[test]
fn wind_card_csv_is_byte_identical() {
    let got = stdout_of(&[WIND_CARD, &["-o", "csv"]].concat());
    let want = include_str!("golden/wind_card_csv.txt");
    assert_eq!(got, want);
}

#[test]
fn wind_card_json_is_byte_identical() {
    let got = stdout_of(&[WIND_CARD, &["-o", "json"]].concat());
    let want = include_str!("golden/wind_card_json.txt");
    assert_eq!(got, want);
}

#[test]
fn compare_table_is_byte_identical() {
    let got = stdout_of(COMPARE);
    let want = include_str!("golden/compare_table.txt");
    assert_eq!(got, want);
}

#[test]
fn compare_csv_is_byte_identical() {
    let got = stdout_of(&[COMPARE, &["-o", "csv"]].concat());
    let want = include_str!("golden/compare_csv.txt");
    assert_eq!(got, want);
}

#[test]
fn compare_json_is_byte_identical() {
    let got = stdout_of(&[COMPARE, &["-o", "json"]].concat());
    let want = include_str!("golden/compare_json.txt");
    assert_eq!(got, want);
}
