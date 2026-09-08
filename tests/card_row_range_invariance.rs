//! MBA-1476: the value a card reports at a range must not depend on the card's `--start`,
//! `--end` or `--step`.
//!
//! That invariant is the whole contract of a DOPE card — `700 yd` means the same drop, drift,
//! velocity and time whichever card the shooter printed it from — and it was broken in shipped
//! code on all four card surfaces plus their bridge equivalents. Each surface sampled its
//! trajectory on a grid whose spacing WAS the requested `--step`, anchored at zero, then picked
//! the sample nearest each requested range and accepted it if it lay within one and a half
//! whole steps. So:
//!
//! * `--start 300 --step 200` put every requested range exactly between two samples, and the
//!   `300 yd` row reported the 200-yard (or 400-yard) numbers;
//! * extending `--end` moved the grid's far end and silently rewrote rows that had not moved;
//! * a range past the end of a load's flight was answered with the last sample, so the final
//!   two rows of a card could be byte-identical 100 or 200 yards apart.
//!
//! Nothing in the output said any of this had happened. There was no test for the invariant
//! itself, which is why it shipped — so the first test here is that property, stated directly,
//! and the rest are the exact reproductions plus the refusal that replaces the substitution.
//!
//! Sampling is now on the engine's own `CARD_SAMPLE_INTERVAL_M` grid (independent of the rows
//! asked for) and the value at a range is interpolated at exactly that range; a range the
//! solved flight does not span is an error naming it.

use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_ballistics")
}

/// The .308-class load the MBA-1476 reproduction was filed with.
const LOAD: &[&str] = &[
    "-v", "2700", "-b", "0.243", "-m", "175", "-d", "0.308", "--drag-model", "g7",
    "--zero-distance", "100", "--sight-height", "1.5", "--temperature", "59", "--pressure",
    "29.92", "--humidity", "50", "--altitude", "0",
];

fn run(args: &[&str]) -> std::process::Output {
    Command::new(bin()).args(args).output().expect("run ballistics")
}

fn stdout_of(args: &[&str]) -> String {
    let out = run(args);
    assert!(
        out.status.success(),
        "`{}` failed: {}",
        args.join(" "),
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).expect("utf8")
}

/// The CSV data lines (header dropped) of a card invocation.
fn csv_rows(args: &[&str]) -> Vec<String> {
    stdout_of(args)
        .lines()
        .skip(1)
        .filter(|l| !l.trim().is_empty())
        .map(str::to_owned)
        .collect()
}

/// The one CSV row whose leading range column is `range`, panicking when the card does not
/// contain it (a card that dropped the row is as wrong as one that faked it).
fn row_at(rows: &[String], range: u32) -> String {
    let key = format!("{range},");
    rows.iter()
        .find(|l| l.starts_with(&key))
        .unwrap_or_else(|| panic!("no {range} yd row among:\n{}", rows.join("\n")))
        .clone()
}

/// One numeric field of a CSV row.
fn field(row: &str, index: usize) -> f64 {
    row.split(',')
        .nth(index)
        .unwrap_or_else(|| panic!("row {row} has no field {index}"))
        .parse()
        .unwrap_or_else(|_| panic!("field {index} of {row} is not a number"))
}

fn range_table(start: u32, end: u32, step: u32) -> Vec<String> {
    let (start, end, step) = (start.to_string(), end.to_string(), step.to_string());
    let mut args: Vec<&str> = vec!["range-table"];
    args.extend_from_slice(LOAD);
    args.extend_from_slice(&[
        "--adjustment-unit", "mil", "--start", &start, "--end", &end, "--step", &step, "--output",
        "csv",
    ]);
    csv_rows(&args)
}

// ---------------------------------------------------------------------------------------
// 1. The property itself: the row at a range is the same row on every card containing it.
// ---------------------------------------------------------------------------------------

/// THE test whose absence let MBA-1476 ship. `700 yd` is solved through six different
/// start/end/step combinations — two of them with a `--start` that is not a multiple of
/// `--step`, which is the case that was broken — and every one must produce the identical row.
#[test]
fn range_table_row_is_invariant_to_start_end_and_step() {
    // (start, end, step), each chosen so that 700 is on the row grid.
    let combos = [
        (100, 1300, 100),
        (300, 1300, 200), // start not a multiple of step: every row was offset
        (300, 1300, 400), // ditto, and a much coarser step
        (500, 700, 100),  // 700 is the LAST row: the far-end collapse case
        (700, 700, 100),  // a one-row card
        (100, 900, 300),  // a third grid alignment entirely
    ];

    let mut reference: Option<(String, (u32, u32, u32))> = None;
    for combo in combos {
        let (start, end, step) = combo;
        let row = row_at(&range_table(start, end, step), 700);
        match &reference {
            None => reference = Some((row, combo)),
            Some((want, from)) => assert_eq!(
                &row, want,
                "the 700 yd row differs between --start {} --end {} --step {} and \
                 --start {start} --end {end} --step {step}",
                from.0, from.1, from.2
            ),
        }
    }
}

/// The same property on `come-ups`. Its Come-up column is a difference between consecutive
/// rows and so legitimately depends on the step; the dial, velocity, energy and time columns
/// do not, and those are what this pins.
#[test]
fn come_ups_row_is_invariant_to_start_and_step() {
    let come_ups = |start: u32, end: u32, step: u32| -> Vec<String> {
        let (start, end, step) = (start.to_string(), end.to_string(), step.to_string());
        let mut args: Vec<&str> = vec!["come-ups"];
        args.extend_from_slice(LOAD);
        args.extend_from_slice(&[
            "--adjustment-unit", "mil", "--start", &start, "--end", &end, "--step", &step, "-o",
            "csv",
        ]);
        csv_rows(&args)
    };

    // range_yd,drop_mil,come_up_mil,velocity_fps,energy_ft-lb,time_s — every field but
    // come_up (index 2) is a property of the range alone.
    let independent = |row: &str| -> Vec<f64> {
        [0usize, 1, 3, 4, 5].iter().map(|&i| field(row, i)).collect()
    };

    let reference = independent(&row_at(&come_ups(100, 1000, 100), 700));
    for (start, end, step) in [(300, 1300, 200), (300, 1100, 400), (700, 700, 100)] {
        assert_eq!(
            independent(&row_at(&come_ups(start, end, step), 700)),
            reference,
            "come-ups 700 yd row changed under --start {start} --end {end} --step {step}"
        );
    }
}

/// The same property on `wind-card`, whose cells are a drift matrix rather than a row of
/// scalars — the nearest-sample lookup ran once per (range, wind speed) cell.
#[test]
fn wind_card_row_is_invariant_to_start_and_step() {
    let wind_card = |start: u32, end: u32, step: u32| -> Vec<String> {
        let (start, end, step) = (start.to_string(), end.to_string(), step.to_string());
        let mut args: Vec<&str> = vec!["wind-card"];
        args.extend_from_slice(LOAD);
        args.extend_from_slice(&[
            "--adjustment-unit", "mil", "--wind-speeds", "5,10,15", "--start", &start, "--end",
            &end, "--step", &step, "-o", "csv",
        ]);
        csv_rows(&args)
    };

    let reference = row_at(&wind_card(100, 1000, 100), 700);
    for (start, end, step) in [(300, 1300, 200), (300, 1100, 400), (700, 700, 100)] {
        assert_eq!(
            row_at(&wind_card(start, end, step), 700),
            reference,
            "wind-card 700 yd row changed under --start {start} --end {end} --step {step}"
        );
    }
}

/// The same property on `compare`, the fourth surface — which had no tolerance check at all,
/// so its nearest-sample pick was accepted however far away it was.
#[test]
fn compare_row_is_invariant_to_start_and_step() {
    let compare = |start: u32, end: u32, step: u32| -> Vec<String> {
        let (start, end, step) = (start.to_string(), end.to_string(), step.to_string());
        csv_rows(&[
            "compare",
            "--load",
            "A:g7:0.243:175:2700",
            "--load",
            "B:g7:0.523:168:2700",
            "--zero-distance",
            "100",
            "--start",
            &start,
            "--end",
            &end,
            "--step",
            &step,
            "-o",
            "csv",
        ])
    };

    // `compare` requires --start < --end, so its single-row case is a two-row card whose
    // first row is the one under test.
    let reference = row_at(&compare(100, 1000, 100), 700);
    for (start, end, step) in [(300, 1300, 200), (300, 1100, 400), (700, 800, 100)] {
        assert_eq!(
            row_at(&compare(start, end, step), 700),
            reference,
            "compare 700 yd row changed under --start {start} --end {end} --step {step}"
        );
    }
}

// ---------------------------------------------------------------------------------------
// 2. The reproduction from the report, with its real numbers.
// ---------------------------------------------------------------------------------------

/// The filed reproduction. On `--start 300 --end 1300`, the 300 yd row read
/// 2162 fps at `--step 100` (correct), 2334 fps at `--step 200` (the 200-yard value) and
/// 1999 fps at `--step 400` (the 400-yard value).
#[test]
fn mba_1476_the_300_yard_row_no_longer_moves_with_step() {
    // Velocity is column 5 of the range-table CSV.
    const VELOCITY: usize = 5;
    let truth = field(&row_at(&range_table(100, 1300, 100), 300), VELOCITY);
    assert!(
        (truth - 2162.0).abs() <= 1.0,
        "the 300 yd velocity on an aligned card should be ~2162 fps, got {truth}"
    );

    for step in [200u32, 400] {
        let got = field(&row_at(&range_table(300, 1300, step), 300), VELOCITY);
        assert!(
            (got - truth).abs() <= 1.0,
            "--start 300 --step {step} reports {got} fps at 300 yd, not {truth}"
        );
    }

    // And the two specific wrong answers are gone, not merely close.
    for (step, wrong, whose) in [(200u32, 2334.0, "200 yd"), (400, 1999.0, "400 yd")] {
        let got = field(&row_at(&range_table(300, 1300, step), 300), VELOCITY);
        assert!(
            (got - wrong).abs() > 50.0,
            "--start 300 --step {step} still reports the {whose} velocity ({got} fps) \
             against the 300 yd row"
        );
    }
}

/// The second half of the reproduction: with `--start 300 --step 200`, extending `--end`
/// from 700 to 900 rewrote the 700 yd row (4.4299 mil / 95.69 in / 1690 fps became
/// 7.2378 mil / 208.45 in / 1404 fps) and made the last two rows byte-identical 200 yards
/// apart.
#[test]
fn mba_1476_extending_end_does_not_rewrite_an_existing_row() {
    let short = range_table(300, 700, 200);
    let long = range_table(300, 900, 200);

    assert_eq!(
        row_at(&short, 700),
        row_at(&long, 700),
        "extending --end from 700 to 900 rewrote the 700 yd row"
    );

    // The trailing-row collapse: 700 and 900 were byte-identical.
    assert_ne!(
        row_at(&long, 700),
        row_at(&long, 900),
        "the 700 and 900 yd rows are identical 200 yards apart"
    );

    // Both rows are also each other's genuine neighbours, not copies of the far end:
    // drop grows monotonically and velocity falls.
    let (drop700, drop900) = (field(&row_at(&long, 700), 1), field(&row_at(&long, 900), 1));
    let (vel700, vel900) = (field(&row_at(&long, 700), 5), field(&row_at(&long, 900), 5));
    assert!(drop900 > drop700, "drop must grow with range ({drop700} -> {drop900})");
    assert!(vel900 < vel700, "velocity must fall with range ({vel700} -> {vel900})");
}

// ---------------------------------------------------------------------------------------
// 3. A row the flight cannot supply is an error, never a substitution.
// ---------------------------------------------------------------------------------------

/// A light, low-BC bullet whose solved flight ends around 871 yd. Asked for rows out to
/// 2000 yd, the shipped code answered 900 yd with the 800 yd sample — the two rows printed
/// byte-identical — and then silently truncated the card. Nothing on screen said the 900 yd
/// line was a copy.
const UNREACHABLE_LOAD: &[&str] = &[
    "-v", "900", "-b", "0.1", "-m", "40", "-d", "0.224", "--drag-model", "g1", "--zero-distance",
    "100",
];

#[test]
fn a_range_the_flight_cannot_reach_is_refused_not_substituted() {
    let mut args: Vec<&str> = vec!["range-table"];
    args.extend_from_slice(UNREACHABLE_LOAD);
    args.extend_from_slice(&["--start", "100", "--end", "2000", "--step", "100", "-o", "csv"]);
    let out = run(&args);

    assert!(
        !out.status.success(),
        "a card asking for ranges past the flight must fail, but it printed:\n{}",
        String::from_utf8_lossy(&out.stdout)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("no trajectory sample at 900 yd"),
        "the error must name the range it could not supply; got: {stderr}"
    );
    // The unreachable row is refused, not printed with borrowed numbers.
    assert!(
        !String::from_utf8_lossy(&out.stdout).contains("\n900,"),
        "a 900 yd row was printed for a flight that does not reach it"
    );
}

/// The same refusal on the wind card, whose unsatisfiable cells used to become `0.0` —
/// a drift of zero is a perfectly plausible-looking number, and it was printed for a range
/// the bullet never reached.
#[test]
fn wind_card_refuses_a_range_the_flight_cannot_reach() {
    let mut args: Vec<&str> = vec!["wind-card"];
    args.extend_from_slice(UNREACHABLE_LOAD);
    args.extend_from_slice(&[
        "--wind-speeds", "10", "--start", "100", "--end", "2000", "--step", "100", "-o", "csv",
    ]);
    let out = run(&args);

    assert!(
        !out.status.success(),
        "the wind card must fail rather than print zero-drift cells, but it printed:\n{}",
        String::from_utf8_lossy(&out.stdout)
    );
    assert!(
        String::from_utf8_lossy(&out.stderr).contains("no trajectory sample at"),
        "the error must name the range it could not supply"
    );
}

// ---------------------------------------------------------------------------------------
// 4. The same three properties on the bridge's card services, which mobile calls directly.
// ---------------------------------------------------------------------------------------

mod bridge {
    use ballistics_engine::card_service::{
        come_ups_v1, range_table_v1, wind_card_v1, CardRequestV1, CardResponseV1,
        CardServiceError,
    };

    fn request(start: f64, end: f64, step: f64) -> CardRequestV1 {
        serde_json::from_value(serde_json::json!({
            "units": "imperial",
            "muzzle_velocity": 2700.0,
            "ballistic_coefficient": 0.243,
            "mass": 175.0,
            "diameter": 0.308,
            "drag_model": "g7",
            "sight_height": 1.5,
            "zero_distance": 100.0,
            "temperature": 59.0,
            "pressure": 29.92,
            "humidity": 50.0,
            "altitude": 0.0,
            "wind_speed": 10.0,
            "wind_direction_deg": 90.0,
            "start": start,
            "end": end,
            "step": step,
            "adjustment_unit": "mil",
            "wind_speeds": [5.0, 10.0, 15.0],
        }))
        .expect("card request")
    }

    fn row_at(card: &CardResponseV1, range: f64) -> String {
        let row = card
            .rows
            .iter()
            .find(|r| (r.range - range).abs() < 1e-6)
            .unwrap_or_else(|| panic!("no {range} row in the card"));
        // Serialized rather than field-by-field so a future column is covered automatically.
        serde_json::to_string(row).expect("row json")
    }

    /// `card.range_table`, `card.come_ups` and `card.wind` all read their rows through the
    /// same sampling path, so the invariant is asserted on all three.
    #[test]
    fn card_services_report_the_same_row_for_a_range_whatever_the_card_asked_for() {
        for (name, service) in [
            ("range_table", range_table_v1 as fn(&CardRequestV1) -> _),
            ("come_ups", come_ups_v1),
            ("wind", wind_card_v1),
        ] {
            let reference = row_at(&service(&request(100.0, 1000.0, 100.0)).unwrap(), 700.0);
            for (start, end, step) in [
                (300.0, 1300.0, 200.0),
                (300.0, 1100.0, 400.0),
                (700.0, 700.0, 100.0),
            ] {
                let got = row_at(&service(&request(start, end, step)).unwrap(), 700.0);
                // come-ups' Come-up column is a row-to-row difference and legitimately
                // depends on the step; every other field is a property of the range.
                let (got, reference) = if name == "come_ups" {
                    (strip_come_up(&got), strip_come_up(&reference))
                } else {
                    (got.clone(), reference.clone())
                };
                assert_eq!(
                    got, reference,
                    "card.{name}'s 700 yd row changed under start {start} / end {end} / \
                     step {step}"
                );
            }
        }
    }

    fn strip_come_up(row_json: &str) -> String {
        let mut value: serde_json::Value = serde_json::from_str(row_json).expect("row json");
        value
            .as_object_mut()
            .expect("row object")
            .remove("come_up");
        value.to_string()
    }

    /// The bridge refuses an unsatisfiable row too, with a message naming the range —
    /// rather than handing a mobile caller a row built from a different range's numbers.
    #[test]
    fn card_services_refuse_a_range_the_flight_cannot_reach() {
        let unreachable: CardRequestV1 = serde_json::from_value(serde_json::json!({
            "units": "imperial",
            "muzzle_velocity": 900.0,
            "ballistic_coefficient": 0.1,
            "mass": 40.0,
            "diameter": 0.224,
            "drag_model": "g1",
            "zero_distance": 100.0,
            "start": 100.0,
            "end": 2000.0,
            "step": 100.0,
            "wind_speeds": [10.0],
        }))
        .expect("card request");

        for (name, service) in [
            ("range_table", range_table_v1 as fn(&CardRequestV1) -> _),
            ("come_ups", come_ups_v1),
            ("wind", wind_card_v1),
        ] {
            match service(&unreachable) {
                Err(CardServiceError::Trajectory(message)) => assert!(
                    message.contains("no trajectory sample at"),
                    "card.{name}: the error must name the range it could not supply, \
                     got: {message}"
                ),
                Err(other) => panic!("card.{name}: unexpected error {other}"),
                Ok(card) => panic!(
                    "card.{name}: returned {} rows for a flight that does not reach 2000 yd",
                    card.rows.len()
                ),
            }
        }
    }
}
