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
//! and the rest are the exact reproductions plus the truncation that replaces the substitution.
//!
//! Sampling is now on the engine's own `CARD_SAMPLE_INTERVAL_M` grid (independent of the rows
//! asked for) and the value at a range is interpolated at exactly that range. A range the
//! solved flight does not span is never fabricated and never substituted: the card ends at
//! the last row the flight reaches and says what it left out, on stderr for a CLI reader and
//! as a structured field on the JSON and bridge surfaces. Only a card that reaches NO row at
//! all is an error — there is nothing to print.

use std::process::Command;

#[path = "support/pdf_text.rs"]
mod pdf_text;

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
// 3. A row the flight cannot supply is never a substitution -- and never the whole card.
// ---------------------------------------------------------------------------------------

/// A light, low-BC bullet whose solved flight ends around 872 yd. Asked for rows out to
/// 2000 yd, the shipped code answered 900 yd with the 800 yd sample — the two rows printed
/// byte-identical — and then silently truncated the card. Nothing on screen said the 900 yd
/// line was a copy.
const UNREACHABLE_LOAD: &[&str] = &[
    "-v", "900", "-b", "0.1", "-m", "40", "-d", "0.224", "--drag-model", "g1", "--zero-distance",
    "100",
];

// ---------------------------------------------------------------------------------------
// 3a. All FIVE card surfaces, as one table (MBA-1478).
//
// Five surfaces share one sampling concern, and the way MBA-1476 shipped was that they had
// drifted into seven call sites with three different tolerance behaviours. `adaptive-card`
// then sat out that fix entirely -- its output was byte-identical before and after -- and
// kept all three of the behaviours the other four shed: it hard-REFUSED a card its load
// outran, it carried its own private `end * 1.02` sampling headroom (an even wider version of
// the fractional-`--end` hole `card_sample_max_range_m` closed), and it quoted a reach that
// moved with `--end` because it read the last grid sample rather than the flight's own
// terminal distance.
//
// So the properties below are asserted against every surface at once, from this table. A test
// that exercised one surface would still pass if another quietly stopped using the shared
// path -- which is exactly what happened -- and that is the failure this table exists to make
// impossible. Adding a sixth card surface means adding it here.
// ---------------------------------------------------------------------------------------

/// The two loads these tests use, in whichever spelling a surface's own flags want.
#[derive(Clone, Copy)]
enum Load {
    /// The .308-class reference load: reaches far past any card asked for below.
    Reference,
    /// [`UNREACHABLE_LOAD`]: a flight that ends around 872 yd, so an ordinary card outruns it.
    Short,
}

impl Load {
    /// The flags the four single-load surfaces take.
    fn flags(self) -> &'static [&'static str] {
        match self {
            Self::Reference => LOAD,
            Self::Short => UNREACHABLE_LOAD,
        }
    }

    /// `compare`'s own `name:drag:bc:mass:velocity` spelling. It needs two loads to compare,
    /// and shares ONE range axis across them, so the pair is chosen to make the axis end for
    /// the intended reason: two long loads for `Reference`, and for `Short` a long one paired
    /// with the .22-class flight, which is then the load that sets the shared end.
    fn compare_flags(self) -> &'static [&'static str] {
        match self {
            Self::Reference => &["--load", "far:g7:0.3:175:2800", "--load", "also:g7:0.243:175:2700"],
            Self::Short => &["--load", "far:g7:0.3:175:2800", "--load", "short:g1:0.1:40:900"],
        }
    }
}

/// One CLI card surface, plus how to ask it for a card.
struct Surface {
    name: &'static str,
    /// The whole argv for a card over `start..=end` in `output` format.
    argv: fn(Load, &str, &str, &str) -> Vec<String>,
    /// The whole argv for a card whose FURTHEST REQUESTED ROW is exactly `end`.
    ///
    /// Spelled per surface because the commands do not accept the same degenerate card:
    /// three of them take a one-row `--start end --end end`, `compare` requires a strictly
    /// forward interval, and `adaptive-card` has no `--step` at all. What they have in common
    /// is the only thing this is used for — a card whose last row sits exactly on `--end`,
    /// which is the row the sampled span has to clear.
    furthest_row_argv: fn(Load, f64, &str) -> Vec<String>,
    /// True when this surface takes [`Load::flags`], i.e. solves the very same flight as the
    /// others that do. `compare` does not: it names its loads in its own
    /// `name:drag:bc:mass:velocity` spelling, which carries no diameter.
    ///
    /// Used to hold every surface that solves ONE flight to ONE reach figure, to the yard.
    /// That equality is what separates the flight's own terminal distance from the last point
    /// of a sampling grid — the two differ by less than a grid cell, so a band would not tell
    /// them apart, and the last grid point is the wrong figure for the reason
    /// `SampledFlight::reach_m` documents.
    solves_the_shared_load: bool,
}

fn owned(parts: &[&str]) -> Vec<String> {
    parts.iter().map(|s| (*s).to_string()).collect()
}

/// Every card surface that samples a trajectory for a range card.
fn surfaces() -> Vec<Surface> {
    /// A stepped card on one of the three single-load, `--step`-taking surfaces.
    /// `wind-card` additionally needs a wind column to have a matrix at all.
    fn stepped(command: &str, load: Load, start: &str, end: &str, step: &str, output: &str)
        -> Vec<String>
    {
        let mut argv = vec![command.to_string()];
        argv.extend(owned(load.flags()));
        if command == "wind-card" {
            argv.extend(owned(&["--wind-speeds", "10"]));
        }
        argv.extend(owned(&["--start", start, "--end", end, "--step", step, "-o", output]));
        argv
    }
    /// A one-row card whose only row is `end`, for the surfaces that accept `start == end`.
    fn one_row(command: &str, load: Load, end: f64, output: &str) -> Vec<String> {
        let end = format!("{end}");
        stepped(command, load, &end, &end, "1", output)
    }
    vec![
        Surface {
            name: "range-table",
            argv: |load, start, end, output| {
                stepped("range-table", load, start, end, "100", output)
            },
            furthest_row_argv: |load, end, output| one_row("range-table", load, end, output),
            solves_the_shared_load: true,
        },
        Surface {
            name: "come-ups",
            argv: |load, start, end, output| {
                stepped("come-ups", load, start, end, "100", output)
            },
            furthest_row_argv: |load, end, output| one_row("come-ups", load, end, output),
            solves_the_shared_load: true,
        },
        Surface {
            name: "wind-card",
            argv: |load, start, end, output| {
                stepped("wind-card", load, start, end, "100", output)
            },
            furthest_row_argv: |load, end, output| one_row("wind-card", load, end, output),
            solves_the_shared_load: true,
        },
        Surface {
            name: "compare",
            argv: |load, start, end, output| {
                let mut argv = vec!["compare".to_string()];
                argv.extend(owned(load.compare_flags()));
                argv.extend(owned(&["--zero-distance", "100"]));
                argv.extend(owned(&[
                    "--start", start, "--end", end, "--step", "100", "-o", output,
                ]));
                argv
            },
            // `compare` refuses `--start == --end`, so its last row is put on `--end` with a
            // half-`end` step from a half-`end` start instead.
            furthest_row_argv: |load, end, output| {
                let (start, step, end) =
                    (format!("{:.4}", end / 2.0), format!("{:.4}", end / 2.0), format!("{end}"));
                let mut argv = vec!["compare".to_string()];
                argv.extend(owned(load.compare_flags()));
                argv.extend(owned(&["--zero-distance", "100"]));
                argv.extend(owned(&[
                    "--start", &start, "--end", &end, "--step", &step, "-o", output,
                ]));
                argv
            },
            solves_the_shared_load: false,
        },
        Surface {
            name: "adaptive-card",
            argv: |load, start, end, output| {
                let mut argv = vec!["adaptive-card".to_string()];
                argv.extend(owned(load.flags()));
                argv.extend(owned(&["--start", start, "--end", end, "-o", output]));
                argv
            },
            // No `--step`: the domain end is always a seeded row, so half the requested end
            // is a start that stays forward and positive however small `end` gets.
            furthest_row_argv: |load, end, output| {
                let (start, end) = (format!("{:.4}", end / 2.0), format!("{end}"));
                let mut argv = vec!["adaptive-card".to_string()];
                argv.extend(owned(load.flags()));
                argv.extend(owned(&["--start", &start, "--end", &end, "-o", output]));
                argv
            },
            solves_the_shared_load: true,
        },
    ]
}

fn run_argv(argv: &[String]) -> std::process::Output {
    let args: Vec<&str> = argv.iter().map(String::as_str).collect();
    run(&args)
}

/// The leading range column of every CSV data row, in order. Parsed as `f64` because
/// `adaptive-card` chooses its own rows and they need not be whole yards.
fn printed_ranges(stdout: &str) -> Vec<f64> {
    stdout
        .lines()
        .skip(1)
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            l.split(',')
                .next()
                .and_then(|c| c.trim().parse::<f64>().ok())
                .unwrap_or_else(|| panic!("row `{l}` does not start with a range"))
        })
        .collect()
}

/// The `reaches only N unit` figure out of a truncation warning or a refusal, normalised.
///
/// A warning continues past it with `; the rows past that ...` and a refusal ends there inside
/// the `Error: "..."` quoting, so the figure is taken as its own two words rather than as
/// whatever happens to follow it.
fn reach_phrase(stderr: &str) -> String {
    let (_, tail) = stderr
        .split_once("reaches only ")
        .unwrap_or_else(|| panic!("no reach stated in: {stderr}"));
    let mut words = tail.split_whitespace();
    let value = words.next().unwrap_or_else(|| panic!("no reach value in: {stderr}"));
    let unit: String = words
        .next()
        .unwrap_or_else(|| panic!("no reach unit in: {stderr}"))
        .chars()
        .take_while(char::is_ascii_alphabetic)
        .collect();
    format!("{value} {unit}")
}

/// The same figure as a number of yards.
fn reach_yards(stderr: &str) -> f64 {
    let phrase = reach_phrase(stderr);
    let (value, unit) = phrase.split_once(' ').unwrap_or_else(|| panic!("bad reach: {phrase}"));
    assert_eq!(unit, "yd", "these cards are imperial");
    value.parse().unwrap_or_else(|_| panic!("bad reach: {phrase}"))
}

/// Where [`Load::Short`]'s flight ends, in yards.
///
/// A band, because `compare` names its loads in its own `name:drag:bc:mass:velocity`
/// spelling, which carries no diameter and so solves a flight a yard away from
/// `UNREACHABLE_LOAD`'s. It is a sanity bound only — the sharp assertion is
/// `solves_the_shared_load`'s exact equality between the four surfaces that DO solve the
/// identical flight.
const SHORT_LOAD_REACH_YD: std::ops::RangeInclusive<f64> = 865.0..=875.0;

/// EVERY card surface prints the rows its load reaches and says what it left out.
///
/// `adaptive-card` was the one that did not: `.22 LR --end 1200` exited 1 with no rows at all,
/// where `come-ups` on the same load printed eighteen and a notice. It went through
/// `HoldCurve`, so its rows were never read from the wrong range (MBA-1476's own defect never
/// touched it) -- it simply refused the whole card, which on a default invocation reads as a
/// broken tool for exactly the ordinary loads a range card is most useful for.
///
/// The notice's `last_row` is asserted against the last row actually PRINTED, which is a
/// property every surface can be held to whatever rows it chose.
#[test]
fn every_card_surface_truncates_rather_than_refusing_a_card_its_load_outruns() {
    let mut shared_load_reach: Option<(&str, f64)> = None;
    for surface in surfaces() {
        let name = surface.name;
        let out = run_argv(&(surface.argv)(Load::Short, "100", "2000", "csv"));
        assert!(
            out.status.success(),
            "{name} must print the rows it reaches, not refuse the card: {}",
            String::from_utf8_lossy(&out.stderr)
        );

        let stdout = String::from_utf8_lossy(&out.stdout);
        let ranges = printed_ranges(&stdout);
        assert!(
            ranges.len() >= 2,
            "{name} printed no usable card:\n{stdout}"
        );
        let last_row = *ranges.last().expect("a last row");

        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains(&format!("card truncated at {last_row:.0} yd")),
            "{name} must name the last row it printed ({last_row}); got: {stderr}"
        );
        assert!(
            stderr.contains("2000 yd was requested"),
            "{name} must name the end that was requested; got: {stderr}"
        );
        // A property of the LOAD, and the last row printed is inside it: a card stops at the
        // furthest row its flight reaches, never past it.
        let reach = reach_yards(&stderr);
        assert!(
            SHORT_LOAD_REACH_YD.contains(&reach),
            "{name} must quote this load's own reach, got {reach} yd"
        );
        assert!(
            last_row <= reach,
            "{name} printed a {last_row} yd row for a flight that stops at {reach} yd:\n{stdout}"
        );

        // ONE flight, ONE reach. Every surface here solves the identical load, so any
        // disagreement between them means one of them is quoting something other than the
        // flight's terminal distance -- which is precisely what reading the figure off a
        // sampling grid does: the last grid point sits up to a whole cell short of the end
        // of the flight, so it is a different, smaller, and quietly wrong number.
        if surface.solves_the_shared_load {
            match shared_load_reach {
                None => shared_load_reach = Some((name, reach)),
                Some((first, want)) => assert_eq!(
                    reach, want,
                    "{name} and {first} solve the same load and must quote the same reach"
                ),
            }
        }
    }
    assert!(
        shared_load_reach.is_some(),
        "the shared-load reach must have been checked at least once"
    );
}

/// `adaptive-card`'s machine-readable half, spelled exactly as the other four surfaces spell
/// it (MBA-1478).
///
/// Its `-o json` is `AdaptiveCardReportV1` serialized verbatim, so the truncation is the one
/// thing added to it: a caller that already reads `truncated` off a range table reads this
/// one the same way, rather than learning a fifth spelling for the same three numbers.
#[test]
fn the_adaptive_cards_json_carries_the_same_truncated_block_as_the_others() {
    let mut args: Vec<&str> = vec!["adaptive-card"];
    args.extend_from_slice(UNREACHABLE_LOAD);
    args.extend_from_slice(&["--start", "100", "--end", "2000", "-o", "json"]);
    let report: serde_json::Value = serde_json::from_str(&stdout_of(&args)).expect("json");

    let truncated = report
        .get("truncated")
        .expect("a truncated adaptive card must say so structurally");
    assert_eq!(truncated["requested_end"], 2000.0);
    let last_row = truncated["last_row"].as_f64().expect("last_row");
    assert!(
        (truncated["reach"].as_f64().expect("reach") - 871.7).abs() < 1.0,
        "reach must be the flight's own terminal distance; got {truncated}"
    );
    assert!(
        truncated["message"].as_str().expect("message").starts_with("card truncated at"),
        "the message must be the shared wording; got {truncated}"
    );

    // The rows the report DID produce all sit inside the truncated domain. `rows[].range`
    // is metres on this surface (`AdaptiveCardReportV1` is metric throughout), while the
    // truncation block is in the card's display unit, as it is everywhere else.
    let rows = report["rows"].as_array().expect("rows");
    assert!(!rows.is_empty(), "the reachable rows must be present");
    for row in rows {
        let range_m = row["range"].as_f64().expect("range");
        assert!(
            range_m / 0.9144 <= last_row + 1e-6,
            "a {range_m} m row is past the truncated end of {last_row} yd"
        );
    }

    // Additive: a card that runs to an end it can reach says nothing at all.
    let mut whole: Vec<&str> = vec!["adaptive-card"];
    whole.extend_from_slice(UNREACHABLE_LOAD);
    whole.extend_from_slice(&["--start", "100", "--end", "800", "-o", "json"]);
    let whole: serde_json::Value = serde_json::from_str(&stdout_of(&whole)).expect("json");
    assert!(
        whole.get("truncated").is_none(),
        "an untruncated adaptive card must be byte-identical to before: {whole}"
    );
}

/// A card whose flight reaches NONE of its rows has nothing to print, so every surface keeps
/// refusing it -- in the one shared wording, naming the first row it could not supply.
#[test]
fn every_card_surface_still_refuses_a_card_with_no_reachable_row() {
    for surface in surfaces() {
        let name = surface.name;
        let out = run_argv(&(surface.argv)(Load::Short, "1000", "2000", "csv"));
        assert!(
            !out.status.success(),
            "{name} has no card to print here, so it must fail:\n{}",
            String::from_utf8_lossy(&out.stdout)
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains("no trajectory sample at 1000 yd"),
            "{name} must name the range it could not supply; got: {stderr}"
        );
        let reach = reach_yards(&stderr);
        assert!(
            SHORT_LOAD_REACH_YD.contains(&reach),
            "{name} must state the load's reach in its refusal, got {reach} yd"
        );
    }
}

/// The rows a card CAN supply are still the shooter's card.
///
/// This test deliberately replaces `a_range_the_flight_cannot_reach_is_refused_not_substituted`,
/// which asserted that a card containing any unreachable row failed outright with no output.
/// Not fabricating the unreachable row was right and is still asserted here; refusing the
/// reachable ones with it was an overcorrection, and on a no-flags invocation — `range-table`
/// defaults `--end` to 1200 yd, which ordinary .22 LR / 9 mm / .45 ACP loads do not reach —
/// it made the command look broken. The card now truncates: every row the flight reaches,
/// none it does not, and a warning naming the requested end and the load's actual reach.
#[test]
fn a_card_prints_the_rows_it_reaches_and_says_what_it_truncated() {
    let mut args: Vec<&str> = vec!["range-table"];
    args.extend_from_slice(UNREACHABLE_LOAD);
    args.extend_from_slice(&["--start", "100", "--end", "2000", "--step", "100", "-o", "csv"]);
    let out = run(&args);

    assert!(
        out.status.success(),
        "the reachable rows must still print; got: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);

    // Present: every row inside the flight.
    for range in [100, 200, 300, 400, 500, 600, 700, 800] {
        assert!(
            stdout.contains(&format!("\n{range},")),
            "the {range} yd row is inside this flight and must be printed:\n{stdout}"
        );
    }
    // Absent: the rows past it — left out, never fabricated from a range the bullet does
    // reach, which is the substitution this whole change removed.
    for range in [900, 1000, 1100, 2000] {
        assert!(
            !stdout.contains(&format!("\n{range},")),
            "a {range} yd row was printed for a flight that does not reach it:\n{stdout}"
        );
    }

    // And the truncation is stated, not silent.
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("card truncated at 800 yd")
            && stderr.contains("2000 yd was requested")
            && stderr.contains("reaches only 872 yd"),
        "the warning must name the last row, the requested end and the load's reach; \
         got: {stderr}"
    );

    // The rows it kept are byte-identical to the same card asked for exactly that end:
    // truncating must not perturb the rows that survive it.
    let mut exact: Vec<&str> = vec!["range-table"];
    exact.extend_from_slice(UNREACHABLE_LOAD);
    exact.extend_from_slice(&["--start", "100", "--end", "800", "--step", "100", "-o", "csv"]);
    assert_eq!(
        stdout,
        stdout_of(&exact),
        "a truncated card must equal the card that asked for its last row exactly"
    );
}

/// The machine-readable half of the same notice. A mobile or Flask caller reading JSON cannot
/// see a stderr line, and a card that quietly returns fewer rows than were asked for, with
/// nothing on the response saying so, is the same class of silent failure as the substituted
/// row this work removed.
#[test]
fn a_truncated_card_carries_a_structured_field_on_json() {
    let mut args: Vec<&str> = vec!["range-table"];
    args.extend_from_slice(UNREACHABLE_LOAD);
    args.extend_from_slice(&["--start", "100", "--end", "2000", "--step", "100", "-o", "json"]);
    let card: serde_json::Value =
        serde_json::from_str(&stdout_of(&args)).expect("range-table json");

    let truncated = card
        .get("truncated")
        .expect("a truncated card must carry a `truncated` block");
    assert_eq!(truncated["requested_end"], 2000.0);
    assert_eq!(truncated["last_row"], 800.0);
    assert!(
        (truncated["reach"].as_f64().expect("reach") - 871.7).abs() < 1.0,
        "reach must be the flight's own terminal distance; got {truncated}"
    );
    assert_eq!(
        card["data"].as_array().expect("rows").len(),
        8,
        "only the reachable rows belong in the card"
    );

    // Additive: a card that runs to its requested end says nothing at all.
    let mut whole: Vec<&str> = vec!["range-table"];
    whole.extend_from_slice(UNREACHABLE_LOAD);
    whole.extend_from_slice(&["--start", "100", "--end", "800", "--step", "100", "-o", "json"]);
    let whole: serde_json::Value = serde_json::from_str(&stdout_of(&whole)).expect("json");
    assert!(
        whole.get("truncated").is_none(),
        "an untruncated card must be byte-identical to before: {whole}"
    );
}

/// A card whose flight reaches NONE of its rows has no card to print, so it stays a refusal
/// naming the first range it could not supply.
#[test]
fn a_card_whose_flight_reaches_no_row_at_all_is_still_refused() {
    let mut args: Vec<&str> = vec!["range-table"];
    args.extend_from_slice(UNREACHABLE_LOAD);
    args.extend_from_slice(&["--start", "1000", "--end", "2000", "--step", "100", "-o", "csv"]);
    let out = run(&args);

    assert!(
        !out.status.success(),
        "there is no card to print here, so this must fail:\n{}",
        String::from_utf8_lossy(&out.stdout)
    );
    assert!(
        String::from_utf8_lossy(&out.stderr).contains("no trajectory sample at 1000 yd"),
        "the error must name the range it could not supply; got: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

/// The same truncation on the wind card, whose unsatisfiable cells used to become `0.0` —
/// a drift of zero is a perfectly plausible-looking number, and it was printed for a range
/// the bullet never reached. It is now the end of the matrix, for every wind column at once.
#[test]
fn wind_card_truncates_rather_than_printing_zero_drift_cells() {
    let mut args: Vec<&str> = vec!["wind-card"];
    args.extend_from_slice(UNREACHABLE_LOAD);
    args.extend_from_slice(&[
        "--wind-speeds", "10,20", "--start", "100", "--end", "2000", "--step", "100", "-o", "csv",
    ]);
    let out = run(&args);

    assert!(
        out.status.success(),
        "the reachable rows must still print; got: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("\n800,"),
        "the 800 yd row is inside this flight:\n{stdout}"
    );
    assert!(
        !stdout.contains("\n900,"),
        "a 900 yd drift row was printed for a flight that does not reach it:\n{stdout}"
    );
    // Every printed row still carries one cell per wind speed.
    for row in stdout.lines().skip(1).filter(|l| !l.trim().is_empty()) {
        assert_eq!(
            row.split(',').count(),
            3,
            "row `{row}` is missing a wind column"
        );
    }
    assert!(
        String::from_utf8_lossy(&out.stderr).contains("card truncated at 800 yd"),
        "the truncation must be stated"
    );
}

/// `come-ups` and `compare` truncate on the same terms, and `compare` names the load that
/// ended the range axis all of its loads share.
#[test]
fn come_ups_and_compare_truncate_on_the_same_terms() {
    let mut come_ups: Vec<&str> = vec!["come-ups"];
    come_ups.extend_from_slice(UNREACHABLE_LOAD);
    come_ups.extend_from_slice(&["--start", "100", "--end", "2000", "--step", "100", "-o", "csv"]);
    let out = run(&come_ups);
    assert!(
        out.status.success(),
        "come-ups must print its reachable rows; got: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("\n800,") && !stdout.contains("\n900,"), "{stdout}");
    assert!(
        String::from_utf8_lossy(&out.stderr).contains("card truncated at 800 yd"),
        "come-ups must state the truncation"
    );

    // compare: one load reaches far, one does not. The card runs to the shorter, and says so.
    let out = run(&[
        "compare", "--load", "far:g7:0.3:175:2800", "--load", "short:g1:0.1:40:900",
        "--zero-distance", "100", "--start", "100", "--end", "2000", "--step", "100", "-o", "csv",
    ]);
    assert!(
        out.status.success(),
        "compare must print its reachable rows; got: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("\n800,") && !stdout.contains("\n900,"), "{stdout}");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("load 'short'") && stderr.contains("card truncated at 800 yd"),
        "compare must name the load that ended the shared range axis; got: {stderr}"
    );
}

/// A card's sampled span must clear its last row by at least one whole grid cell — on EVERY
/// surface, because the span is `hold_curve::card_sample_max_range_m` on every surface.
///
/// The span was `end * 1.1`, which leaves less than one ~1-yard grid cell of headroom once
/// `--end` drops below about 10 yd: the last grid point then lands SHORT of the last requested
/// row, and the row is unreachable through no fault of the load. 36 of the 150 fractional ends
/// between 0.1 and 15.0 yd failed that way — 0.1-0.8, 1.1-1.7, 2.1-2.6 and so on — while every
/// integer end (an exact multiple of the grid) passed, which is what made it easy to miss.
///
/// `adaptive-card` had the same hole and a wider one: its own `end * 1.02` leaves less than a
/// whole cell for every `--end` under ~50 yd. Both are gone because both now ask the same
/// helper how far to sample, and this sweep is what says so — it fails the moment a surface
/// goes back to computing its own headroom.
#[test]
fn a_fractional_end_is_still_inside_the_sampled_span_on_every_surface() {
    for surface in surfaces() {
        let name = surface.name;
        for tenths in 1..=150u32 {
            let end = f64::from(tenths) / 10.0;
            let out = run_argv(&(surface.furthest_row_argv)(Load::Reference, end, "csv"));
            assert!(
                out.status.success(),
                "{name}: --end {end} must be inside its own card's sampled span; got: {}",
                String::from_utf8_lossy(&out.stderr)
            );
        }
    }
}

/// The reach a truncated card quotes is a property of the LOAD, not of the request — on
/// EVERY surface, because every surface reads it off `hold_curve::SampledFlight::reach_m`.
///
/// It used to be the last SAMPLE distance, and the sampled span is derived from `--end`, so
/// the same load's "reaches only ..." figure moved when only `--end` had. It is now the solved
/// flight's own terminal distance, identical across a 5x spread of `--end`.
///
/// `adaptive-card` reported the last point of its `HoldCurve` for the same reason and with the
/// same consequence, until `HoldCurve` started carrying the flight's reach alongside its
/// samples. This is the assertion that catches a surface going back to the grid for it.
#[test]
fn the_reported_reach_does_not_move_with_end_on_any_surface() {
    for surface in surfaces() {
        let name = surface.name;
        let reach_of = |end: &str| -> String {
            let out = run_argv(&(surface.argv)(Load::Short, "100", end, "csv"));
            assert!(
                out.status.success(),
                "{name}: --end {end} must still print its rows"
            );
            reach_phrase(&String::from_utf8(out.stderr).expect("utf8"))
        };

        let baseline = reach_of("1000");
        for end in ["1200", "1500", "2000", "3000", "5000"] {
            assert_eq!(
                reach_of(end),
                baseline,
                "{name}: the load's reach must not depend on the --end that was asked for"
            );
        }
    }
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

    /// The bridge truncates on the same terms as the CLI, and says so in a field a caller can
    /// read — rather than handing a mobile caller a row built from a different range's numbers,
    /// and rather than failing the whole card over its unreachable tail.
    ///
    /// This test deliberately replaces `card_services_refuse_a_range_the_flight_cannot_reach`,
    /// which asserted the whole-card refusal. Not fabricating the row is still asserted here;
    /// refusing the reachable rows with it was the overcorrection being fixed.
    #[test]
    fn card_services_truncate_a_card_the_flight_outruns() {
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
            let card = service(&unreachable)
                .unwrap_or_else(|e| panic!("card.{name}: the reachable rows must survive: {e}"));

            let last = card.rows.last().expect("at least one reachable row").range;
            assert_eq!(last, 800.0, "card.{name}: rows must stop at the flight's reach");
            assert!(
                !card.rows.iter().any(|row| row.range > 800.0),
                "card.{name}: a row past the flight's reach was fabricated"
            );

            let truncation = card
                .truncation
                .unwrap_or_else(|| panic!("card.{name}: a truncated card must say so"));
            assert_eq!(truncation.requested_end, 2000.0);
            assert_eq!(truncation.last_row, 800.0);
            assert!(
                (truncation.reach - 871.7).abs() < 1.0,
                "card.{name}: reach must be the flight's own terminal distance, got {}",
                truncation.reach
            );
        }
    }

    /// A card whose flight reaches none of its rows still fails, with the range named.
    #[test]
    fn card_services_refuse_a_card_with_no_reachable_row_at_all() {
        let hopeless: CardRequestV1 = serde_json::from_value(serde_json::json!({
            "units": "imperial",
            "muzzle_velocity": 900.0,
            "ballistic_coefficient": 0.1,
            "mass": 40.0,
            "diameter": 0.224,
            "drag_model": "g1",
            "zero_distance": 100.0,
            "start": 1000.0,
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
            match service(&hopeless) {
                Err(CardServiceError::Trajectory(message)) => assert!(
                    message.contains("no trajectory sample at 1000 yd"),
                    "card.{name}: the error must name the range it could not supply, \
                     got: {message}"
                ),
                Err(other) => panic!("card.{name}: unexpected error {other}"),
                Ok(card) => panic!(
                    "card.{name}: returned {} rows for a flight that reaches none of them",
                    card.rows.len()
                ),
            }
        }
    }

    /// An untruncated card carries no truncation at all — the field is additive, so a caller
    /// reading a card that ran to its requested end sees exactly what it saw before.
    #[test]
    fn an_untruncated_card_service_response_has_no_truncation() {
        for (name, service) in [
            ("range_table", range_table_v1 as fn(&CardRequestV1) -> _),
            ("come_ups", come_ups_v1),
            ("wind", wind_card_v1),
        ] {
            let card = service(&request(100.0, 1000.0, 100.0)).expect("card");
            assert!(
                card.truncation.is_none(),
                "card.{name}: a card that ran to its requested end must not claim truncation"
            );
        }
    }
}

// ---------------------------------------------------------------------------------------
// 5. The PRINTED card (MBA-1477).
//
// The truncation notice reached three of the four readers: a stderr warning for a CLI reader,
// a `truncated` block in CLI JSON, `CardResponseV1::truncation` on the bridge. The PDF got
// none of them, on either of its two paths -- a fresh solve through `pdf_card_v1`, and a
// reprint of a card an app had saved -- so a shooter's printed range card stopped at 950 yd of
// a 1200 yd request with nothing on the paper explaining why. That is the one card surface
// most likely to be carried to a range with no screen beside it, and before MBA-1476 the same
// request printed a FABRICATED duplicate row rather than stopping short: the PDF has never
// told the truth about a card that outruns its load, it has only changed which way it lied.
//
// These assertions render a document and read the text back out of it. Asserting from the
// source that the renderer was handed a notice would pass just as happily if the footer never
// drew it.
// ---------------------------------------------------------------------------------------

#[cfg(feature = "pdf")]
mod printed {
    use super::pdf_text;
    use ballistics_engine::card_service::{pdf_card_v1, range_table_v1, CardRequestV1, StoredCardV1};

    /// The ticket's own reproduction: a .22-class load asked for a 100..1200 yd card. The
    /// solved flight ends around 872 yd, so 18 rows print and the rest do not.
    fn truncating_request() -> CardRequestV1 {
        serde_json::from_value(serde_json::json!({
            "units": "imperial",
            "muzzle_velocity": 900.0,
            "ballistic_coefficient": 0.1,
            "mass": 40.0,
            "diameter": 0.224,
            "drag_model": "g1",
            "zero_distance": 100.0,
            "start": 100.0,
            "end": 1200.0,
            "step": 50.0,
            "adjustment_unit": "mil",
        }))
        .expect("card request")
    }

    /// The same load and card, asked for exactly the end it can reach.
    fn whole_request() -> CardRequestV1 {
        let mut request = truncating_request();
        request.end = 850.0;
        request
    }

    /// The one line a truncated 100..1200 card must carry, in the wording every card surface
    /// shares (`CardTruncation::printed_note`).
    const EXPECTED_NOTE: &str = "TRUNCATED at 850 yd: 1200 yd requested, this load reaches only";

    /// SOLVE path: `pdf_card_v1(req, None)` prints the notice in the document itself.
    ///
    /// Verified in the ticket by calling this directly: it returned `Ok` with `row_count` 18,
    /// one page, and a `PdfCardV1` that carried no truncation field at all.
    #[test]
    fn a_truncated_card_says_so_on_the_paper() {
        let card = pdf_card_v1(&truncating_request(), None).expect("the reachable rows must print");

        let truncation = card
            .truncation
            .expect("a truncated PDF must carry the same block the on-screen card does");
        assert_eq!(truncation.requested_end, 1200.0);
        assert_eq!(truncation.last_row, 850.0);
        assert!(
            (truncation.reach - 871.7).abs() < 1.0,
            "the reach must be the flight's own terminal distance, got {}",
            truncation.reach
        );

        assert!(
            pdf_text::contains_phrase(&card.pdf_bytes, "solved truncated card", EXPECTED_NOTE),
            "the printed card must state its truncation; it says:\n{}",
            pdf_text::tokens(&card.pdf_bytes, "solved truncated card").join(" ")
        );
    }

    /// And the structured half agrees with the printed half: the rows it says it stopped at
    /// are the rows it actually printed.
    #[test]
    fn the_printed_truncation_describes_the_rows_that_reached_the_paper() {
        let card = pdf_card_v1(&truncating_request(), None).expect("card");
        let screen = range_table_v1(&truncating_request()).expect("card.range_table");

        assert_eq!(
            card.row_count,
            screen.rows.len(),
            "the paper and the screen print the same rows"
        );
        assert_eq!(
            card.truncation.map(|t| t.last_row),
            screen.truncation.map(|t| t.last_row),
            "the paper and the screen truncate at the same row"
        );
        assert_eq!(
            card.truncation.map(|t| t.last_row),
            screen.rows.last().map(|row| row.range),
            "the notice must name the last row that was actually printed"
        );
    }

    /// Additive: a card that ran to its requested end prints no notice at all, so an
    /// untruncated document is exactly the document it was before.
    #[test]
    fn a_card_that_reaches_its_end_prints_no_notice() {
        let card = pdf_card_v1(&whole_request(), None).expect("card");
        assert!(card.truncation.is_none(), "this card reaches every row it asked for");

        let text = pdf_text::tokens(&card.pdf_bytes, "whole card").join(" ");
        assert!(
            !text.contains("TRUNCATED"),
            "an untruncated card must say nothing about truncation:\n{text}"
        );
    }

    /// REPRINT path: a saved card keeps its notice.
    ///
    /// `StoredCardResponseV1` ignores unknown fields by design (a card saved by a newer engine
    /// must still print on an older one), so until it grew a `truncation` field of its own the
    /// block an app had faithfully stored was dropped on the way in and the reprint came out
    /// silently short -- the screen the card was saved from said it was truncated, and the
    /// paper made from that very document did not.
    #[test]
    fn a_saved_truncated_card_still_says_so_when_it_is_reprinted() {
        let request = truncating_request();
        let screen = range_table_v1(&request).expect("card.range_table");
        let stored_document = serde_json::to_value(&screen).expect("the response an app stores");
        assert!(
            stored_document.get("truncation").is_some(),
            "the saved document must carry the block in the first place: {stored_document}"
        );

        let stored: StoredCardV1 = serde_json::from_value(serde_json::json!({
            "card": stored_document,
            "engine_version": "0.0.0-saved",
        }))
        .expect("a stored card round-trips through its own mirror");

        let card = pdf_card_v1(&request, Some(&stored)).expect("a reprint prints");
        assert_eq!(
            card.truncation.map(|t| t.last_row),
            Some(850.0),
            "a reprint must keep the truncation the card was saved with"
        );
        assert!(
            pdf_text::contains_phrase(&card.pdf_bytes, "reprinted truncated card", EXPECTED_NOTE),
            "the reprinted card must state its truncation; it says:\n{}",
            pdf_text::tokens(&card.pdf_bytes, "reprinted truncated card").join(" ")
        );
    }

    /// A stored card with no block reprints with no notice — which is also every card saved
    /// before the field existed, and the honest rendering of "nothing is known about this".
    #[test]
    fn a_saved_untruncated_card_reprints_without_a_notice() {
        let request = whole_request();
        let screen = range_table_v1(&request).expect("card.range_table");
        let stored: StoredCardV1 = serde_json::from_value(serde_json::json!({
            "card": serde_json::to_value(&screen).expect("stored"),
        }))
        .expect("stored card");

        let card = pdf_card_v1(&request, Some(&stored)).expect("a reprint prints");
        assert!(card.truncation.is_none());
        assert!(
            !pdf_text::tokens(&card.pdf_bytes, "reprinted whole card")
                .join(" ")
                .contains("TRUNCATED"),
            "an untruncated reprint must say nothing about truncation"
        );
    }
}

/// `adaptive-card -o pdf`: the fifth surface's own printed card carries the same line.
///
/// It is the surface MBA-1478 taught to truncate at all, so it is also the surface whose PDF
/// could most easily have been left printing a short card in silence.
#[cfg(feature = "pdf")]
#[test]
fn the_adaptive_cards_printed_form_states_its_truncation() {
    let dir = std::env::temp_dir().join(format!("mba1477-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("scratch dir");
    let path = dir.join("adaptive.pdf");

    let mut args: Vec<&str> = vec!["adaptive-card"];
    args.extend_from_slice(UNREACHABLE_LOAD);
    let path_arg = path.to_string_lossy().to_string();
    args.extend_from_slice(&["--start", "100", "--end", "1200", "-o", "pdf", "--output-file", &path_arg]);
    let out = run(&args);
    assert!(
        out.status.success(),
        "adaptive-card must print the rows it reaches: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let bytes = std::fs::read(&path).expect("the card was written");
    let text = pdf_text::tokens(&bytes, "adaptive card").join(" ");
    assert!(
        text.contains("TRUNCATED at 871 yd:") && text.contains("1200 yd requested"),
        "the printed adaptive card must state its truncation; it says:\n{text}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}
