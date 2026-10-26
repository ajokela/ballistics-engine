//! 0.49.0: the native CLI flies a typed temperature or pressure as typed, on every subcommand,
//! in both unit systems, absolute or QNH, on the shot day and the zero day.
//!
//! Before 0.49.0 a pressure within 0.5 hPa of 1013.25 hPa or a temperature within 0.1 °C of
//! 15 °C at an altitude above 1 m was read as "not given" and replaced with ICAO standard air,
//! because the CLI resolved every `Option` flag to its default and then guessed presence back
//! from the value. Presence now comes from the flag (or the CSV cell, or the saved profile's
//! marker) itself. An omitted value still flies ICAO air at the altitude, bit for bit; the
//! station-sentinel golden matrix checks that byte for byte, and these tests pin the
//! owner-visible numbers and the special paths.
//!
//! Hermetic: anything touching saved profiles points `HOME` at a private temp dir.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU32, Ordering};

use serde_json::Value;

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

const METRIC_BUL: &[&str] = &[
    "--bc",
    "0.475",
    "--drag-model",
    "g1",
    "--velocity",
    "823",
    "-m",
    "10.9",
    "-d",
    "7.82",
];
const METRIC_LOAD: &[&str] = &["-v", "823", "-b", "0.475", "-m", "10.9", "-d", "7.82"];
const IMP_BUL: &[&str] = &[
    "--bc",
    "0.475",
    "--drag-model",
    "g1",
    "--velocity",
    "2700",
    "-m",
    "168",
    "-d",
    "0.308",
];
const IMP_LOAD: &[&str] = &["-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308"];

fn private_home() -> PathBuf {
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "bx-station-presence-cli-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn run_in(home: Option<&Path>, args: &[&str]) -> Output {
    let mut cmd = Command::new(BIN);
    if let Some(home) = home {
        cmd.env("HOME", home);
    }
    cmd.args(args).output().expect("spawn ballistics")
}

fn ok(home: Option<&Path>, args: &[&str]) -> Output {
    let out = run_in(home, args);
    assert!(
        out.status.success(),
        "{args:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out
}

fn json_in(home: Option<&Path>, args: &[&str]) -> Value {
    let mut full: Vec<&str> = args.to_vec();
    full.extend(["-o", "json"]);
    let out = ok(home, &full);
    serde_json::from_slice(&out.stdout)
        .unwrap_or_else(|e| panic!("{args:?}: {e}: {}", String::from_utf8_lossy(&out.stdout)))
}

fn json(args: &[&str]) -> Value {
    json_in(None, args)
}

fn stdout_in(home: Option<&Path>, args: &[&str]) -> String {
    String::from_utf8(ok(home, args).stdout).unwrap()
}

fn stderr_of(args: &[&str]) -> String {
    String::from_utf8(ok(None, args).stderr).unwrap()
}

fn cat<'a>(parts: &[&[&'a str]]) -> Vec<&'a str> {
    parts.iter().flat_map(|p| p.iter().copied()).collect()
}

fn assert_bits(actual: f64, expected: f64, what: &str) {
    assert_eq!(
        actual.to_bits(),
        expected.to_bits(),
        "{what}: {actual:?}, expected {expected:?}"
    );
}

/// `trajectory --max-range 600 --auto-zero 100` impact velocity.
fn trajectory_impact(units: &str, extra: &[&str]) -> f64 {
    let load = if units == "metric" {
        METRIC_LOAD
    } else {
        IMP_LOAD
    };
    let args = cat(&[
        &["-u", units, "trajectory"],
        load,
        &["--max-range", "600", "--auto-zero", "100"],
        extra,
    ]);
    json(&args)["impact_velocity"].as_f64().unwrap()
}

/// `zero --target-distance 300` angle in MOA.
fn zero_moa(units: &str, extra: &[&str]) -> f64 {
    let bul = if units == "metric" {
        METRIC_BUL
    } else {
        IMP_BUL
    };
    let args = cat(&[
        &["-u", units, "zero"],
        bul,
        &["--target-distance", "300"],
        extra,
    ]);
    json(&args)["zero_angle_moa"].as_f64().unwrap()
}

/// `come-ups --zero-distance 100 --end 600`, the 600 row's (velocity, drop).
fn come_ups_600(units: &str, extra: &[&str]) -> (f64, f64) {
    let bul = if units == "metric" {
        METRIC_BUL
    } else {
        IMP_BUL
    };
    let args = cat(&[
        &["-u", units, "come-ups"],
        bul,
        &["--zero-distance", "100", "--end", "600"],
        extra,
    ]);
    let doc = json(&args);
    let row = doc["data"]
        .as_array()
        .unwrap()
        .iter()
        .find(|r| r["range"] == serde_json::json!(600.0))
        .expect("600 row")
        .clone();
    (
        row["velocity"].as_f64().unwrap(),
        row["drop"].as_f64().unwrap(),
    )
}

const A: &[&str] = &[
    "--altitude",
    "150",
    "--temperature",
    "20",
    "--pressure",
    "1031",
    "--pressure-type",
    "qnh",
];
const B: &[&str] = &[
    "--altitude",
    "150",
    "--temperature",
    "20",
    "--pressure",
    "1012.9",
];
const C: &[&str] = &[
    "--altitude",
    "5000",
    "--temperature",
    "70",
    "--pressure",
    "29.92",
];
const D: &[&str] = &[
    "--altitude",
    "5000",
    "--pressure",
    "25.0",
    "--temperature",
    "59",
];
const E: &[&str] = &[
    "--altitude",
    "1500",
    "--pressure",
    "1020",
    "--pressure-type",
    "qnh",
    "--temperature",
    "15",
];

/// The station-sentinel spec's canonical repros (section 3), from its oracle with the
/// sentinels off. "Was" is what ab4eb0c printed (ICAO air at the altitude).
#[test]
fn the_canonical_repros_fly_the_typed_air() {
    // A: QNH 1031 hPa at 150 m reduces to 1012.856 hPa. come-ups was 496.40 m/s, 4.976 mil.
    let (v, d) = come_ups_600("metric", A);
    assert_bits(v, 491.61153010153987, "A come-ups velocity");
    assert_bits(d, 5.020390667980392, "A come-ups drop");
    // trajectory and zero already flew a declared QNH; unchanged.
    assert_bits(
        trajectory_impact("metric", A),
        554.3835755312953,
        "A trajectory",
    );
    assert_bits(zero_moa("metric", A), 9.381729277876076, "A zero");

    // B: 1012.9 hPa typed at 150 m. Was 557.8999196740976 / 9.354189955142585 / 496.40.
    assert_bits(
        trajectory_impact("metric", B),
        554.374816239558,
        "B trajectory",
    );
    assert_bits(zero_moa("metric", B), 9.381729277876076, "B zero");
    let (v, d) = come_ups_600("metric", B);
    assert_bits(v, 491.59958797242194, "B come-ups velocity");
    assert_bits(d, 5.0205012760885595, "B come-ups drop");

    // C: 29.92 inHg typed at 5000 ft. Was 1934.971466102828 fps (the ICAO 24.90 inHg).
    assert_bits(
        trajectory_impact("imperial", C),
        1817.3657972174826,
        "C trajectory",
    );
    assert_bits(zero_moa("imperial", C), 8.568663559077757, "C zero");
    let (v, d) = come_ups_600("imperial", C);
    assert_bits(v, 1694.2464172321636, "C come-ups velocity");
    assert_bits(d, 4.312616421347037, "C come-ups drop");

    // D: 59 F typed at 5000 ft on a 25.0 inHg day. Was 1903.376442150263.
    assert_bits(
        trajectory_impact("imperial", D),
        1921.4805042716182,
        "D trajectory",
    );
    assert_bits(zero_moa("imperial", D), 8.382445281546529, "D zero");
    let (v, d) = come_ups_600("imperial", D);
    assert_bits(v, 1826.6761166696906, "D come-ups velocity");
    assert_bits(d, 4.03341406058797, "D come-ups drop");

    // E: 15 C typed at 1500 m on a QNH 1020 day. Was 580.1206903298374.
    assert_bits(
        trajectory_impact("metric", E),
        585.5496949458625,
        "E trajectory",
    );
    assert_bits(zero_moa("metric", E), 9.156169110725575, "E zero");
    let (v, d) = come_ups_600("metric", E);
    assert_bits(v, 533.6128170428726, "E come-ups velocity");
    assert_bits(d, 4.6699881905728855, "E come-ups drop");

    // M1 / M2: mixed presence. The omitted half still flies ICAO.
    let (v, d) = come_ups_600("metric", &["--altitude", "1500", "--temperature", "15"]);
    assert_bits(v, 535.3431242847164, "M1 velocity");
    assert_bits(d, 4.6565449685633675, "M1 drop");
    let (v, d) = come_ups_600("metric", &["--altitude", "150", "--pressure", "1012.9"]);
    assert_bits(v, 486.60940134345, "M2 velocity");
    assert_bits(d, 5.063875736668549, "M2 drop");
}

/// Every subcommand the probe matrix covers: a typed in-band reading differs from the same
/// command with it omitted. `stability` has no sentinel (its omitted default is sea-level air
/// itself), so only its case A and B, whose typed values differ from that default, move.
#[test]
fn every_subcommand_flies_a_typed_in_band_reading() {
    fn commands(units: &'static str) -> Vec<(&'static str, Vec<String>)> {
        let m = units == "metric";
        let (bul, load) = if m {
            (METRIC_BUL, METRIC_LOAD)
        } else {
            (IMP_BUL, IMP_LOAD)
        };
        let s = |v: &[&str]| v.iter().map(|x| x.to_string()).collect::<Vec<_>>();
        let (twist, length) = if m { ("254", "31.5") } else { ("10", "1.24") };
        let (mass, dia, vel) = if m {
            ("10.9", "7.82", "823")
        } else {
            ("168", "0.308", "2700")
        };
        let data = if m {
            "300,0.85;500,2.7"
        } else {
            "300,34;500,107"
        };
        let load_a = if m {
            "A:g1:0.475:10.9:823:7.82"
        } else {
            "A:g1:0.475:168:2700:0.308"
        };
        let load_b = if m {
            "B:g1:0.500:11.66:838:7.82"
        } else {
            "B:g1:0.500:180:2750:0.308"
        };
        vec![
            (
                "mpbr",
                s(&cat(&[
                    &["mpbr"],
                    bul,
                    &["--vital-zone", if m { "20" } else { "8" }],
                ])),
            ),
            (
                "come-ups",
                s(&cat(&[
                    &["come-ups"],
                    bul,
                    &["--zero-distance", "100", "--end", "600"],
                ])),
            ),
            ("lead", s(&cat(&[&["lead"], bul, &["--target-speed", "4"]]))),
            (
                "wind-card",
                s(&cat(&[&["wind-card"], bul, &["--zero-distance", "100"]])),
            ),
            (
                "range-table",
                s(&cat(&[&["range-table"], bul, &["--zero-distance", "100"]])),
            ),
            (
                "stability",
                s(&[
                    "stability",
                    "--twist-rate",
                    twist,
                    "--length",
                    length,
                    "--mass",
                    mass,
                    "--diameter",
                    dia,
                    "--velocity",
                    vel,
                ]),
            ),
            (
                "true-velocity",
                s(&[
                    "true-velocity",
                    "--measured-drop",
                    "8",
                    "--range",
                    "500",
                    "--bc",
                    "0.475",
                    "--mass",
                    mass,
                    "--diameter",
                    dia,
                    "--offline",
                ]),
            ),
            (
                "true-velocity-observed",
                s(&[
                    "true-velocity",
                    "--measured-drop",
                    "3.6",
                    "--range",
                    "500",
                    "--observed",
                    "700:6.6",
                    "--bc",
                    "0.475",
                    "--mass",
                    mass,
                    "--diameter",
                    dia,
                    "--offline",
                ]),
            ),
            (
                "true-velocity-sigma",
                s(&[
                    "true-velocity",
                    "--measured-drop",
                    "3.6",
                    "--range",
                    "500",
                    "--observed",
                    "700:6.6",
                    "--observation-sigma",
                    "0.1",
                    "--bc",
                    "0.475",
                    "--mass",
                    mass,
                    "--diameter",
                    dia,
                    "--offline",
                ]),
            ),
            (
                "true-velocity-chrono",
                s(&[
                    "true-velocity",
                    "--measured-drop",
                    "8",
                    "--range",
                    "500",
                    "--bc",
                    "0.475",
                    "--mass",
                    mass,
                    "--diameter",
                    dia,
                    "--chrono-velocity",
                    vel,
                    "--chrono-distance",
                    "5",
                    "--offline",
                ]),
            ),
            (
                "true-wind",
                s(&cat(&[
                    &["true-wind", "--miss", "500:15"],
                    bul,
                    &["--twist-rate", twist],
                ])),
            ),
            (
                "plan-truing",
                s(&[
                    "plan-truing",
                    "--measurement-resolution",
                    "0.25",
                    "--candidate-ranges",
                    "300,500,700",
                    "--velocity",
                    vel,
                    "--bc",
                    "0.475",
                    "--mass",
                    mass,
                    "--diameter",
                    dia,
                ]),
            ),
            (
                "compare",
                s(&[
                    "compare",
                    "--zero-distance",
                    "100",
                    "--load",
                    load_a,
                    "--load",
                    load_b,
                ]),
            ),
            (
                "estimate-bc",
                s(&[
                    "estimate-bc",
                    "--velocity",
                    vel,
                    "--mass",
                    mass,
                    "--diameter",
                    dia,
                    "--data",
                    data,
                ]),
            ),
            (
                "mark-to-range",
                s(&cat(&[
                    &["mark-to-range", "--mark", "2.0", "--mark", "5.0"],
                    load,
                    &["--zero-distance", "100"],
                ])),
            ),
            (
                "bdc-match",
                s(&cat(&[
                    &[
                        "bdc-match",
                        "--mark-range",
                        "2:300",
                        "--mark-range",
                        "4:500",
                        "--mark-range",
                        "6:700",
                    ],
                    load,
                    &["--reference-mag", "12", "--zero-distance", "100"],
                ])),
            ),
            (
                "optimal-zero",
                s(&cat(&[
                    &[
                        "optimal-zero",
                        "--target",
                        "200:10",
                        "--target",
                        "400:12",
                        "--target",
                        "600:18",
                    ],
                    load,
                ])),
            ),
            (
                "adaptive-card",
                s(&cat(&[
                    &["adaptive-card"],
                    load,
                    &["--zero-distance", "100", "--start", "100", "--end", "600"],
                ])),
            ),
        ]
    }
    let cases: [(&str, &'static str, &[&str], &[&str]); 5] = [
        (
            "A",
            "metric",
            &["--altitude", "150", "--temperature", "20"],
            &["--pressure", "1031", "--pressure-type", "qnh"],
        ),
        (
            "B",
            "metric",
            &["--altitude", "150", "--temperature", "20"],
            &["--pressure", "1012.9"],
        ),
        (
            "C",
            "imperial",
            &["--altitude", "5000", "--temperature", "70"],
            &["--pressure", "29.92"],
        ),
        (
            "D",
            "imperial",
            &["--altitude", "5000", "--pressure", "25.0"],
            &["--temperature", "59"],
        ),
        (
            "E",
            "metric",
            &[
                "--altitude",
                "1500",
                "--pressure",
                "1020",
                "--pressure-type",
                "qnh",
            ],
            &["--temperature", "15"],
        ),
    ];
    std::thread::scope(|scope| {
        let mut jobs = Vec::new();
        for (case, units, base, supplied) in cases {
            for (name, args) in commands(units) {
                jobs.push(scope.spawn(move || {
                    let mut omitted: Vec<String> = vec!["-u".into(), units.into()];
                    omitted.extend(args.iter().cloned());
                    omitted.extend(base.iter().map(|x| x.to_string()));
                    let mut typed = omitted.clone();
                    typed.extend(supplied.iter().map(|x| x.to_string()));
                    fn refs(v: &[String]) -> Vec<&str> {
                        v.iter().map(String::as_str).collect()
                    }
                    let a = stdout_in(None, &refs(&omitted));
                    let b = stdout_in(None, &refs(&typed));
                    let expect_change = name != "stability" || matches!(case, "A" | "B");
                    // compare echoes the pressure it was handed; strip that line.
                    let strip = |s: &str| {
                        s.lines()
                            .filter(|l| !l.to_lowercase().contains("pressure"))
                            .collect::<Vec<_>>()
                            .join("\n")
                    };
                    (case, name, expect_change, strip(&a) != strip(&b))
                }));
            }
        }
        for job in jobs {
            let (case, name, expect_change, changed) = job.join().unwrap();
            assert_eq!(
                changed, expect_change,
                "case {case} {name}: changed={changed}"
            );
        }
    });
}

const TR_MZ: &[&str] = &[
    "-u",
    "metric",
    "trajectory",
    "-v",
    "823",
    "-b",
    "0.475",
    "-m",
    "10.9",
    "-d",
    "7.82",
    "--max-range",
    "600",
    "--auto-zero",
    "500",
    "--temperature",
    "20",
];
const TR_IZ: &[&str] = &[
    "-u",
    "imperial",
    "trajectory",
    "-v",
    "2700",
    "-b",
    "0.475",
    "-m",
    "168",
    "-d",
    "0.308",
    "--max-range",
    "600",
    "--auto-zero",
    "500",
    "--temperature",
    "70",
];
const TR_M: &[&str] = &[
    "-u",
    "metric",
    "trajectory",
    "-v",
    "823",
    "-b",
    "0.475",
    "-m",
    "10.9",
    "-d",
    "7.82",
    "--max-range",
    "600",
    "--auto-zero",
    "100",
];

fn max_height(args: &[&str]) -> f64 {
    json(args)["max_height"].as_f64().unwrap()
}

#[test]
fn the_zero_day_flies_its_own_typed_air_and_inherits_presence() {
    // S1: --zero-pressure 1012.9 at --zero-altitude 150. Was 2.212887919946592 (ICAO).
    assert_bits(
        max_height(&cat(&[
            TR_MZ,
            &["--zero-altitude", "150", "--zero-pressure", "1012.9"],
        ])),
        2.2201995067396814,
        "S1",
    );
    // The zero day inherits a typed shot-day 1013.25 hPa WITH its presence.
    assert_bits(
        max_height(&cat(&[
            TR_MZ,
            &["--pressure", "1013.25", "--zero-altitude", "150"],
        ])),
        2.2204084092114758,
        "S1 inherited",
    );
    let omitted = max_height(&cat(&[TR_MZ, &["--zero-altitude", "150"]]));
    assert_ne!(omitted.to_bits(), 2.2201995067396814_f64.to_bits());
    // Imperial: 29.92 inHg and 59 F on the zero day.
    let base = cat(&[TR_IZ, &["--zero-altitude", "5000"]]);
    assert_ne!(
        max_height(&cat(&[&base, &["--zero-pressure", "29.92"]])).to_bits(),
        max_height(&base).to_bits()
    );
    let base = cat(&[TR_IZ, &["--zero-altitude", "5000", "--zero-pressure", "25"]]);
    assert_ne!(
        max_height(&cat(&[&base, &["--zero-temperature", "59"]])).to_bits(),
        max_height(&base).to_bits()
    );
    // S1b: the temperature sentinel no longer survives the QNH zero-day path.
    let base = cat(&[
        TR_MZ,
        &[
            "--zero-altitude",
            "1500",
            "--zero-pressure",
            "1020",
            "--zero-pressure-type",
            "qnh",
        ],
    ]);
    assert_ne!(
        max_height(&cat(&[&base, &["--zero-temperature", "15"]])).to_bits(),
        max_height(&base).to_bits()
    );
}

#[test]
fn a_density_altitude_with_a_typed_temperature_flies_that_temperature() {
    // S2: the zero day re-entered the sentinel with the DA-derived air. Was 1.5528031629639063.
    assert_bits(
        max_height(&cat(&[
            TR_M,
            &["--density-altitude", "1500", "--temperature", "15"],
        ])),
        1.5526787637373636,
        "S2",
    );
    // S2b: the zero day's own DA. At a 100 m zero the search quantizes nearby airs to one
    // angle, so a typed 15 C now lands with its 15.2 C neighbour (and the ISA value), where
    // ab4eb0c's swap put it at 1.5525016629998145.
    let base = cat(&[TR_M, &["--zero-density-altitude", "1500"]]);
    let typed = max_height(&cat(&[&base, &["--zero-temperature", "15"]]));
    assert_bits(typed, 1.5523778009403533, "S2b");
    assert_bits(
        typed,
        max_height(&cat(&[&base, &["--zero-temperature", "15.2"]])),
        "S2b neighbour",
    );
    // S3: `zero --density-altitude` read a typed 15 C / 59 F as omitted. Was 9.185019829779707
    // and 8.40867320795938.
    assert_bits(
        zero_moa(
            "metric",
            &["--density-altitude", "1500", "--temperature", "15"],
        ),
        9.194199604024202,
        "S3 metric",
    );
    assert_bits(
        zero_moa(
            "imperial",
            &["--density-altitude", "5000", "--temperature", "59"],
        ),
        8.416541585883234,
        "S3 imperial",
    );
    // S8: an explicit --pressure-type absolute changes nothing; the typed value flies.
    let typed = zero_moa(
        "metric",
        &[
            "--altitude",
            "150",
            "--temperature",
            "20",
            "--pressure",
            "1012.9",
        ],
    );
    assert_bits(
        zero_moa(
            "metric",
            &[
                "--altitude",
                "150",
                "--temperature",
                "20",
                "--pressure",
                "1012.9",
                "--pressure-type",
                "absolute",
            ],
        ),
        typed,
        "S8",
    );
    assert_ne!(
        typed.to_bits(),
        zero_moa("metric", &["--altitude", "150", "--temperature", "20"]).to_bits()
    );
}

fn write_csv(home: &Path, name: &str, body: &str) -> String {
    let path = home.join(name);
    std::fs::write(&path, body).unwrap();
    path.to_string_lossy().into_owned()
}

const TR_I100: &[&str] = &[
    "-u",
    "imperial",
    "trajectory",
    "--max-range",
    "600",
    "--auto-zero",
    "100",
];

#[test]
fn location_and_profile_csv_cells_are_readings() {
    let home = private_home();
    let loc = write_csv(
        &home,
        "loc.csv",
        "NAME,ALTITUDE,TEMPERATURE,PRESSURE\nSUP,150,20,1012.9\nOMIT,150,20,\n\
         TEMP59,5000,59,25.0\nTEMPOMIT,5000,,25.0\nP2992,5000,70,29.92\nPOMIT,5000,70,\n",
    );
    // S4: a CSV cell flies exactly as the same value typed as a flag.
    let site = |units: &str, row: &str| {
        let args = cat(&[
            &["-u", units, "trajectory"],
            if units == "metric" {
                METRIC_LOAD
            } else {
                IMP_LOAD
            },
            &[
                "--max-range",
                "600",
                "--auto-zero",
                "100",
                "--location",
                &loc,
                "--site",
                row,
            ],
        ]);
        json(&args)["impact_velocity"].as_f64().unwrap()
    };
    let flags = |units: &str, extra: &[&str]| trajectory_impact(units, extra);
    assert_bits(
        site("metric", "SUP"),
        flags(
            "metric",
            &[
                "--altitude",
                "150",
                "--temperature",
                "20",
                "--pressure",
                "1012.9",
            ],
        ),
        "S4 1012.9",
    );
    assert_ne!(
        site("metric", "SUP").to_bits(),
        site("metric", "OMIT").to_bits()
    );
    assert_bits(
        site("imperial", "TEMP59"),
        flags(
            "imperial",
            &[
                "--altitude",
                "5000",
                "--temperature",
                "59",
                "--pressure",
                "25.0",
            ],
        ),
        "S4 59 F",
    );
    assert_ne!(
        site("imperial", "TEMP59").to_bits(),
        site("imperial", "TEMPOMIT").to_bits()
    );
    assert_bits(
        site("imperial", "P2992"),
        flags(
            "imperial",
            &[
                "--altitude",
                "5000",
                "--temperature",
                "70",
                "--pressure",
                "29.92",
            ],
        ),
        "S4 29.92",
    );
    assert_ne!(
        site("imperial", "P2992").to_bits(),
        site("imperial", "POMIT").to_bits()
    );

    // S6: a --profile CSV ZERO_TEMP cell. Was 2003.7066512450615.
    let prof = write_csv(
        &home,
        "prof.csv",
        "NAME,VELOCITY,BC,MASS,DIAMETER,ZERO_TEMP,ZERO_ALT\nT59,2700,0.475,168,0.308,59,5000\n",
    );
    let doc = json(&[
        "-u",
        "imperial",
        "trajectory",
        "--profile",
        &prof,
        "--profile-row",
        "T59",
        "-m",
        "168",
        "-d",
        "0.308",
        "--max-range",
        "600",
        "--pressure",
        "25",
    ]);
    assert_bits(
        doc["impact_velocity"].as_f64().unwrap(),
        2020.232528009309,
        "S6",
    );
}

/// E9: under a density altitude, a CSV, profile-CSV or marked saved-profile temperature is the
/// explicit DA temperature, exactly like `--temperature`. An empty cell and an unmarked legacy
/// profile are not.
#[test]
fn density_altitude_takes_a_supplied_temperature_from_any_source() {
    let home = private_home();
    let loc = write_csv(
        &home,
        "loc_da.csv",
        "NAME,DA,TEMPERATURE\nDA59,5000,59\nDAOMIT,5000,\nDA40,5000,40\n",
    );
    let flag_form = |extra: &[&str]| {
        stdout_in(
            None,
            &cat(&[
                TR_I100,
                IMP_LOAD,
                &["--density-altitude", "5000"],
                extra,
                &["-o", "json"],
            ]),
        )
    };
    let csv_form = |row: &str, extra: &[&str]| {
        stdout_in(
            None,
            &cat(&[
                TR_I100,
                IMP_LOAD,
                &["--location", &loc, "--site", row],
                extra,
                &["-o", "json"],
            ]),
        )
    };
    assert_eq!(csv_form("DA59", &[]), flag_form(&["--temperature", "59"]));
    assert_eq!(csv_form("DA40", &[]), flag_form(&["--temperature", "40"]));
    assert_eq!(csv_form("DAOMIT", &[]), flag_form(&[]));
    // Zero-day case 3: the inherited DA temperature keeps its presence next to a --zero-pressure.
    assert_eq!(
        csv_form("DA59", &["--zero-pressure", "25"]),
        flag_form(&["--temperature", "59", "--zero-pressure", "25"])
    );

    // A --profile CSV ZERO_TEMP cell.
    let prof = write_csv(&home, "prof.csv", "NAME,VELOCITY,BC,MASS,DIAMETER,ZERO_TEMP\nT59,2700,0.475,168,0.308,59\nTOMIT,2700,0.475,168,0.308,\n");
    let prof_form = |row: &str| {
        stdout_in(
            None,
            &cat(&[
                TR_I100,
                &[
                    "--profile",
                    &prof,
                    "--profile-row",
                    row,
                    "-m",
                    "168",
                    "-d",
                    "0.308",
                    "--density-altitude",
                    "5000",
                    "-o",
                    "json",
                ],
            ]),
        )
    };
    assert_eq!(prof_form("T59"), flag_form(&["--temperature", "59"]));
    assert_eq!(prof_form("TOMIT"), flag_form(&[]));

    // A marked saved profile, and an unmarked (legacy) one.
    let save = |name: &str, extra: &[&str]| {
        ok(
            Some(&home),
            &cat(&[
                &[
                    "-u",
                    "imperial",
                    "profile",
                    "save",
                    name,
                    "--velocity",
                    "2700",
                    "--bc",
                    "0.475",
                    "--drag-model",
                    "g1",
                    "--mass",
                    "168",
                    "--diameter",
                    "0.308",
                    "--density-altitude",
                    "5000",
                ],
                extra,
            ]),
        );
    };
    save("marked", &["--temperature", "59"]);
    let saved = |name: &str| {
        stdout_in(
            Some(&home),
            &cat(&[TR_I100, &["--saved-profile", name, "-o", "json"]]),
        )
    };
    assert_eq!(saved("marked"), flag_form(&["--temperature", "59"]));
    save("legacy", &["--temperature", "59"]);
    let path = home.join(".ballistics/profiles/legacy.json");
    let text = std::fs::read_to_string(&path).unwrap();
    std::fs::write(
        &path,
        text.replace(",\n  \"temperature_supplied\": true", ""),
    )
    .unwrap();
    assert!(!std::fs::read_to_string(&path)
        .unwrap()
        .contains("_supplied"));
    assert_eq!(saved("legacy"), flag_form(&[]));
}

#[test]
fn saved_profiles_fly_a_marked_reading_and_keep_the_legacy_rule_unmarked() {
    let home = private_home();
    let save = |name: &str, extra: &[&str]| {
        ok(
            Some(&home),
            &cat(&[
                &[
                    "-u",
                    "imperial",
                    "profile",
                    "save",
                    name,
                    "--velocity",
                    "2700",
                    "--bc",
                    "0.475",
                    "--drag-model",
                    "g1",
                    "--mass",
                    "168",
                    "--diameter",
                    "0.308",
                    "--altitude",
                    "5000",
                    "--temperature",
                    "70",
                ],
                extra,
            ]),
        );
    };
    let impact = |name: &str| {
        json_in(Some(&home), &cat(&[TR_I100, &["--saved-profile", name]]))["impact_velocity"]
            .as_f64()
            .unwrap()
    };
    save("typed", &["--pressure", "29.92"]);
    let text = std::fs::read_to_string(home.join(".ballistics/profiles/typed.json")).unwrap();
    assert!(text.contains("\"pressure_supplied\": true"), "{text}");
    // The C repro, through a saved profile.
    assert_bits(
        impact("typed"),
        1817.3657972174826,
        "marked 29.92 at 5000 ft",
    );
    // A hand-written legacy file: the same numbers, no marker.
    std::fs::write(
        home.join(".ballistics/profiles/legacy.json"),
        text.replace("\"typed\"", "\"legacy\"")
            .replace(",\n  \"pressure_supplied\": true", ""),
    )
    .unwrap();
    assert_bits(
        impact("legacy"),
        1934.971466102828,
        "unmarked 29.92 at 5000 ft (ICAO)",
    );
    // E7: a typed --altitude 0 beats the profile's 5000 ft.
    let at_zero = json_in(
        Some(&home),
        &cat(&[TR_I100, &["--saved-profile", "typed", "--altitude", "0"]]),
    );
    let flags = trajectory_impact(
        "imperial",
        &[
            "--altitude",
            "0",
            "--temperature",
            "70",
            "--pressure",
            "29.92",
        ],
    );
    assert_bits(
        at_zero["impact_velocity"].as_f64().unwrap(),
        flags,
        "--altitude 0 over the profile",
    );
}

/// E1 / E2: `plan-truing --profile` applies the profile's stored QNH mode, at the profile's
/// altitude. E10: it honours a stored density altitude, as `trajectory --saved-profile` does.
#[test]
fn plan_truing_honours_the_profiles_stored_mode_and_density_altitude() {
    let home = private_home();
    let save = |name: &str, extra: &[&str]| {
        ok(
            Some(&home),
            &cat(&[
                &[
                    "-u",
                    "metric",
                    "profile",
                    "save",
                    name,
                    "--velocity",
                    "823",
                    "--bc",
                    "0.475",
                    "--drag-model",
                    "g1",
                    "--mass",
                    "10.9",
                    "--diameter",
                    "7.82",
                    "--temperature",
                    "20",
                ],
                extra,
            ]),
        );
    };
    let plan = |extra: &[&str]| {
        stdout_in(
            Some(&home),
            &cat(&[
                &[
                    "-u",
                    "metric",
                    "plan-truing",
                    "--measurement-resolution",
                    "0.25",
                    "--candidate-ranges",
                    "300,500,700",
                ],
                extra,
                &["-o", "json"],
            ]),
        )
    };
    let load = [
        "--velocity",
        "823",
        "--bc",
        "0.475",
        "--mass",
        "10.9",
        "--diameter",
        "7.82",
    ];
    save(
        "qnh",
        &[
            "--altitude",
            "150",
            "--pressure",
            "1031",
            "--pressure-type",
            "qnh",
        ],
    );
    // The CLI pre-reduces a display-unit QNH by scaling it with the station/QNH ratio
    // (`apply_pressure_mode`), so the absolute equivalent is computed the same way.
    let station = ballistics_engine::atmosphere::reduce_qnh_to_station_pressure(1031.0, 150.0);
    let reduced = format!("{:?}", 1031.0 * (station / 1031.0));
    assert_eq!(
        plan(&["--profile", "qnh"]),
        plan(&cat(&[
            &load,
            &[
                "--temperature",
                "20",
                "--altitude",
                "150",
                "--pressure",
                &reduced
            ]
        ]))
    );

    save("da", &["--density-altitude", "1500"]);
    let (a, t, p) =
        ballistics_engine::atmosphere::resolve_atmosphere_for_density_altitude(1500.0, Some(20.0));
    let (a, t, p) = (format!("{a:?}"), format!("{t:?}"), format!("{p:?}"));
    assert_eq!(
        plan(&["--profile", "da"]),
        plan(&cat(&[
            &load,
            &["--altitude", &a, "--temperature", &t, "--pressure", &p]
        ]))
    );
    save("plain", &[]);
    assert_ne!(plan(&["--profile", "da"]), plan(&["--profile", "plain"]));
}

/// E3 / E10: `dsf --saved-profile` applies a stored QNH and a stored density altitude.
#[test]
fn dsf_honours_the_profiles_stored_mode_and_density_altitude() {
    let home = private_home();
    let save = |name: &str, extra: &[&str]| {
        ok(
            Some(&home),
            &cat(&[
                &[
                    "-u",
                    "metric",
                    "profile",
                    "save",
                    name,
                    "--velocity",
                    "350",
                    "--bc",
                    "0.475",
                    "--drag-model",
                    "g1",
                    "--mass",
                    "10.9",
                    "--diameter",
                    "7.82",
                    "--zero-distance",
                    "100",
                    "--temperature",
                    "20",
                ],
                extra,
            ]),
        );
    };
    let dsf = |name: &str| {
        let out = run_in(
            Some(&home),
            &[
                "-u",
                "metric",
                "dsf",
                "--saved-profile",
                name,
                "--range",
                "200",
                "--observed-drop",
                "60.0moa",
            ],
        );
        let text = format!(
            "{}{}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
        text.replace(name, "<P>")
    };
    let reduced = format!(
        "{:?}",
        ballistics_engine::atmosphere::reduce_qnh_to_station_pressure(1031.0, 150.0)
    );
    save(
        "qnh",
        &[
            "--altitude",
            "150",
            "--pressure",
            "1031",
            "--pressure-type",
            "qnh",
        ],
    );
    save("reduced", &["--altitude", "150", "--pressure", &reduced]);
    save("raw", &["--altitude", "150", "--pressure", "1031"]);
    let q = dsf("qnh");
    assert_eq!(q, dsf("reduced"));
    assert_ne!(q, dsf("raw"));

    let (a, t, p) =
        ballistics_engine::atmosphere::resolve_atmosphere_for_density_altitude(1500.0, Some(20.0));
    // An explicit temperature is honoured exactly, so the helper's --temperature 20 is `t`.
    assert_eq!(t.to_bits(), 20.0_f64.to_bits());
    save("da", &["--density-altitude", "1500"]);
    save(
        "triple",
        &[
            "--altitude",
            &format!("{a:?}"),
            "--pressure",
            &format!("{p:?}"),
        ],
    );
    let tri = std::fs::read_to_string(home.join(".ballistics/profiles/triple.json")).unwrap();
    assert!(
        tri.contains("\"temperature_supplied\": true")
            && tri.contains("\"pressure_supplied\": true")
    );
    assert_eq!(dsf("da"), dsf("triple"));
}

#[test]
fn a_typed_implausible_station_pressure_now_warns() {
    let warning = "times the ICAO standard";
    // 29.92 inHg typed at 5000 ft is 1.20 times standard there.
    let typed = cat(&[
        &["-u", "imperial", "trajectory"],
        IMP_LOAD,
        &["--altitude", "5000", "--pressure", "29.92"],
    ]);
    assert!(stderr_of(&typed).contains(warning), "{}", stderr_of(&typed));
    let zero = cat(&[
        &["-u", "imperial", "zero"],
        IMP_BUL,
        &[
            "--target-distance",
            "300",
            "--altitude",
            "5000",
            "--pressure",
            "29.92",
        ],
    ]);
    assert!(stderr_of(&zero).contains(warning));
    let zero_day = cat(&[
        TR_IZ,
        &["--zero-altitude", "5000", "--zero-pressure", "29.92"],
    ]);
    assert!(
        stderr_of(&zero_day).contains("zero-day station pressure"),
        "{}",
        stderr_of(&zero_day)
    );
    // An omitted pressure is never warned.
    let omitted = cat(&[
        &["-u", "imperial", "trajectory"],
        IMP_LOAD,
        &["--altitude", "5000"],
    ]);
    assert!(!stderr_of(&omitted).contains(warning));
    let omitted_zero = cat(&[TR_IZ, &["--zero-altitude", "5000"]]);
    assert!(!stderr_of(&omitted_zero).contains(warning));
}
