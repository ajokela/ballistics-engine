//! 0.49.0: the native CLI flies a typed temperature or pressure as typed, on every subcommand,
//! in both unit systems, absolute or QNH, on the shot day and the zero day.
//!
//! Before 0.49.0 a pressure within 0.5 hPa of 1013.25 hPa or a temperature within 0.1 °C of
//! 15 °C at an altitude above 1 m was read as "not given" and replaced with ICAO standard air,
//! because the CLI resolved every `Option` flag to its default and then guessed presence back
//! from the value. Presence now comes from the flag (or the CSV cell, or the saved profile's
//! marker) itself. An omitted value still flies ICAO air at the altitude, bit for bit (checked
//! against ab4eb0c), and these tests pin the user-visible numbers and the special paths.
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

/// The reported cases, captured from ab4eb0c with the two sentinels compiled out (a typed value
/// flies as typed there). "Was" is what ab4eb0c printed (ICAO air at the altitude).
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

/// Write the `hold-corridor` scenarios file and a generated `reticle hold` reticle into `home`,
/// returning their paths.
fn hold_fixtures(home: &Path) -> (String, String) {
    let scenarios = write_csv(
        home,
        "scenarios.json",
        r#"{"version":1,"nominal":"low","scenarios":[{"name":"low","segments":["4:90:1000"]},{"name":"high","segments":["14:90:1000"]},{"name":"switch","segments":["10:90:400","8:270:1000"]}]}"#,
    );
    let reticle = ok(
        None,
        &[
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
        ],
    );
    let reticle = write_csv(
        home,
        "reticle.json",
        &String::from_utf8(reticle.stdout).unwrap(),
    );
    (scenarios, reticle)
}

/// Every subcommand the probe matrix covers: a typed in-band reading differs from the same
/// command with it omitted. `stability` has no sentinel (its omitted default is sea-level air
/// itself), so only its case A and B, whose typed values differ from that default, move.
#[test]
fn every_subcommand_flies_a_typed_in_band_reading() {
    let home = private_home();
    let (scenarios, reticle) = hold_fixtures(&home);
    let commands = |units: &'static str| -> Vec<(&'static str, Vec<String>)> {
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
            (
                "hold-corridor",
                s(&cat(&[
                    &[
                        "hold-corridor",
                        "--scenarios",
                        &scenarios,
                        "--ranges",
                        "200,400,600",
                        "--target",
                        "rect:12x18",
                    ],
                    load,
                    &["--zero-distance", "100"],
                ])),
            ),
            (
                "reticle-hold",
                s(&cat(&[
                    &[
                        "reticle",
                        "hold",
                        "--reticle-json",
                        &reticle,
                        "--range",
                        "500",
                    ],
                    load,
                    &["--zero-distance", "100"],
                ])),
            ),
            (
                "zero-from-angle",
                s(&cat(&[&["zero"], bul, &["--from-angle", "0.1"]])),
            ),
        ]
    };
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

/// One number from a command's JSON, by a `/`-separated path; `-1` is an array's last element.
fn json_at(doc: &Value, path: &str) -> f64 {
    let mut v = doc;
    for seg in path.split('/') {
        v = match v {
            Value::Array(items) if seg == "-1" => items.last().unwrap(),
            Value::Array(items) => &items[seg.parse::<usize>().unwrap()],
            _ => &v[seg],
        };
    }
    v.as_f64()
        .unwrap_or_else(|| panic!("{path} is not a number: {doc}"))
}

/// The zero flies the typed air too, on every subcommand that solves one. At a 100 m zero the
/// zero-angle search lands nearby airs on one angle, so a run whose zero fell back to the old
/// rule (ICAO air) while its flight used the typed air looked right; at 300 m it does not. Each
/// value was captured from ab4eb0c with the two sentinels compiled out (a typed value flies as
/// typed there), and equals this build's output; "was" is what ab4eb0c printed. Case A is QNH
/// 1031 hPa at 150 m (1012.86 hPa station), B is 1012.9 hPa typed at 150 m, both 20 °C; D is
/// 59 °F typed at 5000 ft on a 25.0 inHg day. `zero --from-angle` already flew a declared QNH,
/// so its case A is unchanged.
#[test]
fn every_subcommand_flies_the_typed_air_at_its_zero_too() {
    let home = private_home();
    let (scenarios, reticle) = hold_fixtures(&home);
    let command = |name: &str, metric: bool| -> Vec<String> {
        let (bul, load) = if metric {
            (METRIC_BUL, METRIC_LOAD)
        } else {
            (IMP_BUL, IMP_LOAD)
        };
        let (mass, dia, vel) = if metric {
            ("10.9", "7.82", "823")
        } else {
            ("168", "0.308", "2700")
        };
        let twist = if metric { "254" } else { "10" };
        let (load_a, load_b) = if metric {
            ("A:g1:0.475:10.9:823:7.82", "B:g1:0.500:11.66:838:7.82")
        } else {
            ("A:g1:0.475:168:2700:0.308", "B:g1:0.500:180:2750:0.308")
        };
        let data = if metric {
            "300,0.85;500,2.7"
        } else {
            "300,34;500,107"
        };
        let z: &[&str] = &["--zero-distance", "300"];
        let tv: &[&str] = &[
            "--bc",
            "0.475",
            "--mass",
            mass,
            "--diameter",
            dia,
            "--offline",
            "--zero-distance",
            "300",
        ];
        let args: Vec<&str> = match name {
            "zero-from-angle" => cat(&[&["zero"], bul, &["--from-angle", "0.1"]]),
            "mpbr" => cat(&[
                &["mpbr"],
                bul,
                &["--vital-zone", if metric { "20" } else { "8" }],
            ]),
            "come-ups" => cat(&[&["come-ups"], bul, z, &["--end", "600"]]),
            "lead" => cat(&[&["lead"], bul, &["--target-speed", "4"]]),
            "wind-card" => cat(&[&["wind-card"], bul, z]),
            "range-table" => cat(&[&["range-table"], bul, z]),
            "true-velocity" => cat(&[
                &["true-velocity", "--measured-drop", "8", "--range", "500"],
                tv,
            ]),
            "true-velocity-observed" => cat(&[
                &[
                    "true-velocity",
                    "--measured-drop",
                    "3.6",
                    "--range",
                    "500",
                    "--observed",
                    "700:6.6",
                ],
                tv,
            ]),
            "true-velocity-chrono" => cat(&[
                &[
                    "true-velocity",
                    "--measured-drop",
                    "8",
                    "--range",
                    "500",
                    "--chrono-velocity",
                    vel,
                    "--chrono-distance",
                    "5",
                ],
                tv,
            ]),
            "true-wind" => cat(&[
                &["true-wind", "--miss", "500:15"],
                bul,
                &["--twist-rate", twist],
                z,
            ]),
            "plan-truing" => cat(&[
                &[
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
                ],
                z,
            ]),
            "compare" => cat(&[&["compare", "--load", load_a, "--load", load_b], z]),
            "estimate-bc" => vec![
                "estimate-bc",
                "--velocity",
                vel,
                "--mass",
                mass,
                "--diameter",
                dia,
                "--data",
                data,
                "--zero-range",
                "100",
            ],
            "mark-to-range" => cat(&[
                &["mark-to-range", "--mark", "2.0", "--mark", "5.0"],
                load,
                z,
            ]),
            "bdc-match" => cat(&[
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
                &["--reference-mag", "12"],
                z,
            ]),
            "optimal-zero" => cat(&[
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
            ]),
            "adaptive-card" => cat(&[
                &["adaptive-card"],
                load,
                z,
                &["--start", "100", "--end", "600"],
            ]),
            "hold-corridor" => cat(&[
                &[
                    "hold-corridor",
                    "--scenarios",
                    &scenarios,
                    "--ranges",
                    "200,400,600",
                    "--target",
                    "rect:12x18",
                ],
                load,
                z,
            ]),
            "reticle-hold" => cat(&[
                &[
                    "reticle",
                    "hold",
                    "--reticle-json",
                    &reticle,
                    "--range",
                    "500",
                ],
                load,
                z,
            ]),
            other => panic!("no command {other}"),
        };
        args.into_iter().map(String::from).collect()
    };
    let air = |case: &str| -> (bool, &'static [&'static str]) {
        match case {
            "A" => (
                true,
                &[
                    "--altitude",
                    "150",
                    "--temperature",
                    "20",
                    "--pressure",
                    "1031",
                    "--pressure-type",
                    "qnh",
                ],
            ),
            "B" => (
                true,
                &[
                    "--altitude",
                    "150",
                    "--temperature",
                    "20",
                    "--pressure",
                    "1012.9",
                ],
            ),
            "D" => (
                false,
                &[
                    "--altitude",
                    "5000",
                    "--pressure",
                    "25.0",
                    "--temperature",
                    "59",
                ],
            ),
            other => panic!("no case {other}"),
        }
    };
    // (case, command, value, now, was at ab4eb0c)
    #[rustfmt::skip]
    let pins: &[(&str, &str, &str, f64, f64)] = &[
        ("A", "zero-from-angle", "zero_range", 184.44065946144926, 184.44065946144926),
        ("A", "mpbr", "mpbr", 277.57525170751234, 277.7648122247587),
        ("A", "come-ups", "data/-1/drop", 3.554752524040378, 3.5188299573131094),
        ("A", "lead", "rows/-1/lead", 3.7827576282715856, 3.763753693530174),
        ("A", "wind-card", "data/-1/wind_10", -7.110933393129025, -6.94583022877229),
        ("A", "range-table", "data/-1/drop", 20206.33848586181, 19840.774089725764),
        ("A", "true-velocity", "calculated_drop_mil", 7.483927282826406, 7.485529736545842),
        ("A", "true-velocity-observed", "fitted_muzzle_velocity", 722.3621602163888, 718.6512144304594),
        ("A", "true-velocity-observed", "mv_window_start_m", 565.6397822977608, 569.2620823181188),
        ("A", "true-velocity-observed", "mv_window_end_m", 628.488646997512, 632.5134247979098),
        ("A", "true-velocity-chrono", "adjustment_percent", -44.654661921496206, -44.65095500018336),
        ("A", "true-wind", "effective_crosswind", 2.12881831805145, 2.172083532303539),
        ("A", "plan-truing", "information/condition_number", 645.8374395583315, 490.42237919090815),
        ("A", "compare", "compare/rows/-1/loads/0/drop", 1134.6135560103025, 1124.9489467973397),
        ("A", "estimate-bc", "variants/0/estimated_bc", 0.25, 0.246),
        ("A", "mark-to-range", "marks/-1/range", 693.8915418242184, 697.4353937913447),
        ("A", "bdc-match", "fitted_magnification", 17.046803176810492, 17.234551565451333),
        ("A", "optimal-zero", "optimal_zero", 428.4605521015392, 427.99036890521256),
        ("A", "adaptive-card", "rows/-1/drop_adj", 3.554752524040378, 3.5188299573131094),
        ("A", "hold-corridor", "rows/-1/minimax_elevation_mil", 3.554867259390395, 3.518940460585484),
        ("A", "reticle-hold", "hold/down_mil", 2.2044838121635126, 2.1855361168682452),
        ("B", "zero-from-angle", "zero_range", 184.43974635007845, 184.80632783196094),
        ("B", "mpbr", "mpbr", 277.5736752363296, 277.7648122247587),
        ("B", "come-ups", "data/-1/drop", 3.554863128953085, 3.5188299573131094),
        ("B", "lead", "rows/-1/lead", 3.7828053900044294, 3.763753693530174),
        ("B", "wind-card", "data/-1/wind_10", -7.111347180475116, -6.94583022877229),
        ("B", "range-table", "data/-1/drop", 20207.28174930798, 19840.774089725764),
        ("B", "true-velocity", "calculated_drop_mil", 7.483991955554914, 7.485529736545842),
        ("B", "true-velocity-observed", "fitted_muzzle_velocity", 722.4107728583531, 718.6512144304594),
        ("B", "true-velocity-observed", "mv_window_start_m", 565.6657738930813, 569.2620823181188),
        ("B", "true-velocity-observed", "mv_window_end_m", 628.5175265478681, 632.5134247979098),
        ("B", "true-velocity-chrono", "adjustment_percent", -44.65467119987547, -44.65095500018336),
        ("B", "true-wind", "effective_crosswind", 2.128711796782978, 2.172083532303539),
        ("B", "plan-truing", "information/condition_number", 645.139317487319, 490.42237919090815),
        ("B", "compare", "compare/rows/-1/loads/0/drop", 1134.6473510599828, 1124.9489467973397),
        ("B", "estimate-bc", "variants/0/estimated_bc", 0.25, 0.246),
        ("B", "mark-to-range", "marks/-1/range", 693.8815080866484, 697.4353937913447),
        ("B", "bdc-match", "fitted_magnification", 17.046232235566688, 17.234551565451333),
        ("B", "optimal-zero", "optimal_zero", 428.4547651761321, 427.99036890521256),
        ("B", "adaptive-card", "rows/-1/drop_adj", 3.554863128953085, 3.5188299573131094),
        ("B", "hold-corridor", "rows/-1/minimax_elevation_mil", 3.5549778707511086, 3.518940460585484),
        ("B", "reticle-hold", "hold/down_mil", 2.204551402110831, 2.1855361168682452),
        ("D", "zero-from-angle", "zero_range", 204.55089485001628, 203.9627191953953),
        ("D", "mpbr", "mpbr", 311.3027544989547, 309.4064871371572),
        ("D", "come-ups", "data/-1/drop", 2.8401683241504876, 2.8771593427158986),
        ("D", "lead", "rows/-1/lead", 1.5864756044145347, 1.5966150517244653),
        ("D", "wind-card", "data/-1/wind_10", -2.2555711195619996, -2.3492521052917232),
        ("D", "range-table", "data/-1/drop", 521.2425486288677, 535.0865231210498),
        ("D", "true-velocity", "calculated_drop_mil", 6.186676348196818, 6.299968420836825),
        ("D", "true-velocity-observed", "fitted_muzzle_velocity", 1882.6119671554238, 2140.5495855985573),
        ("D", "true-velocity-observed", "mv_window_start_m", 803.3096554604081, 549.1094520822987),
        ("D", "true-velocity-observed", "mv_window_end_m", 892.5662838448978, 610.1216134247763),
        ("D", "true-velocity-chrono", "adjustment_percent", -44.4931569064095, -44.49484406197996),
        ("D", "true-wind", "effective_crosswind", 7.279270215377554, 7.058519377470756),
        ("D", "plan-truing", "information/condition_number", 258.5000158802868, 278.5062568096683),
        ("D", "compare", "compare/rows/-1/loads/0/drop", 32.530981086706674, 32.909841432326665),
        ("D", "estimate-bc", "variants/0/estimated_bc", 0.159, 0.164),
        ("D", "mark-to-range", "marks/-1/range", 777.9575523622334, 771.9128933586179),
        ("D", "bdc-match", "fitted_magnification", 21.547484910326702, 21.249038846017886),
        ("D", "optimal-zero", "optimal_zero", 423.81437898352885, 424.39776778488186),
        ("D", "adaptive-card", "rows/-1/drop_adj", 2.8401683241504876, 2.8771593427158986),
        ("D", "hold-corridor", "rows/-1/minimax_elevation_mil", 2.8402563969250787, 2.8771198979743846),
        ("D", "reticle-hold", "hold/down_mil", 1.789278277384567, 1.8095631018860499),
    ];
    std::thread::scope(|scope| {
        let jobs: Vec<_> = pins
            .iter()
            .map(|&(case, name, path, now, was)| {
                let (metric, air) = air(case);
                let mut args: Vec<String> = vec![
                    "-u".into(),
                    if metric { "metric" } else { "imperial" }.into(),
                ];
                args.extend(command(name, metric));
                args.extend(air.iter().map(|x| x.to_string()));
                scope.spawn(move || {
                    let refs: Vec<&str> = args.iter().map(String::as_str).collect();
                    let got = json_at(&json(&refs), path);
                    (case, name, path, now, was, got)
                })
            })
            .collect();
        for job in jobs {
            let (case, name, path, now, was, got) = job.join().unwrap();
            assert_bits(
                got,
                now,
                &format!("case {case} {name} {path} (was {was:?})"),
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

/// Under a density altitude, a CSV, profile-CSV or marked saved-profile temperature is the
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
    // A typed --altitude 0 beats the profile's 5000 ft.
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

/// `plan-truing --profile` applies the profile's stored QNH mode, at the profile's altitude, and
/// honours a stored density altitude, as `trajectory --saved-profile` does.
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

/// `dsf --saved-profile` applies a stored QNH and a stored density altitude.
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

/// Copy a saved profile to `to` with its presence markers removed, which is how a profile
/// written before 0.49.0 reads.
fn unmark_profile(home: &Path, from: &str, to: &str) {
    let dir = home.join(".ballistics/profiles");
    let text = std::fs::read_to_string(dir.join(format!("{from}.json"))).unwrap();
    let legacy: String = text
        .lines()
        .filter(|line| !line.contains("_supplied"))
        .map(|line| {
            if line.trim_start().starts_with("\"name\"") {
                format!("  \"name\": \"{to}\",\n")
            } else {
                format!("{line}\n")
            }
        })
        .collect();
    assert!(!legacy.contains("_supplied"), "{legacy}");
    std::fs::write(dir.join(format!("{to}.json")), legacy).unwrap();
}

/// A profile written before 0.49.0 carries no presence markers, but a stored QNH in it is still
/// a QNH: `trajectory --saved-profile` flies its reduced station pressure whatever the value,
/// and so must `plan-truing --profile` and `dsf --saved-profile`. QNH 1031 hPa at 150 m (and
/// 30.4455 inHg at 492.126 ft) reduces to about 1012.86 hPa, inside the legacy band; it used to
/// fly ICAO air there, the same air as an omitted pressure.
#[test]
fn an_unmarked_profiles_stored_qnh_flies_its_reduced_pressure() {
    struct Case {
        units: &'static str,
        load: [&'static str; 8],
        air: [&'static str; 4],
        qnh: &'static str,
        dsf_velocity: &'static str,
        dsf_range: &'static str,
        dsf_drop: &'static str,
    }
    let cases = [
        Case {
            units: "metric",
            load: [
                "--velocity",
                "823",
                "--bc",
                "0.475",
                "--mass",
                "10.9",
                "--diameter",
                "7.82",
            ],
            air: ["--altitude", "150", "--temperature", "20"],
            qnh: "1031",
            dsf_velocity: "350",
            dsf_range: "200",
            dsf_drop: "15moa",
        },
        Case {
            units: "imperial",
            load: [
                "--velocity",
                "2700",
                "--bc",
                "0.475",
                "--mass",
                "168",
                "--diameter",
                "0.308",
            ],
            air: ["--altitude", "492.126", "--temperature", "68"],
            qnh: "30.4455",
            dsf_velocity: "1150",
            dsf_range: "250",
            dsf_drop: "20moa",
        },
    ];
    for case in &cases {
        // A HOME per unit system: dsf writes its point back into the profile it solves.
        let home = private_home();
        let u = case.units;
        let save = |name: &str, velocity: &str, extra: &[&str]| {
            let mut load = case.load;
            load[1] = velocity;
            ok(
                Some(&home),
                &cat(&[
                    &["-u", u, "profile", "save", name, "--drag-model", "g1"],
                    &load,
                    &["--zero-distance", "100"],
                    &case.air,
                    extra,
                ]),
            );
        };
        let typed_qnh = ["--pressure", case.qnh, "--pressure-type", "qnh"];
        let plan = |extra: &[&str]| {
            stdout_in(
                Some(&home),
                &cat(&[
                    &[
                        "-u",
                        u,
                        "plan-truing",
                        "--measurement-resolution",
                        "0.25",
                        "--candidate-ranges",
                        "300,500,700",
                        "-o",
                        "json",
                    ],
                    extra,
                ]),
            )
        };

        // plan-truing --profile, stored mode.
        save("marked", case.load[1], &typed_qnh);
        unmark_profile(&home, "marked", "legacy");
        let typed = plan(&cat(&[&case.load, &case.air, &typed_qnh]));
        let omitted = plan(&cat(&[&case.load, &case.air]));
        assert_ne!(typed, omitted, "{u}: the QNH must move the plan");
        assert_eq!(plan(&["--profile", "marked"]), typed, "{u}: marked");
        assert_eq!(plan(&["--profile", "legacy"]), typed, "{u}: unmarked");

        // plan-truing --profile with the mode typed for the profile's own pressure.
        save("raw", case.load[1], &["--pressure", case.qnh]);
        unmark_profile(&home, "raw", "raw-legacy");
        assert_eq!(
            plan(&["--profile", "raw-legacy", "--pressure-type", "qnh"]),
            typed,
            "{u}: unmarked absolute profile flown as a typed QNH"
        );

        // dsf --saved-profile, stored mode. A subsonic load, where the air shows in the DSF.
        let dsf = |name: &str| {
            let out = ok(
                Some(&home),
                &[
                    "-u",
                    u,
                    "dsf",
                    "--saved-profile",
                    name,
                    "--range",
                    case.dsf_range,
                    "--observed-drop",
                    case.dsf_drop,
                ],
            );
            format!(
                "{}{}",
                String::from_utf8_lossy(&out.stdout),
                String::from_utf8_lossy(&out.stderr)
            )
            .replace(name, "<P>")
        };
        save("dsf-marked", case.dsf_velocity, &typed_qnh);
        unmark_profile(&home, "dsf-marked", "dsf-legacy");
        save("dsf-omitted", case.dsf_velocity, &[]);
        let marked = dsf("dsf-marked");
        assert_ne!(marked, dsf("dsf-omitted"), "{u}: the QNH must move the DSF");
        assert_eq!(dsf("dsf-legacy"), marked, "{u}: unmarked dsf");
    }
}

/// A profile saved without `--pressure` stores the sea-level default (29.92 inHg or 1013.25 hPa)
/// and no marker. Flown as a QNH (stored, or `--pressure-type qnh` typed for it), that default is
/// still an omitted pressure: `plan-truing` reduces it and the old rule decides, exactly as with
/// no profile. At 3 m (10 ft) the reduced default lands inside the legacy band and flies the ICAO
/// pressure there. The default is recognised in either unit system, so a profile loaded in the
/// other one is omitted too. Every pinned value was captured from ab4eb0c. Only a stored value
/// other than the default is a reading in an unmarked profile.
#[test]
fn a_profiles_omitted_pressure_flown_as_a_qnh_stays_omitted() {
    let home = private_home();
    let metric_load = [
        "--velocity",
        "823",
        "--bc",
        "0.475",
        "--mass",
        "10.9",
        "--diameter",
        "7.82",
    ];
    let imperial_load = [
        "--velocity",
        "2700",
        "--bc",
        "0.475",
        "--mass",
        "168",
        "--diameter",
        "0.308",
    ];
    let metric_air = ["--altitude", "3", "--temperature", "20"];
    let imperial_air = ["--altitude", "10", "--temperature", "68"];
    let save = |units: &str, name: &str, load: &[&str], air: &[&str], extra: &[&str]| {
        ok(
            Some(&home),
            &cat(&[&["-u", units, "profile", "save", name], load, air, extra]),
        );
    };
    let qnh = ["--pressure-type", "qnh"];
    save("metric", "a-metric", &metric_load, &metric_air, &[]);
    save("metric", "q-metric", &metric_load, &metric_air, &qnh);
    save("imperial", "a-imperial", &imperial_load, &imperial_air, &[]);
    save(
        "imperial",
        "q-imperial",
        &imperial_load,
        &imperial_air,
        &qnh,
    );
    for name in ["a-metric", "q-metric", "a-imperial", "q-imperial"] {
        let file = home.join(format!(".ballistics/profiles/{name}.json"));
        let text = std::fs::read_to_string(file).unwrap();
        assert!(!text.contains("pressure_supplied"), "{name}: {text}");
    }
    let plan = |units: &str, extra: &[&str]| {
        json_in(
            Some(&home),
            &cat(&[
                &[
                    "-u",
                    units,
                    "plan-truing",
                    "--measurement-resolution",
                    "0.25",
                    "--candidate-ranges",
                    "300,500,700",
                ],
                extra,
            ]),
        )["information"]["condition_number"]
            .as_f64()
            .unwrap()
    };
    // (units, run, condition number captured from ab4eb0c)
    let cases: [(&str, Vec<&str>, f64); 12] = [
        (
            "metric",
            cat(&[&metric_load, &metric_air, &qnh]),
            147.22911163966467,
        ),
        (
            "metric",
            cat(&[&["--profile", "a-metric", "--altitude", "3"], &qnh]),
            147.22911163966467,
        ),
        (
            "metric",
            cat(&[&["--profile", "a-metric"], &qnh]),
            147.22911163966467,
        ),
        ("metric", vec!["--profile", "q-metric"], 147.22911163966467),
        (
            "metric",
            cat(&[&["--profile", "a-imperial", "--altitude", "3"], &qnh]),
            146.77803224676077,
        ),
        (
            "metric",
            vec!["--profile", "q-imperial"],
            146.78702305164347,
        ),
        (
            "imperial",
            cat(&[&imperial_load, &imperial_air, &qnh]),
            158.11027132668036,
        ),
        (
            "imperial",
            cat(&[&["--profile", "a-imperial", "--altitude", "10"], &qnh]),
            158.11027132668036,
        ),
        (
            "imperial",
            cat(&[&["--profile", "a-imperial"], &qnh]),
            158.11027132668036,
        ),
        (
            "imperial",
            vec!["--profile", "q-imperial"],
            158.11027132668036,
        ),
        (
            "imperial",
            cat(&[&["--profile", "a-metric", "--altitude", "10"], &qnh]),
            158.11768242204755,
        ),
        (
            "imperial",
            vec!["--profile", "q-metric"],
            158.11766725855364,
        ),
    ];
    for (units, run, was) in &cases {
        assert_bits(plan(units, run), *was, &format!("{units} {run:?}"));
    }

    // A stored value other than the default is a reading, marker or not: QNH 1013.5 hPa reduces
    // to about 1013.14 hPa at 3 m, inside the band, and flies as reduced.
    save(
        "metric",
        "q-typed",
        &metric_load,
        &metric_air,
        &["--pressure", "1013.5", "--pressure-type", "qnh"],
    );
    unmark_profile(&home, "q-typed", "q-typed-legacy");
    let typed = plan(
        "metric",
        &cat(&[
            &metric_load,
            &metric_air,
            &["--pressure", "1013.5", "--pressure-type", "qnh"],
        ]),
    );
    assert_ne!(typed.to_bits(), 147.22911163966467_f64.to_bits());
    assert_bits(plan("metric", &["--profile", "q-typed"]), typed, "marked");
    assert_bits(
        plan("metric", &["--profile", "q-typed-legacy"]),
        typed,
        "unmarked",
    );
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

/// A one-shot HTTP service that records the request body and answers 500, so `true-velocity
/// --offline-fallback` falls back to its local solve.
#[cfg(feature = "online")]
fn failing_service() -> (String, std::thread::JoinHandle<String>) {
    use std::io::{Read, Write};
    let listener = std::net::TcpListener::bind("127.0.0.1:0").unwrap();
    let url = format!("http://{}", listener.local_addr().unwrap());
    let server = std::thread::spawn(move || {
        let (mut stream, _) = listener.accept().unwrap();
        let mut request = Vec::new();
        let mut buffer = [0_u8; 4096];
        let body_start = loop {
            let n = stream.read(&mut buffer).unwrap();
            assert!(n > 0, "connection closed before the headers ended");
            request.extend_from_slice(&buffer[..n]);
            if let Some(i) = request.windows(4).position(|w| w == b"\r\n\r\n") {
                break i + 4;
            }
        };
        let headers = String::from_utf8_lossy(&request[..body_start]).to_lowercase();
        let length: usize = headers
            .lines()
            .find_map(|l| l.strip_prefix("content-length:"))
            .map(|v| v.trim().parse().unwrap())
            .unwrap_or(0);
        while request.len() < body_start + length {
            let n = stream.read(&mut buffer).unwrap();
            assert!(n > 0, "connection closed before the body ended");
            request.extend_from_slice(&buffer[..n]);
        }
        stream
            .write_all(
                b"HTTP/1.1 500 Internal Server Error\r\nContent-Type: application/json\r\n\
                  Content-Length: 2\r\nConnection: close\r\n\r\n{}",
            )
            .unwrap();
        String::from_utf8(request[body_start..body_start + length].to_vec()).unwrap()
    });
    (url, server)
}

/// `true-velocity` online sends a typed reading with its presence marker and an omitted one
/// exactly as before (no marker key). When the service fails, `--offline-fallback` solves
/// locally with the same presence as `--offline`, so a typed in-band reading flies as typed.
#[cfg(feature = "online")]
#[test]
fn true_velocity_online_sends_presence_and_falls_back_with_it() {
    let home = private_home();
    std::fs::create_dir_all(home.join(".ballistics")).unwrap();
    std::fs::write(
        home.join(".ballistics/tos.json"),
        r#"{"accepted": true, "accepted_version": "2026-01-26", "accepted_at": "test", "terms_hash": "test"}"#,
    )
    .unwrap();
    let base: &[&str] = &[
        "true-velocity",
        "--measured-drop",
        "8",
        "--range",
        "500",
        "--bc",
        "0.475",
        "--zero-distance",
        "300",
        "-o",
        "json",
    ];
    // (units, air flags, expected temperature marker, expected pressure marker)
    type Case = (
        &'static str,
        &'static [&'static str],
        Option<bool>,
        Option<bool>,
    );
    let cases: [Case; 3] = [
        (
            "metric",
            &[
                "--mass",
                "10.9",
                "--diameter",
                "7.82",
                "--altitude",
                "150",
                "--temperature",
                "20",
                "--pressure",
                "1012.9",
            ],
            Some(true),
            Some(true),
        ),
        (
            "imperial",
            &[
                "--mass",
                "168",
                "--diameter",
                "0.308",
                "--altitude",
                "5000",
                "--temperature",
                "59",
            ],
            Some(true),
            None,
        ),
        (
            "metric",
            &["--mass", "10.9", "--diameter", "7.82", "--altitude", "150"],
            None,
            None,
        ),
    ];
    for (units, air, temperature_marker, pressure_marker) in cases {
        let (url, server) = failing_service();
        let online = run_in(
            Some(&home),
            &cat(&[
                &["-u", units],
                base,
                air,
                &["--offline-fallback", "--api-url", &url],
            ]),
        );
        let body: Value = serde_json::from_str(&server.join().unwrap()).unwrap();
        assert!(
            online.status.success(),
            "{air:?}: {}",
            String::from_utf8_lossy(&online.stderr)
        );
        assert!(String::from_utf8_lossy(&online.stderr).contains("Falling back"));
        assert_eq!(
            body.get("temperature_supplied").and_then(Value::as_bool),
            temperature_marker,
            "{air:?}: {body}"
        );
        assert_eq!(
            body.get("pressure_supplied").and_then(Value::as_bool),
            pressure_marker,
            "{air:?}: {body}"
        );
        let offline = stdout_in(
            Some(&home),
            &cat(&[&["-u", units], base, air, &["--offline"]]),
        );
        assert_eq!(
            String::from_utf8(online.stdout).unwrap(),
            offline,
            "{air:?}: the fallback must fly what --offline flies"
        );
    }
}

/// `dsf --saved-profile` flies a marked in-band reading as typed, and keeps the legacy rule for
/// the same value in an unmarked (pre-0.49.0) profile, where it is the omitted default.
#[test]
fn dsf_flies_a_marked_in_band_reading_and_keeps_the_legacy_rule_unmarked() {
    let home = private_home();
    let save = |name: &str, air: &[&str]| {
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
                ],
                air,
            ]),
        );
    };
    let dsf = |name: &str| {
        let out = ok(
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
                "15moa",
            ],
        );
        format!(
            "{}{}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        )
        .replace(name, "<P>")
    };
    for (what, air, without) in [
        (
            "1012.9 hPa at 150 m",
            &[
                "--altitude",
                "150",
                "--temperature",
                "20",
                "--pressure",
                "1012.9",
            ][..],
            &["--altitude", "150", "--temperature", "20"][..],
        ),
        (
            "15 C at 1500 m",
            &[
                "--altitude",
                "1500",
                "--pressure",
                "846",
                "--temperature",
                "15",
            ][..],
            &["--altitude", "1500", "--pressure", "846"][..],
        ),
    ] {
        // dsf writes its point back into the profile, so each profile is solved once.
        let (typed, legacy, omitted) = (
            format!("typed-{}", air[1]),
            format!("legacy-{}", air[1]),
            format!("omitted-{}", air[1]),
        );
        save(&typed, air);
        unmark_profile(&home, &typed, &legacy);
        save(&omitted, without);
        let (typed, legacy, omitted) = (dsf(&typed), dsf(&legacy), dsf(&omitted));
        assert_ne!(typed, omitted, "{what}: the reading must fly");
        assert_eq!(legacy, omitted, "{what}: unmarked keeps the legacy rule");
    }
}
