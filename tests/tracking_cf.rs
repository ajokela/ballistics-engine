//! MBA-1358: scope tracking correction factors from a tall-target test.
//!
//! CF semantics pinned here (the published Litz tall-target convention, direction
//! adjudicated by the coordinator — physics wins):
//!   * the stored CF is `actual measured travel / dialed travel` (0.95 = the scope
//!     under-tracks 5%); dial-unit outputs (mil/moa/smoa/iphy/clicks) are DIVIDED by
//!     the axis's CF exactly once at the shared conversion boundary — an
//!     under-tracking scope needs MORE dial — and raw drop/drift inches never scale;
//!   * truing `--observed`/`--measured-drop` dialed values are MULTIPLIED by the CF
//!     (scope-dial units -> true angular: CF=0.95, dialed 10.526… -> the fit sees 10),
//!     and dial-unit report values are shown back in scope units (÷CF);
//!   * CF = 1.0 (or absent) is byte-identical on every surface;
//!   * factors outside the open interval (0.5, 1.5) are hard errors naming the source;
//!   * `tall-target` prints CF = actual measured travel / dialed travel, no solve.

use serde_json::Value;
use std::process::Command;
use std::sync::atomic::{AtomicU32, Ordering};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

fn run(args: &[&str]) -> std::process::Output {
    Command::new(BIN).args(args).output().expect("run binary")
}

fn run_ok(args: &[&str]) -> String {
    let output = run(args);
    assert!(
        output.status.success(),
        "command failed: {:?}\nstderr: {}",
        args,
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).expect("utf8")
}

fn json_ok(args: &[&str]) -> Value {
    serde_json::from_str(&run_ok(args)).expect("json output")
}

const LOAD: &[&str] = &[
    "-b", "0.243", "--drag-model", "g7", "-v", "2700", "-m", "175", "-d", "0.308",
];

// true-velocity takes no -v (it SOLVES the velocity).
const TV_LOAD: &[&str] = &["-b", "0.243", "--drag-model", "g7", "-m", "175", "-d", "0.308"];

// (a) explicit CF = 1.0 is byte-identical on every dial surface.
#[test]
fn cf_of_one_is_byte_identical() {
    let cases: Vec<(Vec<&str>, Vec<&str>)> = vec![
        (
            [&["come-ups"], LOAD, &["--zero-distance", "100", "--start", "200", "--end", "400", "--step", "100"]].concat(),
            vec!["--elevation-cf", "1.0"],
        ),
        (
            [&["range-table"], LOAD, &["--zero-distance", "100", "--start", "200", "--end", "400", "--step", "100"]].concat(),
            vec!["--elevation-cf", "1.0", "--windage-cf", "1.0"],
        ),
        (
            [&["wind-card"], LOAD, &["--zero-distance", "100", "--start", "200", "--end", "300", "--step", "100"]].concat(),
            vec!["--windage-cf", "1.0"],
        ),
        (
            [&["lead"], LOAD, &["--target-speed", "3", "--start", "200", "--end", "300", "--step", "100"]].concat(),
            vec!["--windage-cf", "1.0"],
        ),
        (
            [&["trajectory"], LOAD, &["--max-range", "500", "--auto-zero", "100", "--target-speed", "10", "-o", "csv", "--full"]].concat(),
            vec!["--elevation-cf", "1.0", "--windage-cf", "1.0"],
        ),
        (
            // --offline: the single-observation path would otherwise try the online
            // API under the default `online` feature (network-dependent test).
            [&["true-velocity"], TV_LOAD, &["--range", "600", "--measured-drop", "4.1", "--offline"]].concat(),
            vec!["--elevation-cf", "1.0"],
        ),
        (
            [&["zero"], LOAD, &["--target-distance", "100"]].concat(),
            vec!["--elevation-cf", "1.0"],
        ),
    ];
    for (base, extra) in cases {
        let mut with_cf = base.clone();
        with_cf.extend(extra.iter());
        assert_eq!(
            run_ok(&with_cf),
            run_ok(&base),
            "explicit CF 1.0 must not change a single output byte: {base:?}"
        );
    }
}

// (b) CF = 0.95 (an under-tracking scope) INCREASES every dial surface's values by
// exactly 1/0.95 — dial-unit outputs divide by the CF.
#[test]
fn cf_divides_every_dial_surface_exactly() {
    let close = |scaled: f64, base: f64, what: &str| {
        assert!(
            (scaled - base / 0.95).abs() <= 1e-9 * base.abs().max(1.0),
            "{what}: {scaled} != {base} / 0.95"
        );
        // direction pin: an under-tracking scope needs MORE dial (sign-aware).
        assert!(
            scaled.abs() >= base.abs(),
            "{what}: |{scaled}| must not shrink vs |{base}| under CF 0.95"
        );
    };

    // come-ups: drop and come_up columns are dial values.
    let base_args =
        [&["come-ups"][..], LOAD, &["--zero-distance", "100", "--start", "200", "--end", "400", "--step", "100", "-o", "json"]].concat();
    let mut cf_args = base_args.clone();
    cf_args.extend(["--elevation-cf", "0.95"]);
    let (b, s) = (json_ok(&base_args), json_ok(&cf_args));
    for (rb, rs) in b["data"].as_array().unwrap().iter().zip(s["data"].as_array().unwrap()) {
        close(rs["drop"].as_f64().unwrap(), rb["drop"].as_f64().unwrap(), "come-ups drop");
        close(rs["come_up"].as_f64().unwrap(), rb["come_up"].as_f64().unwrap(), "come-ups come_up");
    }

    // range-table CSV: drop_mil/wind_mil scale; drop_in/wind_in (raw inches) must NOT.
    let base_args =
        [&["range-table"][..], LOAD, &["--zero-distance", "100", "--start", "200", "--end", "400", "--step", "100", "-o", "csv"]].concat();
    let mut cf_args = base_args.clone();
    cf_args.extend(["--elevation-cf", "0.95", "--windage-cf", "0.95"]);
    let (b, s) = (run_ok(&base_args), run_ok(&cf_args));
    for (lb, ls) in b.lines().skip(1).zip(s.lines().skip(1)) {
        let fb: Vec<&str> = lb.split(',').collect();
        let fs: Vec<&str> = ls.split(',').collect();
        // header: range_yd,drop_in,drop_mil,wind_in,wind_mil,...
        assert_eq!(fb[1], fs[1], "raw drop inches must never scale");
        assert_eq!(fb[3], fs[3], "raw wind inches must never scale");
        // drop_mil prints at 3 decimals — compare with print-rounding tolerance.
        let (db, ds): (f64, f64) = (fb[2].parse().unwrap(), fs[2].parse().unwrap());
        assert!(
            (ds - db / 0.95).abs() < 0.002,
            "range-table drop_mil: {ds} vs {db} / 0.95"
        );
        // wind values print at 2 decimals — compare with print tolerance.
        let (wb, ws): (f64, f64) = (fb[4].parse().unwrap(), fs[4].parse().unwrap());
        assert!(
            (ws - wb / 0.95).abs() < 0.02,
            "range-table wind_mil: {ws} vs {wb} / 0.95"
        );
    }

    // wind-card: every wind_N cell is a dial value.
    let base_args =
        [&["wind-card"][..], LOAD, &["--zero-distance", "100", "--start", "200", "--end", "300", "--step", "100", "-o", "json"]].concat();
    let mut cf_args = base_args.clone();
    cf_args.extend(["--windage-cf", "0.95"]);
    let (b, s) = (json_ok(&base_args), json_ok(&cf_args));
    for (rb, rs) in b["data"].as_array().unwrap().iter().zip(s["data"].as_array().unwrap()) {
        for speed in ["wind_5", "wind_10", "wind_15", "wind_20"] {
            close(rs[speed].as_f64().unwrap(), rb[speed].as_f64().unwrap(), "wind-card drift");
        }
    }

    // lead: lead_mil/lead_moa scale (dialed quantity); linear lead must NOT.
    let base_args =
        [&["lead"][..], LOAD, &["--target-speed", "3", "--start", "200", "--end", "300", "--step", "100", "-o", "json"]].concat();
    let mut cf_args = base_args.clone();
    cf_args.extend(["--windage-cf", "0.95"]);
    let (b, s) = (json_ok(&base_args), json_ok(&cf_args));
    for (rb, rs) in b["rows"].as_array().unwrap().iter().zip(s["rows"].as_array().unwrap()) {
        close(rs["lead_mil"].as_f64().unwrap(), rb["lead_mil"].as_f64().unwrap(), "lead_mil");
        close(rs["lead_moa"].as_f64().unwrap(), rb["lead_moa"].as_f64().unwrap(), "lead_moa");
        assert_eq!(
            rs["lead"].as_f64().unwrap().to_bits(),
            rb["lead"].as_f64().unwrap().to_bits(),
            "linear lead distance must never scale"
        );
    }

    // trajectory mover Ring: mover_ring_mil scales, mover_ring_m (raw meters) must NOT.
    let base_args =
        [&["trajectory"][..], LOAD, &["--max-range", "500", "--auto-zero", "100", "--target-speed", "10", "-o", "json", "--full"]].concat();
    let mut cf_args = base_args.clone();
    cf_args.extend(["--windage-cf", "0.95"]);
    let (b, s) = (json_ok(&base_args), json_ok(&cf_args));
    let (pb, ps) = (b["trajectory"].as_array().unwrap(), s["trajectory"].as_array().unwrap());
    let mut checked = 0;
    for (rb, rs) in pb.iter().zip(ps) {
        if let (Some(mb), Some(ms)) = (rb["mover_ring_mil"].as_f64(), rs["mover_ring_mil"].as_f64()) {
            close(ms, mb, "mover_ring_mil");
            assert_eq!(
                rs["mover_ring_m"].as_f64().unwrap().to_bits(),
                rb["mover_ring_m"].as_f64().unwrap().to_bits(),
                "mover_ring_m is raw meters and must never scale"
            );
            checked += 1;
        }
    }
    assert!(checked > 3, "expected several ring points, got {checked}");

    // zero: the MOA/mrad dial outputs scale; the degrees bore-angle echo must NOT.
    let base_args = [&["zero"][..], LOAD, &["--target-distance", "100", "-o", "json"]].concat();
    let mut cf_args = base_args.clone();
    cf_args.extend(["--elevation-cf", "0.95"]);
    let (b, s) = (json_ok(&base_args), json_ok(&cf_args));
    close(
        s["zero_angle_moa"].as_f64().unwrap(),
        b["zero_angle_moa"].as_f64().unwrap(),
        "zero_angle_moa",
    );
    close(
        s["zero_angle_mrad"].as_f64().unwrap(),
        b["zero_angle_mrad"].as_f64().unwrap(),
        "zero_angle_mrad",
    );
    close(
        s["sight_adjustment_moa"].as_f64().unwrap(),
        b["sight_adjustment_moa"].as_f64().unwrap(),
        "sight_adjustment_moa",
    );
    assert_eq!(
        s["zero_angle_degrees"].as_f64().unwrap().to_bits(),
        b["zero_angle_degrees"].as_f64().unwrap().to_bits(),
        "the degrees echo is a bore angle and must never scale"
    );
}

// (b2) clicks quantize the CF-corrected (divided) angle.
#[test]
fn cf_corrects_clicks_through_the_quantizer() {
    let base_args = [
        &["come-ups"][..],
        LOAD,
        &[
            "--zero-distance", "100", "--start", "200", "--end", "600", "--step", "100",
            "-o", "json", "--adjustment-unit", "clicks", "--elevation-click-value", "0.1mil",
        ],
    ]
    .concat();
    let mut cf_args = base_args.clone();
    cf_args.extend(["--elevation-cf", "0.95"]);

    // The uncorrected mil values (same run, mil unit) predict the corrected click
    // counts: clicks = round((mil / 0.95) / 0.1).
    let mil_args = [
        &["come-ups"][..],
        LOAD,
        &["--zero-distance", "100", "--start", "200", "--end", "600", "--step", "100", "-o", "json"],
    ]
    .concat();
    let mil = json_ok(&mil_args);
    let clicks = json_ok(&cf_args);
    for (rm, rc) in mil["data"]
        .as_array()
        .unwrap()
        .iter()
        .zip(clicks["data"].as_array().unwrap())
    {
        let expected = (rm["drop"].as_f64().unwrap() / 0.95 / 0.1).round();
        assert_eq!(
            rc["drop"].as_f64().unwrap(),
            expected,
            "clicks must quantize the CF-scaled angle at {} yd",
            rm["range"]
        );
    }
}

// (c) truing conversion direction: CF=0.95, dialed 10.526… -> the fit sees 10 true.
#[test]
fn truing_multiplies_dialed_observations_by_the_cf() {
    // Single observation: 10.0/0.95 dialed mils with CF 0.95 must fit the same
    // velocity as a 10-mil run without a CF (dialed × CF = true).
    let with_cf = [
        &["true-velocity"][..],
        TV_LOAD,
        &[
            "--range", "600", "--measured-drop", "10.526315789473685",
            "--elevation-cf", "0.95", "--offline", "-o", "json",
        ],
    ]
    .concat();
    let without = [
        &["true-velocity"][..],
        TV_LOAD,
        &["--range", "600", "--measured-drop", "10", "--offline", "-o", "json"],
    ]
    .concat();
    let (cf_doc, base_doc) = (json_ok(&with_cf), json_ok(&without));
    // 10.526…*0.95 can differ from 10.0 by an ulp, and the bisection's bracket
    // collapses at 0.5 fps — 1 fps still pins the direction unambiguously (a fit of
    // 10.53 dialed mils without the conversion lands hundreds of fps away).
    let (v_cf, v_base) = (
        cf_doc["effective_velocity"].as_f64().unwrap(),
        base_doc["effective_velocity"].as_f64().unwrap(),
    );
    assert!(
        (v_cf - v_base).abs() <= 1.0,
        "the fit must consume 10.526… × 0.95 = 10 true mils: {v_cf} vs {v_base}"
    );
    // ...while the echo shows what the shooter dialed (scope units).
    let echoed = cf_doc["measured_drop_mil"].as_f64().unwrap();
    assert!(
        (echoed - 10.526315789473685).abs() < 1e-9,
        "echo must stay the dialed 10.53 mils, got {echoed}"
    );
    // and the model's drop renders in scope units too: base calculated / 0.95.
    let (calc_cf, calc_base) = (
        cf_doc["calculated_drop_mil"].as_f64().unwrap(),
        base_doc["calculated_drop_mil"].as_f64().unwrap(),
    );
    assert!(
        (calc_cf - calc_base / 0.95).abs() < 1e-3,
        "calculated drop must render in scope units: {calc_cf} vs {calc_base} / 0.95"
    );

    // Multi-observation: converting every dialed drop reproduces the un-CF'd fit.
    // 4.1/0.95 = 4.315789…, 7.9/0.95 = 8.315789… dialed with CF 0.95 == (4.1, 7.9) true.
    let with_cf = [
        &["true-velocity"][..],
        TV_LOAD,
        &[
            "--range", "600", "--measured-drop", "4.315789473684211",
            "--observed", "800:8.31578947368421",
            "--elevation-cf", "0.95", "-o", "json",
        ],
    ]
    .concat();
    let without = [
        &["true-velocity"][..],
        TV_LOAD,
        &["--range", "600", "--measured-drop", "4.1", "--observed", "800:7.9", "-o", "json"],
    ]
    .concat();
    let (cf_doc, base_doc) = (json_ok(&with_cf), json_ok(&without));
    let (mv_cf, mv_base) = (
        cf_doc["fitted_muzzle_velocity"].as_f64().unwrap(),
        base_doc["fitted_muzzle_velocity"].as_f64().unwrap(),
    );
    assert!(
        (mv_cf - mv_base).abs() < 1.0,
        "multi-obs fit must consume the converted (true) drops: {mv_cf} vs {mv_base}"
    );
    let (bc_cf, bc_base) = (
        cf_doc["fitted_bc"].as_f64().unwrap(),
        base_doc["fitted_bc"].as_f64().unwrap(),
    );
    assert!((bc_cf - bc_base).abs() < 1e-3, "{bc_cf} vs {bc_base}");
    // The report echoes the dialed (scope-unit) observations back.
    let echoed = cf_doc["observations"][0]["observed_drop_mil"].as_f64().unwrap();
    assert!(
        (echoed - 4.315789473684211).abs() < 1e-9,
        "multi-obs echo must show the entered dialed value, got {echoed}"
    );

    // Linear 'in' drops are tape measurements — the CF must NOT divide them.
    let with_cf = [
        &["true-velocity"][..],
        TV_LOAD,
        &[
            "--range", "600", "--measured-drop", "95.0", "--observed", "800:225.0",
            "--drop-unit", "in", "--elevation-cf", "0.95", "-o", "json",
        ],
    ]
    .concat();
    let without = [
        &["true-velocity"][..],
        TV_LOAD,
        &[
            "--range", "600", "--measured-drop", "95.0", "--observed", "800:225.0",
            "--drop-unit", "in", "-o", "json",
        ],
    ]
    .concat();
    assert_eq!(
        run_ok(&with_cf),
        run_ok(&without),
        "linear inch observations must be untouched by the CF"
    );
}

// (d) validation: 0.4 and 1.6 are hard errors, on flags and on profile load.
#[test]
fn validation_rejects_out_of_band_factors() {
    for bad in ["0.4", "1.6"] {
        let args = [
            &["come-ups"][..],
            LOAD,
            &["--zero-distance", "100", "--start", "200", "--end", "300", "--step", "100", "--elevation-cf", bad],
        ]
        .concat();
        let output = run(&args);
        assert!(!output.status.success(), "--elevation-cf {bad} must be rejected");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains("--elevation-cf") && stderr.contains("0.5") && stderr.contains("1.5"),
            "error must name the flag and the band: {stderr}"
        );
    }

    // Profile field validated ON LOAD, naming the field.
    static N: AtomicU32 = AtomicU32::new(0);
    let home = std::env::temp_dir().join(format!(
        "bx-cf-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&home).unwrap();
    let save = Command::new(BIN)
        .env("HOME", &home)
        .args([
            "profile", "save", "cf-test", "--velocity", "2700", "--bc", "0.475",
            "--drag-model", "g1", "--mass", "168", "--diameter", "0.308",
        ])
        .output()
        .expect("profile save");
    assert!(save.status.success());
    let path = home.join(".ballistics").join("profiles").join("cf-test.json");
    let mut profile: Value = serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
    profile["elevation_cf"] = Value::from(1.6);
    std::fs::write(&path, serde_json::to_string_pretty(&profile).unwrap()).unwrap();

    let output = Command::new(BIN)
        .env("HOME", &home)
        .args([
            "come-ups", "--profile", "cf-test", "--zero-distance", "100",
            "--start", "200", "--end", "300", "--step", "100",
        ])
        .output()
        .expect("come-ups run");
    assert!(!output.status.success(), "an out-of-band stored CF must fail on load");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("elevation_cf"),
        "load error must name the profile field: {stderr}"
    );

    // A stored valid CF applies, and the CLI flag overrides it.
    profile["elevation_cf"] = Value::from(1.1);
    std::fs::write(&path, serde_json::to_string_pretty(&profile).unwrap()).unwrap();
    let stored = Command::new(BIN)
        .env("HOME", &home)
        .args([
            "come-ups", "--profile", "cf-test", "--zero-distance", "100",
            "--start", "200", "--end", "300", "--step", "100", "-o", "json",
        ])
        .output()
        .expect("come-ups stored cf");
    assert!(stored.status.success());
    let flagged = Command::new(BIN)
        .env("HOME", &home)
        .args([
            "come-ups", "--profile", "cf-test", "--zero-distance", "100",
            "--start", "200", "--end", "300", "--step", "100", "-o", "json",
            "--elevation-cf", "1.0",
        ])
        .output()
        .expect("come-ups flag override");
    assert!(flagged.status.success());
    let stored: Value = serde_json::from_slice(&stored.stdout).unwrap();
    let flagged: Value = serde_json::from_slice(&flagged.stdout).unwrap();
    let (ds, df) = (
        stored["data"][0]["drop"].as_f64().unwrap(),
        flagged["data"][0]["drop"].as_f64().unwrap(),
    );
    assert!(
        (ds - df / 1.1).abs() < 1e-9,
        "stored CF 1.1 must apply as a divisor ({ds}) and the 1.0 flag must override it ({df})"
    );
}

// (e) tall-target helper: CF = actual measured travel / dialed travel, pure arithmetic.
#[test]
fn tall_target_computes_actual_over_dialed() {
    // 36.0 in at 100 yd is exactly 10 mil (1 mil = 3.6 in / 100 yd): CF = 1.0000.
    let out = run_ok(&["tall-target", "--dialed", "10", "--measured", "36", "--range", "100"]);
    assert!(out.contains("Correction factor (actual / dialed): 1.0000"), "{out}");

    // 37.8 in = 10.5 mil: CF = 1.0500, and the printout states the division rule.
    let out = run_ok(&["tall-target", "--dialed", "10", "--measured", "37.8", "--range", "100"]);
    assert!(out.contains("Correction factor (actual / dialed): 1.0500"), "{out}");
    assert!(
        out.contains("dial solutions are divided by it"),
        "the helper must state the application direction: {out}"
    );

    // 34.2 in = 9.5 mil: the under-tracking worked example, CF = 0.9500.
    let out = run_ok(&["tall-target", "--dialed", "10", "--measured", "34.2", "--range", "100"]);
    assert!(out.contains("Correction factor (actual / dialed): 0.9500"), "{out}");

    // MOA uses the locked printed-table factor (3438): 20.94 in at 100 yd ≈ 20 MOA.
    let out = run_ok(&[
        "tall-target", "--dialed", "20", "--measured", "20.94", "--range", "100", "--unit", "moa",
    ]);
    assert!(out.contains("Correction factor (actual / dialed): 0.9999"), "{out}");

    // Out-of-band results warn instead of pretending to be usable.
    let out = run_ok(&["tall-target", "--dialed", "10", "--measured", "72", "--range", "100"]);
    assert!(out.contains("outside the accepted (0.5, 1.5) band"), "{out}");

    // No solve: clicks is not an angular dial unit here.
    let output = run(&["tall-target", "--dialed", "10", "--measured", "36", "--range", "100", "--unit", "clicks"]);
    assert!(!output.status.success());
}
